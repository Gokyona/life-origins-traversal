"""Stage-3d: spatial version of evolve.py (independent-model cross-check).

Identical to evolve.py (levels L1..L4, heritable quality q, neutral markers)
except for space:
  * the CELLS slots form a 2D torus grid (16 x 8 for 128 slots; 8 x 4 at
    fidelity 0); each slot has 8 neighbours
  * the second daughter goes to a neighbouring slot: a dead/empty neighbour
    first, otherwise a random non-dividing neighbour. When two dividing cells
    pick the same slot in the same step, a random one wins and the other's
    second daughter is lost (its parent still divides, so its contents halve)
  * leak: every molecule moves to a random neighbouring slot with probability
    lk per step (all molecules leaving one cell in a step go to the same,
    randomly chosen neighbour). lk is a new condition (linear, 0 .. LK_MAX)
Selection is therefore local: a fast lineage spreads only into its neighbourhood.

(evolve.py description follows)

Stage-3c: Darwinian evolution (level L4) on the stage-2b protocell model.

Extends protocell.py / knockouts.py. Every cell carries two heritable numbers:
  q   catalytic quality of its machinery (starts at 1). It speeds up the whole
      cell metabolism with diminishing returns, h(q) = 2q / (1 + q): h(1) = 1,
      h(0) = 0, h -> 2 as q grows, so quality can at most double the rate.
      h multiplies gene, parasite and membrane synthesis alike (parasites use
      the same machinery; membrane keeps pace, so homeostasis is not broken
      merely by faster growth). q is clipped to [0, Q_MAX].
  q'_1..q'_M  M_NEUTRAL neutral markers: inherited, mutated and clipped exactly
      like q but with no effect.
At division each daughter gets its parent's q + N(0, eta) and, independently,
each q'_m + N(0, eta). eta (heritable variation per division) is a new condition.
Everything else is the protocell model (optionally with a knockout, as in
knockouts.py: "no_selection", "equal_split", ...).

Detector L4 (computable per simulation): L3 holds AND the mean q of functional
cells at the end beats the neutral markers of the same cells:
  z = (mean q - mean over m of mean q'_m) / sd over m of mean q'_m  > L4_Z
and mean q - mean q' > L4_DELTA. The markers share the genealogy, so without
selection mean q behaves like one more marker and z follows a t distribution
with M_NEUTRAL - 1 degrees of freedom (scaled by sqrt(1 + 1/M)); the markers'
spread measures how much genetic drift this particular population had, which
a single paired marker cannot (one neutral marker gave 26-39% false positives
at a fixed threshold, validate_evolve.py v1).

Conditions, each u in [0,1]:
  0..7  S, sigma, mu, K, d, food, km, dm   as in protocell.py
  8     eta   heritable variation per division      0 .. 0.2
  9     lk    leak to neighbours per molecule/step  0 .. 0.01
  10,11 z1, z2 distractors
"""
import math
import os
import sys

import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stage2"))
import protocell as pc  # noqa: E402
from protocell import (CELLS, STEPS, DT, G_MIN, W_MIN, BURST, HOMEO_LO, HOMEO_HI, LYSIS_MAX,  # noqa: F401
                       K_MAX, FIDELITY, _growth)

DIM = 12
NAMES = ["S", "sigma", "mu", "K", "d", "food", "km", "dm", "eta", "lk", "z1", "z2"]
MONOTONE = [0, -1, -1, 0, -1, 1, 0, -1, 0, 0, 0, 0]
LEVEL = 4
KNOCKOUT = None
Q_MAX = 10.0
ETA_MAX = 0.2
LK_MAX = 0.01
GRID = {128: (16, 8), 32: (8, 4)}
L4_DELTA = 0.15
L4_Z = 3.0
M_NEUTRAL = 8
NEUTRAL_Q = False      # validation switch: q has no effect either (pure double-neutral control)


def decode(u):
    S, sigma, mu, K, d, food, km, dm = pc.decode(torch.cat([u[:, :8], u[:, 10:12]], 1))
    eta = ETA_MAX * u[:, 8]
    lk = LK_MAX * u[:, 9]
    return S, sigma, mu, K, d, food, km, dm, eta, lk


def neighbours(cells, device):
    """[cells, 8] slot indices of the 8 torus neighbours of every slot."""
    W, H = GRID[cells]
    idx = torch.arange(cells, device=device)
    x, y = idx % W, idx // W
    out = [((x + dx) % W) + W * ((y + dy) % H) for dx in (-1, 0, 1) for dy in (-1, 0, 1) if dx or dy]
    return torch.stack(out, 1)


def h(q):
    return 2 * q / (1 + q)


@torch.no_grad()
def simulate(u, generator=None, level=None, return_stats=False, fidelity=1, knockout="module",
             trace_every=0, neutral_q=None):
    """One population per row of u in [0,1]^12. Returns bool [N] (detector level LEVEL, 1..4)."""
    ko = KNOCKOUT if knockout == "module" else knockout
    neutral_q = NEUTRAL_Q if neutral_q is None else neutral_q
    level = level or LEVEL
    cells, steps = FIDELITY[fidelity]
    g_min = G_MIN * steps / FIDELITY[1][1]
    dev = u.device
    B = u.shape[0]
    S, sigma, mu, K, d, food, km, dm, eta, lk = decode(u)
    nbr = neighbours(cells, dev)                                     # [C, 8]
    active = (torch.arange(K_MAX, device=dev).view(1, K_MAX) < K.view(B, 1)).view(B, 1, K_MAX)
    Kf = K.float().view(B, 1)
    col = lambda v: v.view(B, 1)
    col3 = lambda v: v.view(B, 1, 1)
    rate = col(food * DT)
    rnd = lambda *shape: torch.rand(*shape, device=dev, generator=generator)
    nrm = lambda *shape: torch.randn(*shape, device=dev, generator=generator)
    pois = lambda lam: torch.poisson(lam, generator=generator)
    binom = lambda n, p: torch.binomial(n, p.expand_as(n).contiguous(), generator=generator)

    X = torch.zeros(B, cells, K_MAX + 1, device=dev)
    per_gene = torch.clamp(torch.floor(S / (2 * K.float())), min=1.0).view(B, 1, 1)
    X[..., :K_MAX] = torch.where(active, per_gene, torch.zeros_like(per_gene)).expand(B, cells, K_MAX)
    M = (S / 2).view(B, 1).expand(B, cells).contiguous()
    q = torch.ones(B, cells, device=dev)
    qn = torch.ones(B, cells, M_NEUTRAL, device=dev)
    ndiv = torch.zeros(B, device=dev)
    div2 = torch.zeros(B, device=dev)
    lys2 = torch.zeros(B, device=dev)
    flat = (torch.arange(B, device=dev) * cells).view(B, 1)
    trace = []

    for step in range(steps):
        second_half = step >= steps // 2
        w, has_all = _growth(X, active, Kf)
        hq = torch.ones_like(q) if neutral_q else h(q)
        tot = X.sum(-1)
        crowd = (1.0 - tot / (BURST * M.clamp(min=1.0))).clamp(min=0.0)
        g = (w * crowd * hq * rate).unsqueeze(-1)
        new_genes = pois(X[..., :K_MAX] * g)
        errors = binom(new_genes, col3(mu))
        new_par = pois(X[..., K_MAX:] * g * col3(sigma))
        X[..., :K_MAX] += new_genes - errors
        X[..., K_MAX:] += new_par + errors.sum(-1, keepdim=True)
        X -= binom(X, col3(d))
        if bool((lk > 0).any()):
            out_ = binom(X, col3(lk))
            dest = torch.gather(nbr.unsqueeze(0).expand(B, cells, 8), 2,
                                (rnd(B, cells, 1) * 8).long().clamp(max=7)).squeeze(2)    # [B, C]
            X = X - out_
            X = X.scatter_add(1, dest.unsqueeze(-1).expand_as(out_), out_)
        M = M + pois(col(km) * X[..., 0] * rate * hq * has_all.float())
        M = M - binom(M, col(dm))
        tot = X.sum(-1)
        lysed = (tot > BURST * M.clamp(min=1.0)) & (tot > 0)
        if bool(lysed.any()):
            X = torch.where(lysed.unsqueeze(-1), torch.zeros_like(X), X)
            M = torch.where(lysed, torch.zeros_like(M), M)
            if second_half:
                lys2 += lysed.sum(1).float()
        if trace_every and step % trace_every == 0:
            _, ha = _growth(X, active, Kf)
            nfa = ha.float().sum(1).clamp(min=1)
            trace.append((step, ((q * ha).sum(1) / nfa).mean().item(),
                          ((qn.mean(-1) * ha).sum(1) / nfa).mean().item()))
        div = M >= col(S)
        if not bool(div.any()):
            continue
        if ko == "equal_split":
            half = torch.floor(X / 2)
            d1 = half + (((X - 2 * half) > 0) & (rnd(*X.shape) < 0.5)).float()
        else:
            d1 = binom(X, torch.full((1, 1, 1), 0.5, device=dev))
        d2 = X - d1
        Mh = M / 2
        # inheritance: both daughters mutate q and q' independently
        e = col(eta)
        q1, q2 = (q + e * nrm(B, cells)).clamp(0, Q_MAX), (q + e * nrm(B, cells)).clamp(0, Q_MAX)
        e3 = e.unsqueeze(-1)
        n1 = (qn + e3 * nrm(B, cells, M_NEUTRAL)).clamp(0, Q_MAX)
        n2 = (qn + e3 * nrm(B, cells, M_NEUTRAL)).clamp(0, Q_MAX)
        # local placement: pick a neighbour (dead/empty first unless no_selection; never a dividing one)
        nb = nbr.unsqueeze(0).expand(B, cells, 8)
        score = rnd(B, cells, 8)
        if ko != "no_selection":
            _, alive = _growth(X, active, Kf)
            score = score + 2.0 * (~torch.gather(alive, 1, nb.reshape(B, -1)).view(B, cells, 8)).float()
        score = score - 3.0 * torch.gather(div, 1, nb.reshape(B, -1)).view(B, cells, 8).float()
        target = torch.gather(nb, 2, score.argmax(2, keepdim=True)).squeeze(2)      # [B, C]
        # conflicts: one winner per target slot (random priority); losers' second daughter is lost
        prio = torch.where(div, rnd(B, cells), torch.full((B, cells), -1.0, device=dev))
        best = torch.full((B, cells), -1.0, device=dev).scatter_reduce(1, target, prio, reduce="amax")
        win = div & (torch.gather(best, 1, target) == prio) & ~torch.gather(div, 1, target)
        X = torch.where(div.unsqueeze(-1), d1, X)
        M = torch.where(div, Mh, M)
        q = torch.where(div, q1, q)
        qn = torch.where(div.unsqueeze(-1), n1, qn)
        src = (flat + torch.arange(cells, device=dev).view(1, cells))[win]
        dst = (flat + target)[win]
        X.view(B * cells, -1)[dst] = d2.view(B * cells, -1)[src]
        M.view(-1)[dst] = Mh.reshape(-1)[src]
        q.view(-1)[dst] = q2.reshape(-1)[src]
        qn.view(B * cells, -1)[dst] = n2.reshape(B * cells, -1)[src]
        n = div.sum(1).float()
        ndiv += n
        if second_half:
            div2 += n

    w, has_all = _growth(X, active, Kf)
    functional = has_all & (w >= W_MIN)
    nf = functional.float().sum(1)
    frac_alive = nf / cells
    gens = ndiv / cells
    ratio = torch.where(functional, X.sum(-1) / M.clamp(min=1.0), torch.zeros_like(M)).sum(1) / nf.clamp(min=1)
    lysis_frac = lys2 / div2.clamp(min=1)
    mq = torch.where(functional, q, torch.zeros_like(q)).sum(1) / nf.clamp(min=1)
    mqm = torch.where(functional.unsqueeze(-1), qn, torch.zeros_like(qn)).sum(1) / nf.clamp(min=1).unsqueeze(-1)
    mqn = mqm.mean(-1)
    z = (mq - mqn) / mqm.std(-1).clamp(min=1e-6)
    l1 = functional.any(1)
    l2 = (frac_alive >= 0.5) & (gens >= g_min)
    l3 = l2 & (ratio >= HOMEO_LO) & (ratio <= HOMEO_HI) & (lysis_frac < LYSIS_MAX)
    l4 = l3 & ((mq - mqn) > L4_DELTA) & (z > L4_Z)
    out = {1: l1, 2: l2, 3: l3}.get(level, l4)
    if return_stats:
        st = {"frac_alive": frac_alive, "gens": gens, "ratio": ratio, "lysis_frac": lysis_frac,
              "mean_q": mq, "mean_q_neutral": mqn, "dq": mq - mqn, "z": z, "l1": l1, "l2": l2, "l3": l3, "l4": l4}
        if trace_every:
            st["trace"] = trace
        return out, st
    return out
