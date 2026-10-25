"""Stage-3a: mechanism knockouts of the stage-2b protocell model.

Same interface and conditions as protocell.py (DIM, NAMES, MONOTONE, LEVEL,
simulate(u, generator, level, return_stats, fidelity)), plus a module constant
KNOCKOUT (or the knockout= argument) selecting one mechanism to remove:

  None              the stage-2b model unchanged
  "equal_split"     division halves every molecule type exactly (odd counts are
                    rounded up or down at random): no assortment noise, so the
                    stochastic corrector has no variation between daughters
                    to select on
  "no_compartments" one well-mixed soup per population instead of CELLS cells:
                    same total starting material (CELLS * S/2 molecules),
                    crowding against a fixed capacity CELLS * S as in a cell,
                    and a chemostat outflow that thins the soup back to 1.5x
                    its starting size (the mean fill of a cell population), so
                    molecules keep turning over; no lysis, membrane or division. Detectors: the soup
                    is "functional" like a cell; L2 needs it functional and
                    G_MIN generations, where a generation is a doubling of the
                    starting material (new copies / starting size); L3 is
                    impossible (no protocells) and reported as False
  "no_selection"    the second daughter replaces a uniformly random non-dividing
                    slot, ignoring whether it holds a dead/empty cell (removes
                    the bias that lets working cells take over dead slots)
  "no_coupling"     division at content size S as in scm.py; membrane ignored
                    (no capacity, no crowding, no lysis). L3 is not defined
                    and reported as False
  "no_decay"        d = dm = 0 regardless of the conditions
  "K1"              K = 1 regardless of the conditions (one gene does all)
"""
import os
import sys

import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stage2"))
import protocell as base  # noqa: E402
from protocell import (DIM, NAMES, MONOTONE, CELLS, STEPS, DT, G_MIN, W_MIN, BURST, HOMEO_LO, HOMEO_HI,  # noqa: F401
                       LYSIS_MAX, K_MAX, FIDELITY, decode, _growth)

LEVEL = 3
KNOCKOUT = None
KNOCKOUTS = [None, "equal_split", "no_compartments", "no_selection", "no_coupling", "no_decay", "K1"]


@torch.no_grad()
def simulate(u, generator=None, level=None, return_stats=False, fidelity=1, knockout="module"):
    """One population per row of u in [0,1]^10. Returns bool [N] (detector level LEVEL)."""
    ko = KNOCKOUT if knockout == "module" else knockout
    level = level or LEVEL
    cells, steps = FIDELITY[fidelity]
    g_min = G_MIN * steps / FIDELITY[1][1]
    dev = u.device
    B = u.shape[0]
    S, sigma, mu, K, d, food, km, dm = decode(u)
    if ko == "no_decay":
        d, dm = torch.zeros_like(d), torch.zeros_like(dm)
    if ko == "K1":
        K = torch.ones_like(K)
    soup = ko == "no_compartments"
    n_slots = 1 if soup else cells
    active = (torch.arange(K_MAX, device=dev).view(1, K_MAX) < K.view(B, 1)).view(B, 1, K_MAX)
    Kf = K.float().view(B, 1)
    col = lambda v: v.view(B, 1)
    col3 = lambda v: v.view(B, 1, 1)
    rate = col(food * DT)
    rnd = lambda *shape: torch.rand(*shape, device=dev, generator=generator)
    pois = lambda lam: torch.poisson(lam, generator=generator)
    binom = lambda n, p: torch.binomial(n, p.expand_as(n).contiguous(), generator=generator)

    scale = cells if soup else 1          # the soup holds all cells' starting material
    X = torch.zeros(B, n_slots, K_MAX + 1, device=dev)
    per_gene = torch.clamp(torch.floor(S / (2 * K.float())), min=1.0).view(B, 1, 1) * scale
    X[..., :K_MAX] = torch.where(active, per_gene, torch.zeros_like(per_gene)).expand(B, n_slots, K_MAX)
    start_size = X.sum(-1).mean(1)        # for soup generations
    M = (S / 2 * scale).view(B, 1).expand(B, n_slots).contiguous()
    ndiv = torch.zeros(B, device=dev)
    div2 = torch.zeros(B, device=dev)
    lys2 = torch.zeros(B, device=dev)
    new_total = torch.zeros(B, device=dev)
    flat = (torch.arange(B, device=dev) * n_slots).view(B, 1)

    for step in range(steps):
        second_half = step >= steps // 2
        w, has_all = _growth(X, active, Kf)
        tot = X.sum(-1)
        if ko == "no_coupling":
            crowd = torch.ones_like(tot)
        else:
            crowd = (1.0 - tot / (BURST * M.clamp(min=1.0))).clamp(min=0.0)
        g = (w * crowd * rate).unsqueeze(-1)
        new_genes = pois(X[..., :K_MAX] * g)
        errors = binom(new_genes, col3(mu))
        new_par = pois(X[..., K_MAX:] * g * col3(sigma))
        X[..., :K_MAX] += new_genes - errors
        X[..., K_MAX:] += new_par + errors.sum(-1, keepdim=True)
        X -= binom(X, col3(d))
        if soup:
            new_total += (new_genes.sum(-1) + new_par.sum(-1)).sum(1)
            # chemostat outflow: thin the soup back to 1.5x its starting size (the mean fill of
            # a cell population between S/2 and S), so molecules keep turning over
            keep = (1.5 * start_size / X.sum(-1).sum(1).clamp(min=1)).clamp(max=1.0)
            X = binom(X, keep.view(B, 1, 1))
            continue                      # no membrane, no lysis, no division

        if ko != "no_coupling":
            M = M + pois(col(km) * X[..., 0] * rate * has_all.float())
            M = M - binom(M, col(dm))
            tot = X.sum(-1)
            lysed = (tot > BURST * M.clamp(min=1.0)) & (tot > 0)
            if bool(lysed.any()):
                X = torch.where(lysed.unsqueeze(-1), torch.zeros_like(X), X)
                M = torch.where(lysed, torch.zeros_like(M), M)
                if second_half:
                    lys2 += lysed.sum(1).float()
            div = M >= col(S)
        else:
            div = X.sum(-1) >= col(S)
        if not bool(div.any()):
            continue
        if ko == "equal_split":
            half = torch.floor(X / 2)
            odd = (X - 2 * half) > 0
            d1 = half + (odd & (rnd(*X.shape) < 0.5)).float()
        else:
            d1 = binom(X, torch.full((1, 1, 1), 0.5, device=dev))
        d2 = X - d1
        Mh = M / 2
        if ko == "no_selection":
            score = rnd(B, n_slots)
        else:
            _, alive = _growth(X, active, Kf)
            score = rnd(B, n_slots) + 2.0 * (~alive).float()     # dead / empty slots first
        score = torch.where(div, torch.full_like(score, -1.0), score)
        order = torch.argsort(score, dim=1, descending=True)
        rank = (torch.cumsum(div.long(), 1) - 1).clamp(min=0)
        target = torch.gather(order, 1, rank)
        X = torch.where(div.unsqueeze(-1), d1, X)
        M = torch.where(div, Mh, M)
        src = (flat + torch.arange(n_slots, device=dev).view(1, n_slots))[div]
        dst = (flat + target)[div]
        X.view(B * n_slots, -1)[dst] = d2.view(B * n_slots, -1)[src]
        M.view(-1)[dst] = Mh.reshape(-1)[src]
        n = div.sum(1).float()
        ndiv += n
        if second_half:
            div2 += n

    w, has_all = _growth(X, active, Kf)
    functional = has_all & (w >= W_MIN)
    nf = functional.float().sum(1)
    frac_alive = nf / n_slots
    gens = new_total / start_size.clamp(min=1) if soup else ndiv / n_slots
    ratio = torch.where(functional, X.sum(-1) / M.clamp(min=1.0), torch.zeros_like(M)).sum(1) / nf.clamp(min=1)
    lysis_frac = lys2 / div2.clamp(min=1)
    l1 = functional.any(1)
    l2 = (frac_alive >= 0.5) & (gens >= g_min)
    if ko in ("no_compartments", "no_coupling"):
        l3 = torch.zeros_like(l2)
    else:
        l3 = l2 & (ratio >= HOMEO_LO) & (ratio <= HOMEO_HI) & (lysis_frac < LYSIS_MAX)
    out = {1: l1, 2: l2}.get(level, l3)
    if return_stats:
        # between-cell variation: variance of the gene-0 fraction (of all molecules) across functional cells
        f0 = X[..., 0] / X.sum(-1).clamp(min=1.0)
        fm = torch.where(functional, f0, torch.zeros_like(f0)).sum(1) / nf.clamp(min=1)
        var_g0 = torch.where(functional, (f0 - fm.unsqueeze(1)) ** 2, torch.zeros_like(f0)).sum(1) / (nf - 1).clamp(min=1)
        return out, {"frac_alive": frac_alive, "gens": gens, "ratio": ratio, "lysis_frac": lysis_frac,
                     "var_g0": var_g0, "l1": l1, "l2": l2, "l3": l3}
    return out
