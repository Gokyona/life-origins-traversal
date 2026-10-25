"""Stage-2b test problem: stochastic-corrector protocells whose division follows from
their own chemistry (after Szathmary & Demeter 1987, with Ganti's chemoton idea that
template replication, metabolism and membrane growth are coupled).

Extends scm.py. Each cell now also carries a membrane amount M (lipid units):
  * gene type 0 is the membrane enzyme (one of the K essential genes; with K = 1
    the single gene does both jobs). Membrane synthesis per step
    ~ Poisson(km * x_0 * food * dt) in a working cell (all K genes present);
    every lipid unit turns over with probability dm per step
  * capacity: a cell holds about M molecules comfortably (a newborn cell has
    M0 = S/2 lipid units and S/2 molecules). Replication slows with crowding,
    w_eff = w * max(0, 1 - tot / (BURST * M)); a cell with tot > BURST * M lyses
    (dies, contents and membrane lost)
  * division is triggered ONLY by the membrane: when M >= S (= 2 * M0) the cell
    splits, each daughter gets M/2 and Binomial(1/2) of every molecule type;
    one daughter keeps the slot, the other replaces a dead/empty cell first,
    otherwise a random non-dividing cell (as in scm.py)
  * everything else (gene/parasite copying, mutation to parasite, degradation d)
    is as in scm.py, with w replaced by w_eff
Too slow a membrane -> crowded cells, slow division and bursting; too fast ->
vesicles split before their contents double and become under-filled (content /
capacity falls below HOMEO_LO; at S >= 16 genes are rarely lost outright). So
there is a window in km for L3.

Detector levels (module constant LEVEL, default 3). A cell is *functional* when
it carries all K gene types and w >= W_MIN (lysed cells are empty, so excluded).
  L1: at the end some cell is functional
  L2: at the end >= 50% of cells are functional, and the population has gone
      through >= G_MIN divisions per cell slot on average
  L3: L2, plus homeostasis: the mean content/capacity ratio (tot / M) of
      functional cells at the end lies in [HOMEO_LO, HOMEO_HI], and lysis events
      in the second half of the run are < LYSIS_MAX of division events there.
      (Lysis is counted per division rather than per cell death: with fixed
      slots, every division also removes a cell, so divisions are the natural
      count of deaths by replacement. The mean replaces the median because
      nanmedian is not reliable on XPU.)

Choices made during validation: DT = 0.1 (with 0.05, low-food populations could
not reach G_MIN in STEPS steps and P(L3) over the space was only 5%); km spans
0.1 .. 30 so that both edges of the L3 window fall inside the space.

Conditions, each u in [0,1]:
  0 S      membrane division size (lipid units), log   8 .. 512
  1 sigma  parasite replication advantage              1 .. 3      (favourable -)
  2 mu     mutation-to-parasite probability, log       1e-3 .. 0.3 (favourable -)
  3 K      number of essential gene types, integer     1 .. 3
  4 d      degradation probability per step            0 .. 0.02   (favourable -)
  5 food   growth rate (resource supply), log          0.25 .. 1   (favourable +)
  6 km     membrane synthesis rate per enzyme, log     0.1 .. 30   (window)
  7 dm     membrane turnover probability per step      0 .. 0.02   (favourable -)
  8 z1     distractor
  9 z2     distractor
"""
import math
import torch

DIM = 10
NAMES = ["S", "sigma", "mu", "K", "d", "food", "km", "dm", "z1", "z2"]
MONOTONE = [0, -1, -1, 0, -1, 1, 0, -1, 0, 0]
LEVEL = 3
CELLS = 128      # cell slots per population
STEPS = 1000     # tau-leap steps per simulation
DT = 0.1         # max gene growth per step = k * food * DT (k = 1)
G_MIN = 15       # divisions per cell slot required for L2
W_MIN = 0.3      # growth factor below which a cell counts as crippled
BURST = 2.0      # a cell lyses when its content exceeds BURST * capacity
HOMEO_LO, HOMEO_HI = 0.25, 1.5
LYSIS_MAX = 0.1
K_MAX = 3


def decode(u):
    S = 8.0 * 64.0 ** u[:, 0]
    sigma = 1.0 + 2.0 * u[:, 1]
    mu = 10 ** (-3.0 + math.log10(300.0) * u[:, 2])
    K = (1 + torch.floor(u[:, 3] * K_MAX)).clamp(max=K_MAX).long()
    d = 0.02 * u[:, 4]
    food = 0.25 * 4.0 ** u[:, 5]
    km = 10 ** (-1.0 + math.log10(300.0) * u[:, 6])
    dm = 0.02 * u[:, 7]
    return S, sigma, mu, K, d, food, km, dm


def _growth(X, active, Kf):
    """Cell growth factor w [B,C] and has-all-genes mask from counts X [B,C,K_MAX+1]."""
    tot = X.sum(-1)
    genes = X[..., :K_MAX]
    has_all = ((genes > 0) | ~active).all(-1) & (tot > 0)
    frac = genes / tot.clamp(min=1.0).unsqueeze(-1)
    logf = torch.where(active, torch.log(frac.clamp(min=1e-12)), torch.zeros_like(frac))
    w = Kf * torch.exp(logf.sum(-1) / Kf)
    return torch.where(has_all, w, torch.zeros_like(w)), has_all


@torch.no_grad()
def simulate(u, generator=None, level=None, return_stats=False):
    """One population per row of u in [0,1]^10. Returns bool [N] (detector level LEVEL)."""
    level = level or LEVEL
    dev = u.device
    B = u.shape[0]
    S, sigma, mu, K, d, food, km, dm = decode(u)
    active = (torch.arange(K_MAX, device=dev).view(1, K_MAX) < K.view(B, 1)).view(B, 1, K_MAX)
    Kf = K.float().view(B, 1)
    col = lambda v: v.view(B, 1)
    col3 = lambda v: v.view(B, 1, 1)
    rate = col(food * DT)
    rnd = lambda *shape: torch.rand(*shape, device=dev, generator=generator)
    pois = lambda lam: torch.poisson(lam, generator=generator)
    binom = lambda n, p: torch.binomial(n, p.expand_as(n).contiguous(), generator=generator)

    # start: every cell newborn (M0 = S/2 lipids, S/2 molecules), genes balanced
    X = torch.zeros(B, CELLS, K_MAX + 1, device=dev)
    per_gene = torch.clamp(torch.floor(S / (2 * K.float())), min=1.0).view(B, 1, 1)
    X[..., :K_MAX] = torch.where(active, per_gene, torch.zeros_like(per_gene)).expand(B, CELLS, K_MAX)
    M = (S / 2).view(B, 1).expand(B, CELLS).contiguous()
    ndiv = torch.zeros(B, device=dev)
    div2 = torch.zeros(B, device=dev)
    lys2 = torch.zeros(B, device=dev)
    flat = (torch.arange(B, device=dev) * CELLS).view(B, 1)

    for step in range(STEPS):
        second_half = step >= STEPS // 2
        w, has_all = _growth(X, active, Kf)
        tot = X.sum(-1)
        crowd = (1.0 - tot / (BURST * M.clamp(min=1.0))).clamp(min=0.0)
        g = (w * crowd * rate).unsqueeze(-1)                            # [B,C,1]
        new_genes = pois(X[..., :K_MAX] * g)
        errors = binom(new_genes, col3(mu))
        new_par = pois(X[..., K_MAX:] * g * col3(sigma))
        X[..., :K_MAX] += new_genes - errors
        X[..., K_MAX:] += new_par + errors.sum(-1, keepdim=True)
        X -= binom(X, col3(d))

        # membrane: synthesis by gene 0 in working cells, turnover
        M = M + pois(col(km) * X[..., 0] * rate * has_all.float())
        M = M - binom(M, col(dm))

        # lysis: content beyond what the membrane can hold
        tot = X.sum(-1)
        lysed = (tot > BURST * M.clamp(min=1.0)) & (tot > 0)
        if bool(lysed.any()):
            X = torch.where(lysed.unsqueeze(-1), torch.zeros_like(X), X)
            M = torch.where(lysed, torch.zeros_like(M), M)
            if second_half:
                lys2 += lysed.sum(1).float()

        # division: triggered by the membrane only
        div = M >= col(S)
        if not bool(div.any()):
            continue
        d1 = binom(X, torch.full((1, 1, 1), 0.5, device=dev))
        d2 = X - d1
        Mh = M / 2
        _, alive = _growth(X, active, Kf)
        score = rnd(B, CELLS) + 2.0 * (~alive).float()     # dead / empty slots first
        score = torch.where(div, torch.full_like(score, -1.0), score)
        order = torch.argsort(score, dim=1, descending=True)
        rank = (torch.cumsum(div.long(), 1) - 1).clamp(min=0)
        target = torch.gather(order, 1, rank)
        X = torch.where(div.unsqueeze(-1), d1, X)
        M = torch.where(div, Mh, M)
        src = (flat + torch.arange(CELLS, device=dev).view(1, CELLS))[div]
        dst = (flat + target)[div]
        X.view(B * CELLS, -1)[dst] = d2.view(B * CELLS, -1)[src]
        M.view(-1)[dst] = Mh.reshape(-1)[src]
        n = div.sum(1).float()
        ndiv += n
        if second_half:
            div2 += n

    w, has_all = _growth(X, active, Kf)
    functional = has_all & (w >= W_MIN)
    nf = functional.float().sum(1)
    frac_alive = nf / CELLS
    gens = ndiv / CELLS
    ratio = torch.where(functional, X.sum(-1) / M.clamp(min=1.0), torch.zeros_like(M)).sum(1) / nf.clamp(min=1)
    lysis_frac = lys2 / div2.clamp(min=1)
    l1 = functional.any(1)
    l2 = (frac_alive >= 0.5) & (gens >= G_MIN)
    l3 = l2 & (ratio >= HOMEO_LO) & (ratio <= HOMEO_HI) & (lysis_frac < LYSIS_MAX)
    out = {1: l1, 2: l2}.get(level, l3)
    if return_stats:
        return out, {"frac_alive": frac_alive, "gens": gens, "ratio": ratio, "lysis_frac": lysis_frac,
                     "l1": l1, "l2": l2, "l3": l3}
    return out
