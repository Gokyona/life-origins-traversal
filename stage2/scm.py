"""Stage-2a test problem: the stochastic corrector model (after Szathmary & Demeter 1987).

A population of C protocells (fixed number of slots). Each cell holds counts of
K essential gene types (mutually dependent templates) and one parasite type.
Per tau-leap step:
  * cell growth factor w = K * geometric mean of the K gene fractions (fractions of
    all molecules in the cell, parasites included); w = 0 if any gene type is
    missing (the cell is dead), w = 1 for a perfectly balanced, parasite-free cell
  * gene copies ~ Poisson(k * x_i * w * food * dt); each new copy is a parasite
    with probability mu (replication error destroys the gene's function)
  * parasite copies ~ Poisson(sigma * k * x_p * w * food * dt): parasites use the
    cell's machinery and add nothing to it
  * every molecule is degraded with probability d per step
  * a cell whose total reaches the division size S splits, each type
    Binomial(1/2); one daughter keeps the slot, the other replaces a dead or
    empty cell if there is one, otherwise a random cell (never a dividing one)

Selection acts between cells (balanced, parasite-poor cells grow and divide
faster), drift acts at division (small cells lose gene types by chance), and
within-cell competition favours parasites (sigma > 1). This is the corrector
mechanism of Szathmary & Demeter (1987).

Detector levels (module constant LEVEL). A cell is *functional* when it carries
all K gene types and its growth factor w >= W_MIN (not crippled by parasites).
  L1: at the end some cell is functional
  L2: at the end >= 50% of cells are functional, and the population has gone
      through >= G_MIN divisions per cell slot on average (heredity kept across
      many generations)

Deviations from Szathmary & Demeter (1987), kept deliberately simple:
  * fixed number of cell slots with replacement (Moran-like) instead of a
    growing / culled population
  * discrete tau-leap steps instead of an exact stochastic simulation
  * the growth factor is a geometric mean of gene fractions (their metabolic
    coupling is model-specific); parasites arise only by mutation of gene copies
  * K and the parasite advantage are free conditions (the original fixes them)

Conditions, each u in [0,1]:
  0 S      division size (molecules), log              8 .. 512
  1 sigma  parasite replication advantage              1 .. 3     (favourable -)
  2 mu     mutation-to-parasite probability, log       1e-3 .. 0.3 (favourable -)
  3 K      number of essential gene types, integer     1 .. 3
  4 d      degradation probability per step            0 .. 0.02  (favourable -)
  5 food   growth rate (resource supply), log          0.25 .. 1  (favourable +)
  6 z1     distractor
  7 z2     distractor
"""
import math
import torch

DIM = 8
NAMES = ["S", "sigma", "mu", "K", "d", "food", "z1", "z2"]
MONOTONE = [0, -1, -1, 0, -1, 1, 0, 0]
LEVEL = 2
CELLS = 128      # cell slots per population
STEPS = 1000     # tau-leap steps per simulation
DT = 0.05        # max gene growth per step = k * food * DT (k = 1)
G_MIN = 15       # divisions per cell slot required for L2
W_MIN = 0.3      # growth factor below which a cell counts as crippled
K_MAX = 3


def decode(u):
    S = 8.0 * 64.0 ** u[:, 0]
    sigma = 1.0 + 2.0 * u[:, 1]
    mu = 10 ** (-3.0 + math.log10(300.0) * u[:, 2])
    K = (1 + torch.floor(u[:, 3] * K_MAX)).clamp(max=K_MAX).long()
    d = 0.02 * u[:, 4]
    food = 0.25 * 4.0 ** u[:, 5]
    return S, sigma, mu, K, d, food


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
    """One population per row of u in [0,1]^8. Returns bool [N] (detector level LEVEL)."""
    level = level or LEVEL
    dev = u.device
    B = u.shape[0]
    S, sigma, mu, K, d, food = decode(u)
    active = (torch.arange(K_MAX, device=dev).view(1, K_MAX) < K.view(B, 1)).view(B, 1, K_MAX)
    Kf = K.float().view(B, 1)
    col = lambda v: v.view(B, 1)
    rate = col(food * DT)
    rnd = lambda *shape: torch.rand(*shape, device=dev, generator=generator)
    pois = lambda lam: torch.poisson(lam, generator=generator)
    binom = lambda n, p: torch.binomial(n, p.expand_as(n).contiguous(), generator=generator)

    # start: every cell at half the division size, genes balanced, no parasites
    X = torch.zeros(B, CELLS, K_MAX + 1, device=dev)
    per_gene = torch.clamp(torch.floor(S / (2 * K.float())), min=1.0).view(B, 1, 1)
    X[..., :K_MAX] = torch.where(active, per_gene, torch.zeros_like(per_gene)).expand(B, CELLS, K_MAX)
    ndiv = torch.zeros(B, device=dev)
    flat = (torch.arange(B, device=dev) * CELLS).view(B, 1)

    for _ in range(STEPS):
        w, has_all = _growth(X, active, Kf)
        g = (w * rate).unsqueeze(-1)                                   # [B,C,1]
        new_genes = pois(X[..., :K_MAX] * g)
        errors = binom(new_genes, col(mu).view(B, 1, 1))
        new_par = pois(X[..., K_MAX:] * g * col(sigma).view(B, 1, 1))
        X[..., :K_MAX] += new_genes - errors
        X[..., K_MAX:] += new_par + errors.sum(-1, keepdim=True)
        X -= binom(X, col(d).view(B, 1, 1))

        tot = X.sum(-1)
        div = tot >= col(S)
        if not bool(div.any()):
            continue
        d1 = binom(X, torch.full((1, 1, 1), 0.5, device=dev))
        d2 = X - d1
        _, alive = _growth(X, active, Kf)
        score = rnd(B, CELLS) + 2.0 * (~alive).float()     # dead / empty slots first
        score = torch.where(div, torch.full_like(score, -1.0), score)
        order = torch.argsort(score, dim=1, descending=True)
        rank = (torch.cumsum(div.long(), 1) - 1).clamp(min=0)
        target = torch.gather(order, 1, rank)
        X = torch.where(div.unsqueeze(-1), d1, X)
        src = (flat + torch.arange(CELLS, device=dev).view(1, CELLS))[div]
        dst = (flat + target)[div]
        Xf = X.view(B * CELLS, -1)
        Xf[dst] = d2.view(B * CELLS, -1)[src]
        ndiv += div.sum(1).float()

    w, has_all = _growth(X, active, Kf)
    functional = has_all & (w >= W_MIN)
    frac_alive = functional.float().mean(1)
    gens = ndiv / CELLS
    l1 = functional.any(1)
    l2 = (frac_alive >= 0.5) & (gens >= G_MIN)
    out = l1 if level == 1 else l2
    if return_stats:
        return out, {"frac_alive": frac_alive, "gens": gens, "l1": l1, "l2": l2}
    return out
