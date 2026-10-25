"""Stage-1b test problem: a harder, 10-condition version of the protocell chemostat.

Same core as stage1/model.py (food F, master replicator R, mutants M), plus:
parasites whose fitness relative to the master is a condition, a tunable base
hydrolysis rate, wet/dry cycling, and two distractor dimensions. The life
region is thinner and the boundary more curved than in stage 1.

  0 Da     log10 replication/dilution ratio      [-0.5, 2.0]   (favourable +)
  1 muL    error number mu*L                       [ 0.0, 1.5]   (favourable -)
  2 phi    log10 food inflow per step              [ 0.5, 2.5]   (favourable +)
  3 s      log10 spontaneous formation rate        [-8.5,-4.5]   (favourable +)
  4 T      temperature (normalised)                [ 0.0, 1.0]   (non-monotone)
  5 h0     log10 base hydrolysis rate              [-3.5,-1.5]   (favourable -)
  6 sigma  parasite fitness / master fitness       [ 0.2, 1.2]   (favourable -)
  7 A      wet/dry cycle amplitude                 [ 0.0, 1.0]   (unknown)
  8 z1     distractor
  9 z2     distractor
"""
import math
import torch

DIM = 10
NAMES = ["Da", "muL", "phi", "s", "T", "h0", "sigma", "A", "z1", "z2"]
MONOTONE = [1, -1, 1, 1, 0, -1, -1, 0, 0, 0]

STEPS = 800
DILUTION = 0.01
F_SCALE = 1000.0
R_SUCCESS = 30
CYCLE = 100


@torch.no_grad()
def simulate(u, generator=None, steps=STEPS):
    """One stochastic trajectory per row of u in [0,1]^10. Returns bool [N]."""
    da = 10 ** (-0.5 + 2.5 * u[:, 0])
    mul = 1.5 * u[:, 1]
    phi = 10 ** (0.5 + 2.0 * u[:, 2])
    s = 10 ** (-8.5 + 4.0 * u[:, 3])
    t = u[:, 4]
    h0 = 10 ** (-3.5 + 2.0 * u[:, 5])
    sigma = 0.2 + 1.0 * u[:, 6]
    amp = u[:, 7]

    k = da * DILUTION * torch.exp(-((t - 0.5) ** 2) / (2 * 0.25 ** 2))
    q = torch.exp(-mul)
    hyd = h0 * torch.exp(3.0 * t)

    n = u.shape[0]
    F = phi / DILUTION
    R = torch.zeros(n, device=u.device)
    M = torch.zeros(n, device=u.device)
    pois = (lambda r: torch.poisson(r.clamp(min=0.0), generator=generator))
    for step in range(steps):
        # dry phase (c>0): copying speeds up, hydrolysis slows; wet phase the reverse
        c = math.sin(2 * math.pi * step / CYCLE)
        kk = k * (1 + 0.8 * amp * c)
        loss = (DILUTION + hyd * (1 - 0.8 * amp * c)).clamp(max=1.0)
        rep_r = torch.minimum(pois(kk * R * F / F_SCALE), F)
        faithful = torch.minimum(pois(rep_r * q), rep_r)
        rep_m = torch.minimum(pois(sigma * kk * M * F / F_SCALE), F - rep_r)
        spont = pois(s * F)
        F = (F - rep_r - rep_m - spont).clamp(min=0.0)
        R = R + faithful + spont
        M = M + (rep_r - faithful) + rep_m
        R = R - torch.minimum(pois(R * loss), R)
        M = M - torch.minimum(pois(M * loss), M)
        F = F - torch.minimum(pois(F * DILUTION), F) + pois(phi)
    background = s * (phi / DILUTION) / (DILUTION + hyd)
    return R >= R_SUCCESS + 5.0 * background
