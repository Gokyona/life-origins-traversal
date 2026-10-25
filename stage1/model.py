"""Minimal stochastic "protocell chemostat" model (stage-0 test problem).

A well-mixed flow reactor holds food F, a master replicator R and junk
mutants M.  R can appear spontaneously from food (the rare "origin" event),
copies itself using food with fidelity Q = exp(-muL), and loses copies to
dilution and temperature-dependent hydrolysis.  Mutants copy at half speed.
A run counts as "life emerged" if, at the end, R holds at least R_SUCCESS
copies plus five times the abiotic background (what spontaneous formation
alone would sustain), so that only self-replication can pass.

The six conditions live in the unit cube u in [0,1]^6 and map to:
  0 Da    log10 replication/dilution ratio      [-0.5, 2.0]   (favourable +)
  1 muL   error number mu*L                       [ 0.0, 1.5]   (favourable -)
  2 phi   log10 food inflow per step              [ 0.5, 2.5]   (favourable +)
  3 s     log10 spontaneous formation rate        [-8.5,-4.5]   (favourable +)
  4 T     temperature (normalised)                [ 0.0, 1.0]   (non-monotone)
  5 z     distractor: has no effect on the model  [ 0.0, 1.0]
"""
import math
import torch

DIM = 6
NAMES = ["Da", "muL", "phi", "s", "T", "z"]
# +1: larger is better, -1: smaller is better, 0: unknown / not declared
MONOTONE = [1, -1, 1, 1, 0, 0]

STEPS = 800
DILUTION = 0.01
F_SCALE = 1000.0
R_SUCCESS = 30


def decode(u):
    """Unit-cube conditions -> physical parameters (all tensors of shape [N])."""
    da = 10 ** (-0.5 + 2.5 * u[:, 0])
    mul = 1.5 * u[:, 1]
    phi = 10 ** (0.5 + 2.0 * u[:, 2])
    s = 10 ** (-8.5 + 4.0 * u[:, 3])
    t = u[:, 4]
    return da, mul, phi, s, t


def _pois(rate):
    return torch.poisson(rate.clamp(min=0.0))


@torch.no_grad()
def simulate(u, generator=None, steps=STEPS):
    """Run one stochastic trajectory per row of u. Returns bool tensor [N]."""
    da, mul, phi, s, t = decode(u)
    k = da * DILUTION * torch.exp(-((t - 0.5) ** 2) / (2 * 0.25 ** 2))
    q = torch.exp(-mul)
    hyd = 0.001 * torch.exp(3.0 * t)
    loss = (DILUTION + hyd).clamp(max=1.0)

    n = u.shape[0]
    F = phi / DILUTION  # start at food steady state, no replicators
    R = torch.zeros(n, device=u.device)
    M = torch.zeros(n, device=u.device)
    pois = (lambda r: torch.poisson(r.clamp(min=0.0), generator=generator))
    for _ in range(steps):
        rep_r = torch.minimum(pois(k * R * F / F_SCALE), F)
        faithful = torch.minimum(pois(rep_r * q), rep_r)
        rep_m = torch.minimum(pois(0.5 * k * M * F / F_SCALE), F - rep_r)
        spont = pois(s * F)
        F = (F - rep_r - rep_m - spont).clamp(min=0.0)
        R = R + faithful + spont
        M = M + (rep_r - faithful) + rep_m
        R = R - torch.minimum(pois(R * loss), R)
        M = M - torch.minimum(pois(M * loss), M)
        F = F - torch.minimum(pois(F * DILUTION), F) + pois(phi)
    # life = replication amplifies R far above the abiotic background s*F/loss
    background = s * (phi / DILUTION) / loss
    return R >= R_SUCCESS + 5.0 * background
