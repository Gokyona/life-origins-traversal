"""Binary polymer model with random catalysis (Kauffman 1986; Hordijk & Steel 2004).

Molecules are all binary strings of length 1..n. Every split of a string ab into
a + b is one reversible reaction (ligation a + b -> ab / cleavage ab -> a + b).
The food set is every string of length <= t. Each molecule catalyses each
reaction independently with probability p; the catalysis level
f = p * (number of reactions) is the mean number of reactions a molecule catalyses.
A fraction `dens` of the possible reactions exists at all.

A run "produces life" (level L1: a self-sustaining autocatalytic network) when the
instance contains a RAF set: reflexively autocatalytic (every reaction catalysed
by a molecule the network itself makes) and food-generated. Detected exactly with
the Hordijk-Steel maxRAF reduction, batched on the GPU.

Theory: for n up to ~20 a RAF appears with high probability once f is about 1-2
(Hordijk & Steel 2004; Mossel & Steel 2005), so the life boundary is known in
advance and the engine's map can be checked against it.

Conditions, each u in [0,1]:
  0 n     maximum polymer length, integer             5 .. 10
  1 f     log10 catalysis level                      -1 .. 1    (favourable +)
  2 t     food set: strings up to length t, integer   1 .. 3     (favourable +)
  3 dens  fraction of possible reactions present     0.3 .. 1   (favourable +)
  4 c     fraction of molecules able to catalyse      0.1 .. 1   (f is kept fixed, so
          catalysis concentrates on fewer molecules)
  5 z1    distractor
  6 z2    distractor
"""
import functools
import math
import torch

DIM = 7
NAMES = ["n", "f", "t", "dens", "c", "z1", "z2"]
MONOTONE = [0, 1, 1, 1, 0, 0, 0]
N_MIN, N_MAX = 5, 10
T_MIN, T_MAX = 1, 3
K_SLOTS = 10  # max catalysts stored per reaction; Poisson mean is <= ~3.2 here


def decode(u):
    n = (N_MIN + torch.floor(u[:, 0] * (N_MAX - N_MIN + 1))).clamp(max=N_MAX).long()
    f = 10 ** (-1.0 + 2.0 * u[:, 1])
    t = (T_MIN + torch.floor(u[:, 2] * (T_MAX - T_MIN + 1))).clamp(max=T_MAX).long()
    dens = 0.3 + 0.7 * u[:, 3]
    c = 0.1 + 0.9 * u[:, 4]
    return n, f, t, dens, c


def _index(s):
    """Molecule id of a binary string: strings ordered by length, then value."""
    return (1 << len(s)) - 2 + int(s, 2)


@functools.lru_cache(maxsize=None)
def network(n, device):
    """Reaction list for max length n: tensors a, b, ab of molecule ids, and M."""
    a, b, ab = [], [], []
    for L in range(2, n + 1):
        for v in range(1 << L):
            s = format(v, f"0{L}b")
            for cut in range(1, L):
                a.append(_index(s[:cut]))
                b.append(_index(s[cut:]))
                ab.append(_index(s))
    M = (1 << (n + 1)) - 2
    t = lambda z: torch.tensor(z, dtype=torch.long, device=device)
    return t(a), t(b), t(ab), M


def _max_raf(food, present, cat, a, b, ab, M):
    """Hordijk-Steel reduction for a batch. Returns bool [B]: maxRAF non-empty.

    food [B,M] bool, present [B,R] bool (reaction exists and has >= 1 catalyst),
    cat [B,R,K] long catalyst ids (-1 = empty slot).
    """
    B, R = present.shape
    active = present.clone()
    cat_valid = cat >= 0
    cat_idx = cat.clamp(min=0)
    while True:
        # closure of the food set under the active (reversible) reactions
        avail = food.clone()
        while True:
            lig = active & avail[:, a] & avail[:, b]
            cle = active & avail[:, ab]
            new = avail.to(torch.uint8)  # scatter kernels do not take bool
            new.scatter_reduce_(1, ab.expand(B, R), lig.to(torch.uint8), reduce="amax")
            cle = cle.to(torch.uint8)
            new.scatter_reduce_(1, a.expand(B, R), cle, reduce="amax")
            new.scatter_reduce_(1, b.expand(B, R), cle, reduce="amax")
            new = new.bool()
            if torch.equal(new, avail):
                break
            avail = new
        # keep reactions that are catalysed by, and have reactants in, the closure
        catalysed = (torch.gather(avail, 1, cat_idx.view(B, -1)).view(B, R, -1) & cat_valid).any(2)
        reach = (avail[:, a] & avail[:, b]) | avail[:, ab]
        new_active = active & catalysed & reach
        if torch.equal(new_active, active):
            return active.any(1)
        active = new_active


@torch.no_grad()
def simulate(u, generator=None, max_batch=None):
    """One random catalysis assignment per row of u in [0,1]^7. Returns bool [N]."""
    dev = u.device
    n, f, t, dens, c = decode(u)
    out = torch.zeros(u.shape[0], dtype=torch.bool, device=dev)
    for nv in torch.unique(n).tolist():
        a, b, ab, M = network(nv, dev)
        R = a.numel()
        # keep the per-batch catalyst table around ~50M entries
        cap = max_batch or max(1, int(5e7 // (R * K_SLOTS)))
        for tv in torch.unique(t[n == nv]).tolist():
            rows = torch.nonzero((n == nv) & (t == tv)).squeeze(1)
            food_n = (1 << (tv + 1)) - 2
            for i in range(0, rows.numel(), cap):
                r = rows[i:i + cap]
                Bn = r.numel()
                lam = (f[r] / R * M).unsqueeze(1).expand(Bn, R)  # mean catalysts per reaction
                k = torch.poisson(lam.contiguous(), generator=generator).clamp(max=K_SLOTS).long()
                exists = torch.rand(Bn, R, device=dev, generator=generator) < dens[r].unsqueeze(1)
                # catalysts come from a random subset of round(c*M) molecules per instance
                perm = torch.argsort(torch.rand(Bn, M, device=dev, generator=generator), 1)
                m_c = (c[r] * M).round().clamp(min=1).view(Bn, 1, 1)
                pick = (torch.rand(Bn, R, K_SLOTS, device=dev, generator=generator) * m_c).long()
                ids = torch.gather(perm, 1, pick.view(Bn, -1)).view(Bn, R, K_SLOTS)
                slot = torch.arange(K_SLOTS, device=dev)
                cat = torch.where(slot < k.unsqueeze(2), ids, torch.full_like(ids, -1))
                present = exists & (k > 0)
                food = torch.zeros(Bn, M, dtype=torch.bool, device=dev)
                food[:, :food_n] = True
                out[r] = _max_raf(food, present, cat, a, b, ab, M)
    return out


if __name__ == "__main__":
    # sanity check against theory: P(RAF) vs catalysis level f, food t=2, all reactions
    # present, every molecule a potential catalyst
    import argparse, time
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--seeds", type=int, default=200)
    a_ = ap.parse_args()
    dev = torch.device(a_.device)
    g = torch.Generator(device=dev).manual_seed(0)
    fs = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
    print("P(RAF), food t=2, dens=1;  rows n, columns f =", fs)
    for nv in range(N_MIN, N_MAX + 1):
        t0 = time.time()
        row = []
        for fv in fs:
            u = torch.zeros(a_.seeds, DIM, device=dev)
            u[:, 0] = (nv - N_MIN + 0.5) / (N_MAX - N_MIN + 1)
            u[:, 1] = (math.log10(fv) + 1) / 2
            u[:, 2] = (2 - T_MIN + 0.5) / (T_MAX - T_MIN + 1)
            u[:, 3] = 1.0
            u[:, 4] = 1.0
            row.append(simulate(u, generator=g).float().mean().item())
        print(f"n={nv:2d} " + " ".join(f"{p:5.2f}" for p in row) + f"   ({time.time()-t0:.1f}s)")
