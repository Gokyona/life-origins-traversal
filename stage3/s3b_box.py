"""Stage-3b: minimal sufficient condition boxes for L1, L2, L3 on the protocell model.

1. Map: COMPASS (engine4_soft, BUDGET runs, --fast) per level; final map-maker
   (fit_es_multi, hidden 256) on all its points.
2. Greedy box: start with every condition free on [0,1]. Each step samples N_BOX
   points uniformly in the current box; for every condition and every
   contiguous interval on a grid of GRID bins inside its current range, the
   points of the current sample that fall in the interval are a uniform sample
   of the narrowed box, so its purity (share with map P > 0.5) is read off
   directly. The chosen constraint maximises purity gain per unit of log-volume
   given up, (purity_new - purity) / -log(width_new / width_old), i.e. it buys
   the most purity for the least narrowing. Repeat until purity >= TARGET.
   The grid has 12 bins so the integer K boundaries (1/3, 2/3) lie on it.
3. Confirm each final box by direct simulation: N_CONFIRM random points inside
   the box x CONFIRM_SEEDS runs; report the share with P > 0.5 and the mean P.
4. Between-cell variation check for the equal_split knockout: variance of the
   gene-0 fraction across functional cells at the end, baseline vs equal_split.
Writes s3b.json.
"""
import argparse, json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "stage2"))
sys.path.insert(2, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "stage1b")))
import bench1b  # noqa: E402
import knockouts as ko  # noqa: E402
import surrogate  # noqa: E402
from validate_knockouts import to_u, REF  # noqa: E402

bench1b.model = ko
bench1b.FAST = True
ko.KNOCKOUT = None
GRID = 12
N_BOX = 50000
FROM_U = {
    "S": lambda u: 8 * 64 ** u, "sigma": lambda u: 1 + 2 * u, "mu": lambda u: 1e-3 * 300 ** u,
    "K": lambda u: min(1 + math.floor(u * 3 - 1e-9), 3) if u > 0 else 1, "d": lambda u: 0.02 * u,
    "food": lambda u: 0.25 * 4 ** u, "km": lambda u: 0.1 * 300 ** u, "dm": lambda u: 0.02 * u,
    "z1": lambda u: u, "z2": lambda u: u,
}


def predict(net, u, dev):
    return surrogate.predict_multi(net, torch.tensor(np.asarray(u, np.float32), device=dev)).cpu().numpy()[0]


def greedy_box(net, dev, target, rng):
    lo, hi = np.zeros(ko.DIM), np.ones(ko.DIM)
    steps = []
    while True:
        U = lo + rng.random((N_BOX, ko.DIM)) * (hi - lo)
        life = predict(net, U, dev) > 0.5
        purity = float(life.mean())
        if purity >= target or len(steps) >= 30:
            return lo, hi, purity, steps
        best = None
        for j in range(ko.DIM):
            g = np.round(np.linspace(0, 1, GRID + 1), 6)
            g = g[(g >= lo[j] - 1e-9) & (g <= hi[j] + 1e-9)]
            for ia in range(len(g)):
                for ib in range(ia + 1, len(g)):
                    a, b = g[ia], g[ib]
                    if a <= lo[j] + 1e-9 and b >= hi[j] - 1e-9:
                        continue
                    m = (U[:, j] >= a) & (U[:, j] <= b)
                    if m.sum() < 200:
                        continue
                    pn = float(life[m].mean())
                    shrink = -math.log((b - a) / (hi[j] - lo[j]))
                    score = (pn - purity) / shrink
                    if best is None or score > best[0]:
                        best = (score, j, a, b, pn)
        if best is None or best[0] <= 0:
            return lo, hi, purity, steps
        _, j, a, b, pn = best
        lo[j], hi[j] = a, b
        steps.append({"condition": ko.NAMES[j], "u_range": [float(a), float(b)],
                      "value_range": [FROM_U[ko.NAMES[j]](a), FROM_U[ko.NAMES[j]](b)], "purity_after": pn})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--budget", type=int, default=32000)
    ap.add_argument("--confirm", type=int, default=2000)
    ap.add_argument("--confirm-seeds", type=int, default=16)
    ap.add_argument("--out", default="s3b.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    t0 = time.time()
    res = {"budget": a.budget, "grid": GRID, "levels": {}}
    for L in (1, 2, 3):
        ko.LEVEL = L
        t = time.time()
        x, s, n, _ = bench1b.m_engine3(bench1b.Sim(dev, 2000 + 10 * L, 1 << 15), a.budget, 0, dev, soft=True)
        net = surrogate.fit_es_multi(torch.tensor(x, device=dev), torch.tensor(s, device=dev),
                                     torch.tensor(n, device=dev), [len(x)], seed=99, hidden=256, max_steps=10000)
        print(f"\nL{L}: map from {len(x)} runs in {time.time()-t:.0f}s", flush=True)
        lvl = {}
        for target in (0.95, 0.90):
            lo, hi, purity, steps = greedy_box(net, dev, target, np.random.default_rng(L))
            vol = float(np.prod(hi - lo))
            # direct confirmation inside the box
            g = torch.Generator(device=dev).manual_seed(500 + L)
            U = lo + np.random.default_rng(600 + L).random((a.confirm, ko.DIM)) * (hi - lo)
            ur = torch.tensor(np.repeat(U, a.confirm_seeds, 0), dtype=torch.float32, device=dev)
            y = torch.cat([ko.simulate(ur[i:i + (1 << 15)], generator=g, level=L)
                           for i in range(0, len(ur), 1 << 15)]).float().view(a.confirm, a.confirm_seeds).mean(1)
            y = y.cpu().numpy()
            free = [nm for j, nm in enumerate(ko.NAMES) if lo[j] == 0 and hi[j] == 1]
            box = {nm: [FROM_U[nm](lo[j]), FROM_U[nm](hi[j])] for j, nm in enumerate(ko.NAMES)
                   if not (lo[j] == 0 and hi[j] == 1)}
            lvl[f"target_{int(target * 100)}"] = {
                "steps": steps, "map_purity": purity, "volume": vol, "box_values": box, "free": free,
                "confirm_share_p_gt_half": float(np.mean(y > 0.5)), "confirm_mean_p": float(y.mean())}
            print(f"  target {target:.2f}: {len(steps)} constraints, volume {vol:.3f}, map purity {purity:.3f}, "
                  f"direct check: share P>0.5 {np.mean(y > 0.5):.3f}, mean P {y.mean():.3f}", flush=True)
            for st in steps:
                v0, v1 = st["value_range"]
                print(f"    {st['condition']:>5} in [{v0:.4g}, {v1:.4g}]  -> purity {st['purity_after']:.3f}")
            print(f"    free: {', '.join(free)}", flush=True)
        res["levels"][f"L{L}"] = lvl
        json.dump(res, open(a.out, "w"))

    # between-cell variance of the gene-0 fraction: baseline vs equal_split
    print("\nBetween-cell variance of gene-0 fraction among functional cells (mean over 200 populations)")
    var = {}
    g = torch.Generator(device=dev).manual_seed(7)
    for label, kw in (("reference", REF), ("S=256 sigma=2", dict(REF, S=256, sigma=2.0))):
        u = torch.tensor([to_u(**kw)] * 200, dtype=torch.float32, device=dev)
        for k in (None, "equal_split"):
            _, st = ko.simulate(u, generator=g, return_stats=True, knockout=k, level=3)
            v = float(st["var_g0"].mean())
            var[f"{label} | {k}"] = {"var_g0": v, "P_L3": float(st["l3"].float().mean()),
                                     "functional": float(st["frac_alive"].mean())}
            print(f"  {label:>14} {str(k):>12}: var {v:.5f} (sd {math.sqrt(v):.4f}), P(L3) {st['l3'].float().mean():.2f}, "
                  f"functional {st['frac_alive'].mean():.2f}", flush=True)
    res["between_cell_variance"] = var
    res["seconds"] = time.time() - t0
    json.dump(res, open(a.out, "w"))
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
