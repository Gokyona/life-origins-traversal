"""Does the engine's map recover the known RAF catalysis threshold?

For each slice (max length n = 5..10, food length t = 1..3, all reactions present,
every molecule a potential catalyst, distractors 0.5) find f*, the catalysis level
where P(RAF) = 0.5:
  * direct: brute-force simulation on a fine log-f grid (reference),
  * map:    the shared map-maker (fit_es_multi, hidden 256) trained on a method's
            first B simulations,
for engine4_soft and sobol_1seed at several budgets. Literature: f* about 1-2
for food t=2 (Hordijk & Steel 2004; Hordijk, Kauffman & Steel 2011).
"""
import argparse, json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "stage1b")))
import bench1b  # noqa: E402
import raf  # noqa: E402
import surrogate  # noqa: E402

bench1b.model = raf
bench1b.FAST = True

NS = list(range(raf.N_MIN, raf.N_MAX + 1))
TS = [1, 2, 3]
FGRID = np.exp(np.linspace(math.log(0.5), math.log(10.0), 49))


def grid_u(nv, tv, fs, dev):
    u = torch.zeros(len(fs), raf.DIM, device=dev)
    u[:, 0] = (nv - raf.N_MIN + 0.5) / (raf.N_MAX - raf.N_MIN + 1)
    u[:, 1] = torch.tensor((np.log10(fs) + 1) / 2, dtype=torch.float32, device=dev)
    u[:, 2] = (tv - raf.T_MIN + 0.5) / (raf.T_MAX - raf.T_MIN + 1)
    u[:, 3] = 1.0
    u[:, 4] = 1.0
    u[:, 5:] = 0.5
    return u


def crossing(fs, p):
    """f where p first rises through 0.5 (log-linear interpolation); nan if never."""
    for i in range(1, len(fs)):
        if p[i - 1] < 0.5 <= p[i]:
            w = (0.5 - p[i - 1]) / (p[i] - p[i - 1])
            return float(math.exp(math.log(fs[i - 1]) + w * (math.log(fs[i]) - math.log(fs[i - 1]))))
    return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--budget", type=int, default=32000)
    ap.add_argument("--direct-seeds", type=int, default=400)
    ap.add_argument("--out", default="threshold_check.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    budgets = [b for b in (8000, 16000, 32000, 64000) if b <= a.budget]
    slices = [(nv, tv) for tv in TS for nv in NS]
    direct, maps = {}, {}

    t0 = time.time()
    g = torch.Generator(device=dev).manual_seed(777)
    se = 1.0 / math.sqrt(a.direct_seeds)  # binomial s.e. at p = 0.5
    for sl in slices:
        u = grid_u(*sl, FGRID, dev).repeat_interleave(a.direct_seeds, 0)
        p = raf.simulate(u, generator=g).float().view(len(FGRID), a.direct_seeds).mean(1).cpu().numpy()
        # rough 95% interval for f*: where the curve crosses 0.5 -/+ 2 s.e.
        direct[sl] = {"p": p.tolist(), "f_star": crossing(FGRID, p),
                      "ci": (crossing(FGRID, p + 2 * se), crossing(FGRID, p - 2 * se))}
    print(f"direct simulation done in {time.time()-t0:.0f}s", flush=True)

    methods = {
        "engine4_soft": lambda sim, r: bench1b.m_engine3(sim, a.budget, r, dev, soft=True),
        "sobol_1seed": lambda sim, r: bench1b.m_points(sim, a.budget, r, 1, True),
    }
    for name, fn in methods.items():
        maps[name] = []
        for r in range(a.reps):
            t = time.time()
            x, s, n, _ = fn(bench1b.Sim(dev, 2000 + r, 1 << 20), r)
            cum = np.cumsum(n)
            ks = [int(np.searchsorted(cum, b, side="right")) for b in budgets]
            ens = surrogate.fit_es_multi(torch.tensor(x, device=dev), torch.tensor(s, device=dev),
                                         torch.tensor(n, device=dev), ks, seed=7 + r, hidden=256,
                                         max_steps=10000)
            entry = {}
            for sl in slices:
                p = surrogate.predict_multi(ens, grid_u(*sl, FGRID, dev)).cpu().numpy()
                entry[sl] = {b: crossing(FGRID, p[i]) for i, b in enumerate(budgets)}
            maps[name].append(entry)
            print(f"{name} rep {r}: {time.time()-t:.0f}s", flush=True)

    print("\nf* where P(RAF)=0.5  (dens=1, c=1, distractors 0.5); map columns: mean ± sd over reps")
    print("| n | t | direct sim (95% CI) | "
          + " | ".join(f"{m} {b // 1000}k" for m in methods for b in budgets) + " |")
    print("|---|---|---|" + "---|" * (len(methods) * len(budgets)))
    for sl in slices:
        d = direct[sl]
        cells = [f"{d['f_star']:.3f} ({d['ci'][0]:.3f}-{d['ci'][1]:.3f})"]
        for m in methods:
            for b in budgets:
                v = [e[sl][b] for e in maps[m]]
                cells.append(f"{np.nanmean(v):.3f} ± {np.nanstd(v):.3f}")
        print(f"| {sl[0]} | {sl[1]} | " + " | ".join(cells) + " |")
    print()
    for m in methods:
        for b in budgets:
            err = np.array([abs(math.log10(e[sl][b] / direct[sl]["f_star"]))
                            for e in maps[m] for sl in slices])
            miss = int(np.isnan(err).sum())
            note = f"; {miss} slice-reps without a crossing" if miss else ""
            print(f"{m} {b // 1000}k: mean |log10 error| of f* {np.nanmean(err):.4f} "
                  f"(= {100 * (10 ** np.nanmean(err) - 1):.1f}% in f){note}")
    print("\nreference P(RAF) curves, t=2 (f:p)")
    for nv in (6, 8, 10):
        p = direct[(nv, 2)]["p"]
        print(f"n={nv}: " + ", ".join(f"{f:.2f}:{q:.2f}" for f, q in zip(FGRID[::3], p[::3])))
    key = lambda sl: f"n{sl[0]}_t{sl[1]}"
    json.dump({"f_grid": FGRID.tolist(), "budgets": budgets,
               "direct": {key(k): v for k, v in direct.items()},
               "map": {m: [{key(k): v for k, v in e.items()} for e in reps] for m, reps in maps.items()}},
              open(a.out, "w"))


if __name__ == "__main__":
    main()
