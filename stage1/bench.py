"""Stage-1 benchmark: does ordered testing find the life/no-life boundary
with fewer simulations than random search?

Ground truth: held-out random conditions, each simulated GT_SEEDS times.
Every method gets the same surrogate (distance-weighted kNN) to turn its
tested points into a map; we compare accuracy against the number of
stochastic runs it spent.
"""
import argparse, json, time
import numpy as np
import torch
from scipy.stats import qmc

import model
from engine import OrderedEngine, Runner, knn_predict


def pick_device(name):
    if name != "auto":
        return torch.device(name)
    if hasattr(torch, "xpu") and torch.xpu.is_available():
        return torch.device("xpu")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def ground_truth(device, n_test, seeds, seed=12345, chunk=200_000):
    rng = np.random.default_rng(seed)
    x = rng.random((n_test, model.DIM)).astype(np.float32)
    run = Runner(model.simulate, device, seed)
    reps = np.repeat(x, seeds, axis=0)
    ys = [run(reps[i:i + chunk]) for i in range(0, len(reps), chunk)]
    p = np.concatenate(ys).reshape(n_test, seeds).mean(1)
    return x, p


HARD = None  # mask of stochastic-but-decidable test conditions


def accuracy(pred_p, true_p):
    ok = (pred_p > 0.5) == (true_p > 0.5)
    return (float(ok.mean()), float(ok[HARD].mean()))


def run_random(device, test_x, test_p, budget, seeds_per_point, rep, sobol=False):
    n_pts = budget // seeds_per_point
    if sobol:
        x = qmc.Sobol(d=model.DIM, scramble=True, seed=1000 + rep).random(n_pts).astype(np.float32)
    else:
        x = np.random.default_rng(1000 + rep).random((n_pts, model.DIM)).astype(np.float32)
    run = Runner(model.simulate, device, 2000 + rep)
    reps = np.repeat(x, seeds_per_point, axis=0)
    y = np.concatenate([run(reps[i:i + 200_000]) for i in range(0, len(reps), 200_000)])
    s = y.reshape(n_pts, seeds_per_point).sum(1)
    tq = torch.tensor(test_x, device=device)
    curve = []
    m = 256 // seeds_per_point * 4 if seeds_per_point < 64 else 4
    while m <= n_pts:
        tx = torch.tensor(x[:m], device=device)
        tp = torch.tensor((s[:m] + 0.5) / (seeds_per_point + 1.0), dtype=torch.float32, device=device)
        tw = torch.full((m,), float(seeds_per_point), device=device)
        acc = accuracy(knn_predict(tx, tp, tw, tq).cpu().numpy(), test_p)
        curve.append((m * seeds_per_point, acc))
        m = int(m * 1.25) + 1
    return curve, {}


def run_engine(device, test_x, test_p, budget, rep, **kw):
    eng = OrderedEngine(model.simulate, model.DIM, model.MONOTONE, device, seed=rep, **kw)
    tq = torch.tensor(test_x, device=device)
    curve = []
    while eng.run.runs < budget:
        eng.step()
        curve.append((eng.run.runs, accuracy(eng.predict(tq).cpu().numpy(), test_p)))
    extra = dict(points_tested=len(eng.archive.x), inferred=eng.inferred_count,
                 spot_checks=eng.spot_checks, violations=eng.violations)
    return curve, extra


METHODS = {
    "random_1seed": lambda d, tx, tp, b, r: run_random(d, tx, tp, b, 1, r),
    "random_8seed": lambda d, tx, tp, b, r: run_random(d, tx, tp, b, 8, r),
    "sobol_8seed": lambda d, tx, tp, b, r: run_random(d, tx, tp, b, 8, r, sobol=True),
    "sobol_sprt": lambda d, tx, tp, b, r: run_engine(d, tx, tp, b, r, use_active=False, use_monotone=False),
    "engine_no_monotone": lambda d, tx, tp, b, r: run_engine(d, tx, tp, b, r, use_monotone=False),
    "engine_full": lambda d, tx, tp, b, r: run_engine(d, tx, tp, b, r),
    "sobol_1seed": lambda d, tx, tp, b, r: run_random(d, tx, tp, b, 1, r, sobol=True),
    "engine2_no_monotone": lambda d, tx, tp, b, r: run_engine(d, tx, tp, b, r, use_sprt=False, fixed_seeds=1, use_monotone=False, n_init=1024, batch=512),
    "engine2_full": lambda d, tx, tp, b, r: run_engine(d, tx, tp, b, r, use_sprt=False, fixed_seeds=1, n_init=1024, batch=512),
}


def runs_to_reach(curve, target, which=0):
    for runs, acc in curve:
        if acc[which] >= target:
            return runs
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="auto")
    ap.add_argument("--budget", type=int, default=64_000)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--n-test", type=int, default=2000)
    ap.add_argument("--gt-seeds", type=int, default=128)
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--out", default="results.json")
    a = ap.parse_args()
    dev = pick_device(a.device)
    print("device:", dev, flush=True)

    t0 = time.time()
    test_x, test_p = ground_truth(dev, a.n_test, a.gt_seeds)
    global HARD
    HARD = ((test_p >= 0.1) & (test_p <= 0.4)) | ((test_p >= 0.6) & (test_p <= 0.9))
    print(f"hard set (0.1<=p<=0.4 or 0.6<=p<=0.9): {HARD.sum()} conditions", flush=True)
    print(f"ground truth: {a.n_test} conditions x {a.gt_seeds} seeds in {time.time()-t0:.0f}s; "
          f"life fraction {np.mean(test_p > 0.5):.3f}", flush=True)

    results = {"device": str(dev), "budget": a.budget, "reps": a.reps, "n_test": a.n_test,
               "gt_seeds": a.gt_seeds, "life_fraction": float(np.mean(test_p > 0.5)), "methods": {}}
    for name in a.methods.split(","):
        runs = []
        for r in range(a.reps):
            t = time.time()
            curve, extra = METHODS[name](dev, test_x, test_p, a.budget, r)
            print(f"{name} rep {r}: final acc {curve[-1][1][0]:.4f} / hard {curve[-1][1][1]:.4f} at {curve[-1][0]} runs "
                  f"({time.time()-t:.0f}s) {extra}", flush=True)
            runs.append({"curve": curve, **extra})
        results["methods"][name] = runs
        json.dump(results, open(a.out, "w"))

    # summary at matched budgets (mean over reps of accuracy interpolated at checkpoints)
    def acc_at(curve, runs, which):
        best = None
        for r_, acc in curve:
            if r_ <= runs:
                best = acc[which]
        return best
    cps = [c for c in (2000, 4000, 8000, 16000, 32000, 64000, 128000) if c <= a.budget]
    lines = []
    for which, label in ((0, "all test conditions"), (1, "hard (stochastic) conditions")):
        lines.append(f"\nAccuracy on {label}, by number of simulation runs spent")
        lines.append("| method | " + " | ".join(f"{c:,}" for c in cps) + " |")
        lines.append("|---|" + "---|" * len(cps))
        for name, reps in results["methods"].items():
            cells = []
            for c in cps:
                v = [acc_at(r["curve"], c, which) for r in reps]
                cells.append("-" if any(x is None for x in v) else f"{np.mean(v):.3f}")
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
    summary = "\n".join(lines)
    print(summary)
    results["summary"] = summary
    json.dump(results, open(a.out, "w"))


if __name__ == "__main__":
    main()
