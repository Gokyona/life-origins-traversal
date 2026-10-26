"""Stage-3e(b): straddle level-set selection versus COMPASS disagreement selection.

Same engine loop as bench1b.m_engine3 (FAST map-maker, 5-member bootstrap ensemble,
2048 Sobol start, rounds growing by 1/3, 30% Sobol exploration per round). Only the
rule that picks the other 70% of each round differs:
  engine3        top-k of the ensemble sd                         (stage 1b)
  engine4_soft   sample with probability ~ sd^2                    (stage 1d default, COMPASS)
  straddle       top-k of 1.96 sd - |mean - 0.5|                   (Bryan et al. 2005)
  straddle_soft  sample with probability ~ max(0, straddle score)^2
plus sobol_1seed as the passive baseline. Benchmark: the 2b protocell model at
level L3 (same ground truth settings and seeds as bench_protocell.py), scored by
bench1b.evaluate.
"""
import argparse, json, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "lot-scm"))
sys.path.insert(2, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "lot-run", "stage1b")))
import bench1b  # noqa: E402
import protocell as model  # noqa: E402
import surrogate  # noqa: E402

bench1b.model = model
bench1b.FAST = True
model.LEVEL = 3


def m_select(sim, budget, rep, device, rule, soft, explore=0.3, n_init=2048, ens=5, pool_bits=19):
    rng = np.random.default_rng(3000 + rep)
    pool = qmc.Sobol(d=model.DIM, scramble=True, seed=1000 + rep).random(1 << pool_bits).astype(np.float32)
    used = np.zeros(len(pool), bool)
    nxt = 0
    xs, ss = [], []

    def take_seq(m):
        nonlocal nxt
        idx = []
        while len(idx) < m and nxt < len(pool):
            if not used[nxt]:
                idx.append(nxt)
            nxt += 1
        return np.array(idx, int)

    def test(idx):
        used[idx] = True
        xs.append(pool[idx])
        ss.append(sim(pool[idx]).astype(np.float32))

    test(take_seq(n_init))
    rounds = 0
    pool_t = torch.tensor(pool, device=device)
    while sum(len(v) for v in xs) < budget and nxt < len(pool):
        total = sum(len(v) for v in xs)
        batch = min(max(1024, total // 3), budget - total)
        idx_seq = take_seq(int(batch * explore))
        used[idx_seq] = True
        x = torch.tensor(np.concatenate(xs), device=device)
        s = torch.tensor(np.concatenate(ss), device=device)
        preds = surrogate.predict_multi(surrogate.fit_es_multi(
            x, s, torch.ones_like(s), [len(x)] * ens, seed=rep * 100 + rounds * 10, hidden=128, max_steps=3000,
            bootstrap=True), pool_t)
        sd = preds.std(0).cpu().numpy()
        if rule == "sd":
            score = sd
        else:   # straddle
            score = 1.96 * sd - np.abs(preds.mean(0).cpu().numpy() - 0.5)
        score[used] = -np.inf
        m = batch - len(idx_seq)
        if soft:
            w = np.where(np.isfinite(score) & (score > 0), score, 0.0) ** 2
            k = min(m, int((w > 0).sum()))
            idx_act = rng.choice(len(pool), size=k, replace=False, p=w / w.sum()) if k > 0 else np.zeros(0, int)
            if k < m:   # not enough positive scores: fill with the best remaining
                rest = np.argsort(-np.where(w > 0, -np.inf, score))[:m - k]
                idx_act = np.concatenate([idx_act, rest])
        else:
            idx_act = np.argsort(-score)[:m]
        test(np.concatenate([idx_seq, idx_act]).astype(int))
        rounds += 1
    x, s = np.concatenate(xs)[:budget], np.concatenate(ss)[:budget]
    return x, s, np.ones(len(x), np.float32), {"rounds": rounds}


METHODS = {
    "sobol_1seed": lambda sim, a, r, dev: bench1b.m_points(sim, a.budget, r, 1, True),
    "engine3": lambda sim, a, r, dev: m_select(sim, a.budget, r, dev, "sd", False),
    "engine4_soft": lambda sim, a, r, dev: m_select(sim, a.budget, r, dev, "sd", True),
    "straddle": lambda sim, a, r, dev: m_select(sim, a.budget, r, dev, "straddle", False),
    "straddle_soft": lambda sim, a, r, dev: m_select(sim, a.budget, r, dev, "straddle", True),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--budget", type=int, default=64000)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--n-test", type=int, default=4000)
    ap.add_argument("--gt-seeds", type=int, default=64)
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--out", default="s3e_straddle.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    t0 = time.time()
    test_x, test_p = bench1b.ground_truth(bench1b.Sim(dev, 12345, 1 << 15), a.n_test, a.gt_seeds)
    hard = ((test_p >= 0.1) & (test_p <= 0.4)) | ((test_p >= 0.6) & (test_p <= 0.9))
    print(f"ground truth {a.n_test} x {a.gt_seeds} in {time.time()-t0:.0f}s; life {np.mean(test_p > 0.5):.3f}; "
          f"hard {hard.sum()}", flush=True)
    res = {"budget": a.budget, "reps": a.reps, "methods": {}}
    for name in a.methods.split(","):
        res["methods"][name] = []
        for r in range(a.reps):
            t = time.time()
            x, s, n, extra = METHODS[name](bench1b.Sim(dev, 2000 + r, 1 << 15), a, r, dev)
            t_method = time.time() - t
            curve = bench1b.evaluate(x, s, n, test_x, test_p, hard, dev, a.budget, r)
            res["methods"][name].append({"curve": curve, "seconds": time.time() - t, "method_seconds": t_method, **extra})
            print(f"{name} rep {r}: {curve[-1]} ({time.time()-t:.0f}s, method {t_method:.0f}s)", flush=True)
            json.dump(res, open(a.out, "w"))
    cps = [c for c in bench1b.CHECKPOINTS if c <= a.budget]
    for col, label in ((1, "all"), (2, "hard")):
        print(f"\nAccuracy ({label}), mean ± sd over {a.reps} reps")
        print("| method | " + " | ".join(f"{c // 1000}k" for c in cps) + " |")
        print("|---|" + "---|" * len(cps))
        for name, reps in res["methods"].items():
            cells = []
            for c in cps:
                v = [next((pt[col] for pt in rr["curve"] if abs(pt[0] - c) <= 8), None) for rr in reps]
                cells.append("-" if any(q is None for q in v) else f"{np.mean(v):.3f} ± {np.std(v):.3f}")
            print(f"| {name} | " + " | ".join(cells) + " |")
    print("\nMethod seconds per rep (simulation + selection, excluding scoring); engine overhead = method - sobol method")
    base = np.mean([rr["method_seconds"] for rr in res["methods"].get("sobol_1seed", [{"method_seconds": 0}])])
    for name, reps in res["methods"].items():
        ms = [rr["method_seconds"] for rr in reps]
        print(f"  {name}: {np.mean(ms):.0f} ± {np.std(ms):.0f} s; overhead vs sobol {np.mean(ms) - base:.0f} s")
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
