"""Stage-1b benchmark on the 10-condition model.

Changes from stage 1, following its negative result:
  * harder problem (10 conditions, ~11% life, 2 distractors)
  * one shared neural-network map-maker for all methods (surrogate.py)
  * engine3: one run per condition, Sobol coverage, and boundary refinement
    driven by disagreement between bootstrap networks (epistemic uncertainty)
    instead of p(1-p), which chased intrinsically random regions
Each method produces an ordered list of tested conditions; accuracy is
measured on prefixes of that list at fixed numbers of simulation runs.
"""
import argparse, json, time
import numpy as np
import torch
from scipy.stats import qmc

import model10 as model
import surrogate

# map-maker settings; --surrogate es switches both to early-stopping fits (stage 1c)
EVAL_FIT = lambda x, s, n, seed: surrogate.fit(x, s, n, seed=seed)
ACQ_FIT = lambda x, s, n, seed: surrogate.fit(x, s, n, seed=seed, bootstrap=True)

CHECKPOINTS = [2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000]


def pick_device(name):
    if name != "auto":
        return torch.device(name)
    if hasattr(torch, "xpu") and torch.xpu.is_available():
        return torch.device("xpu")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


class Sim:
    """Counts runs; every run draws fresh random numbers from one seeded stream."""

    def __init__(self, device, seed, chunk):
        self.device, self.chunk = device, chunk
        self.gen = torch.Generator(device=device)
        self.gen.manual_seed(seed)

    def __call__(self, x_np):
        out = []
        for i in range(0, len(x_np), self.chunk):
            x = torch.tensor(x_np[i:i + self.chunk], dtype=torch.float32, device=self.device)
            out.append(model.simulate(x, generator=self.gen).cpu().numpy())
        return np.concatenate(out) if out else np.zeros(0, dtype=bool)


def ground_truth(sim, n_test, seeds, seed=12345):
    x = np.random.default_rng(seed).random((n_test, model.DIM)).astype(np.float32)
    y = sim(np.repeat(x, seeds, axis=0)).reshape(n_test, seeds)
    return x, y.mean(1)


# ---- methods: each returns arrays (x [K,D], s [K], n [K]) in testing order ----

def m_points(sim, budget, rep, seeds, sobol):
    k = budget // seeds
    if sobol:
        x = qmc.Sobol(d=model.DIM, scramble=True, seed=1000 + rep).random(k).astype(np.float32)
    else:
        x = np.random.default_rng(1000 + rep).random((k, model.DIM)).astype(np.float32)
    y = sim(np.repeat(x, seeds, axis=0)).reshape(k, seeds)
    return x, y.sum(1).astype(np.float32), np.full(k, seeds, np.float32), {}


def m_engine3(sim, budget, rep, device, explore=0.3, n_init=2048, ens=5, pool_bits=19,
              anneal_from=None, explore_end=0.7, soft=False):
    # stage 1d: anneal_from = runs after which the explore share rises linearly to explore_end;
    # soft = draw refinement points with probability ~ disagreement^2 instead of strict top-k
    rng = np.random.default_rng(3000 + rep)
    sob = qmc.Sobol(d=model.DIM, scramble=True, seed=1000 + rep)
    pool = sob.random(1 << pool_bits).astype(np.float32)
    used = np.zeros(len(pool), bool)
    seen = set()  # rule 1: canonical keys of every condition already simulated
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
        keep = []
        for i in idx:
            key = tuple((pool[i] * 4096).astype(np.int64))
            if key not in seen:
                seen.add(key)
                keep.append(i)
        used[idx] = True
        if keep:
            xb = pool[keep]
            xs.append(xb)
            ss.append(sim(xb).astype(np.float32))

    test(take_seq(n_init))
    rounds = 0
    pool_t = torch.tensor(pool, device=device)
    while sum(len(a) for a in xs) < budget and nxt < len(pool):
        total = sum(len(a) for a in xs)
        batch = min(max(1024, total // 3), budget - total)
        frac = explore
        if anneal_from is not None and total > anneal_from:
            frac = explore + (explore_end - explore) * min(1.0, (total - anneal_from) / max(1, budget - anneal_from))
        n_exp = int(batch * frac)
        idx_seq = take_seq(n_exp)
        used[idx_seq] = True
        x = torch.tensor(np.concatenate(xs), device=device)
        s = torch.tensor(np.concatenate(ss), device=device)
        n = torch.ones_like(s)
        preds = torch.stack([surrogate.predict(ACQ_FIT(x, s, n, rep * 100 + rounds * 10 + e), pool_t)
                             for e in range(ens)])
        score = preds.std(0).cpu().numpy()
        score[used] = -1.0
        m = batch - len(idx_seq)
        if soft:
            w = np.where(score > 0, score, 0.0) ** 2
            idx_act = rng.choice(len(pool), size=min(m, int((w > 0).sum())), replace=False, p=w / w.sum())
        else:
            idx_act = np.argsort(-score)[:m]
        test(np.concatenate([idx_seq, idx_act]))
        rounds += 1
    x = np.concatenate(xs)[:budget]
    s = np.concatenate(ss)[:budget]
    return x, s, np.ones(len(x), np.float32), {"rounds": rounds}


def evaluate(x, s, n, test_x, test_p, hard, device, budget, rep):
    cum = np.cumsum(n)
    tq = torch.tensor(test_x, device=device)
    curve = []
    for c in CHECKPOINTS:
        if c > budget:
            break
        k = int(np.searchsorted(cum, c, side="right"))
        if k < 10:
            continue
        net = EVAL_FIT(torch.tensor(x[:k], device=device), torch.tensor(s[:k], device=device),
                       torch.tensor(n[:k], device=device), 7 + rep)
        p = surrogate.predict(net, tq).cpu().numpy()
        ok = (p > 0.5) == (test_p > 0.5)
        curve.append((int(cum[k - 1]), float(ok.mean()), float(ok[hard].mean())))
    return curve


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="auto")
    ap.add_argument("--budget", type=int, default=128000)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--n-test", type=int, default=8000)
    ap.add_argument("--gt-seeds", type=int, default=128)
    ap.add_argument("--chunk", type=int, default=1 << 20, help="max simulations per GPU batch")
    ap.add_argument("--methods", default="random_1seed,sobol_1seed,sobol_8seed,engine3,engine3_explore50")
    ap.add_argument("--surrogate", choices=["fixed", "es"], default="fixed",
                    help="fixed: 600-step fits (stage 1b); es: early-stopping fits (stage 1c)")
    ap.add_argument("--out", default="results_1b.json")
    a = ap.parse_args()
    if a.surrogate == "es":
        global EVAL_FIT, ACQ_FIT
        EVAL_FIT = lambda x, s, n, seed: surrogate.fit_es(x, s, n, seed=seed, hidden=256, max_steps=10000)
        ACQ_FIT = lambda x, s, n, seed: surrogate.fit_es(x, s, n, seed=seed, hidden=128, max_steps=3000,
                                                         bootstrap=True)
    dev = pick_device(a.device)
    print("device:", dev, flush=True)

    t0 = time.time()
    test_x, test_p = ground_truth(Sim(dev, 12345, a.chunk), a.n_test, a.gt_seeds)
    hard = ((test_p >= 0.1) & (test_p <= 0.4)) | ((test_p >= 0.6) & (test_p <= 0.9))
    print(f"ground truth: {a.n_test} x {a.gt_seeds} runs in {time.time()-t0:.0f}s; life fraction "
          f"{np.mean(test_p > 0.5):.3f}; hard set {hard.sum()}", flush=True)

    methods = {
        "random_1seed": lambda sim, r: m_points(sim, a.budget, r, 1, False),
        "sobol_1seed": lambda sim, r: m_points(sim, a.budget, r, 1, True),
        "sobol_8seed": lambda sim, r: m_points(sim, a.budget, r, 8, True),
        "engine3": lambda sim, r: m_engine3(sim, a.budget, r, dev),
        "engine3_explore50": lambda sim, r: m_engine3(sim, a.budget, r, dev, explore=0.5),
        "engine4_anneal": lambda sim, r: m_engine3(sim, a.budget, r, dev, anneal_from=32000),
        "engine4_soft": lambda sim, r: m_engine3(sim, a.budget, r, dev, soft=True),
        "engine4_both": lambda sim, r: m_engine3(sim, a.budget, r, dev, anneal_from=32000, soft=True),
    }
    res = {"surrogate": a.surrogate, "device": str(dev), "budget": a.budget, "reps": a.reps, "n_test": a.n_test, "gt_seeds": a.gt_seeds,
           "life_fraction": float(np.mean(test_p > 0.5)), "hard": int(hard.sum()), "methods": {}}
    for name in a.methods.split(","):
        res["methods"][name] = []
        for r in range(a.reps):
            t = time.time()
            x, s, n, extra = methods[name](Sim(dev, 2000 + r, a.chunk), r)
            curve = evaluate(x, s, n, test_x, test_p, hard, dev, a.budget, r)
            res["methods"][name].append({"curve": curve, **extra})
            print(f"{name} rep {r}: {curve[-1]} ({time.time()-t:.0f}s) {extra}", flush=True)
            json.dump(res, open(a.out, "w"))

    cps = [c for c in CHECKPOINTS if c <= a.budget]
    lines = []
    for col, label in ((1, "all test conditions"), (2, "hard (stochastic) conditions")):
        lines += [f"\nAccuracy on {label} (mean of {a.reps} reps), by simulation runs",
                  "| method | " + " | ".join(f"{c // 1000}k" for c in cps) + " |", "|---|" + "---|" * len(cps)]
        for name, reps in res["methods"].items():
            cells = []
            for c in cps:
                v = [next((pt[col] for pt in rr["curve"] if abs(pt[0] - c) <= 8), None) for rr in reps]
                cells.append("-" if any(x is None for x in v) else f"{np.mean(v):.3f} ± {np.std(v):.3f}")
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
    res["summary"] = "\n".join(lines)
    print(res["summary"])
    json.dump(res, open(a.out, "w"))


if __name__ == "__main__":
    main()
