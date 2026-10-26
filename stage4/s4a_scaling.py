"""Stage-4a: does COMPASS (engine4_soft, FAST) still beat Sobol when the condition space
has 20-30 dimensions instead of 10?

The 2b protocell model at level L3 (lot-s3/protocell.py, as in s3e_straddle.py) is
wrapped so that DIM = 10 + extra, the extra inputs appended after the 10 original ones.
Two padding types:
  distractor  the extra inputs are ignored by the model (pure noise dimensions)
  weak        each extra input j multiplies food by (1 - 0.15 * (u_j - 0.5) * c_j), with
              fixed random signs c_j = +-1 (rng seed 777). One extra dim moves food by at
              most +-7.5%; their product spreads importance over all extra dims. The
              multiplier is applied exactly by shifting u_5 in log4 space
              (food = 0.25 * 4 ** u_5), which may leave [0, 1], so food can go slightly
              below 0.25 or above 1.
Methods: sobol_1seed (bench1b.m_points) and engine4_soft (m_select 'sd', soft=True,
copied from s3e_straddle.py). Scoring: bench1b.evaluate (FAST). Ground truth: the
first 10 coordinates of the test conditions are the same as stage 2b
(rng(12345).random((n, 10))), the extra ones come from rng(12346); same Sim seed.

  python s4a_scaling.py --configs distractor20,weak20,distractor30,weak30 --reps 3
"""
import argparse, json, math, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lot-s3"))
sys.path.insert(1, os.path.join(HERE, "..", "lot-scm"))
sys.path.insert(2, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "lot-run", "stage1b")))
import bench1b  # noqa: E402
import protocell  # noqa: E402
import surrogate  # noqa: E402

bench1b.FAST = True
protocell.LEVEL = 3
WEAK_AMP = 0.15


class Padded:
    """protocell with `extra` appended inputs; looked up by bench1b as `model`."""

    def __init__(self, kind, extra):
        self.kind, self.extra = kind, extra
        self.DIM = protocell.DIM + extra
        self.signs = np.random.default_rng(777).choice([-1.0, 1.0], size=extra).astype(np.float32)

    def simulate(self, u, generator=None):
        base = u[:, :protocell.DIM].clone()
        if self.kind == "weak" and self.extra > 0:
            c = torch.tensor(self.signs, device=u.device)
            mult = (1.0 - WEAK_AMP * (u[:, protocell.DIM:] - 0.5) * c).prod(1)
            base[:, 5] = base[:, 5] + torch.log(mult) / math.log(4.0)
        return protocell.simulate(base, generator=generator)


model = None   # set per config


def ground_truth(sim, n_test, seeds):
    x0 = np.random.default_rng(12345).random((n_test, protocell.DIM))
    x1 = np.random.default_rng(12346).random((n_test, model.DIM - protocell.DIM))
    x = np.concatenate([x0, x1], 1).astype(np.float32)
    y = sim(np.repeat(x, seeds, axis=0)).reshape(n_test, seeds)
    return x, y.mean(1)


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
    t_fit = 0.0
    pool_t = torch.tensor(pool, device=device)
    while sum(len(v) for v in xs) < budget and nxt < len(pool):
        total = sum(len(v) for v in xs)
        batch = min(max(1024, total // 3), budget - total)
        idx_seq = take_seq(int(batch * explore))
        used[idx_seq] = True
        tf = time.time()
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
        t_fit += time.time() - tf
        test(np.concatenate([idx_seq, idx_act]).astype(int))
        rounds += 1
    x, s = np.concatenate(xs)[:budget], np.concatenate(ss)[:budget]
    return x, s, np.ones(len(x), np.float32), {"rounds": rounds, "select_seconds": t_fit}


METHODS = {
    "sobol_1seed": lambda sim, a, r, dev: bench1b.m_points(sim, a.budget, r, 1, True),
    "engine4_soft": lambda sim, a, r, dev: m_select(sim, a.budget, r, dev, "sd", True),
}


def at(curve, c, col):
    return next((pt[col] for pt in curve if abs(pt[0] - c) <= 8), None)


def reach(cps, means, target):
    """Budget at which the mean curve first reaches target (log-linear interpolation;
    the first checkpoint if it is already there, shown as <=)."""
    for i, (c, v) in enumerate(zip(cps, means)):
        if v >= target - 1e-12:
            if i == 0:
                return c
            c0, v0 = cps[i - 1], means[i - 1]
            f = (target - v0) / (v - v0)
            return float(np.exp(np.log(c0) + f * (np.log(c) - np.log(c0))))
    return None


def summarize(res, budget):
    cps = [c for c in bench1b.CHECKPOINTS if c <= budget]
    lines = []
    for cfg, cr in res["configs"].items():
        lines.append(f"\n## {cfg} (D = {cr['dim']}, hard set {cr['hard']}, life {cr['life']:.3f})")
        for col, label in ((1, "all"), (2, "hard")):
            lines.append(f"\nAccuracy ({label}), mean ± sd over reps")
            lines.append("| method | " + " | ".join(f"{c // 1000}k" for c in cps) + " | reaches sobol@max at |")
            lines.append("|---|" + "---|" * (len(cps) + 1))
            sob = cr["methods"].get("sobol_1seed")
            target = None
            if sob and all(at(rr["curve"], cps[-1], col) is not None for rr in sob):
                target = np.mean([at(rr["curve"], cps[-1], col) for rr in sob])
            for name, reps in cr["methods"].items():
                cells, means = [], []
                for c in cps:
                    v = [at(rr["curve"], c, col) for rr in reps]
                    if any(q is None for q in v):
                        cells.append("-"); means.append(-1.0)
                    else:
                        cells.append(f"{np.mean(v):.3f} ± {np.std(v):.3f}"); means.append(float(np.mean(v)))
                r_ = reach(cps, means, target) if target is not None else None
                lines.append(f"| {name} | " + " | ".join(cells) + f" | {'not reached' if r_ is None else ('<= ' if r_ == cps[0] else '') + f'{r_ / 1000:.1f}k'} |")
        base = np.mean([rr["method_seconds"] for rr in cr["methods"].get("sobol_1seed", [{"method_seconds": 0}])])
        lines.append("\nMethod seconds per rep (simulation + selection, excluding scoring)")
        for name, reps in cr["methods"].items():
            ms = [rr["method_seconds"] for rr in reps]
            sel = [rr.get("select_seconds", 0.0) for rr in reps]
            lines.append(f"  {name}: {np.mean(ms):.0f} ± {np.std(ms):.0f} s; overhead vs sobol "
                         f"{np.mean(ms) - base:.0f} s; surrogate fit+score {np.mean(sel):.0f} s")
    return "\n".join(lines)


def main():
    global model
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--budget", type=int, default=64000)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--n-test", type=int, default=4000)
    ap.add_argument("--gt-seeds", type=int, default=64)
    ap.add_argument("--configs", default="distractor20,weak20,distractor30,weak30")
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--out", default="s4a_scaling.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    t0 = time.time()
    res = {"budget": a.budget, "reps": a.reps, "n_test": a.n_test, "gt_seeds": a.gt_seeds,
           "weak_amp": WEAK_AMP, "configs": {}}
    for cfg in a.configs.split(","):
        kind, dim = cfg.rstrip("0123456789"), int(cfg[len(cfg.rstrip("0123456789")):])
        model = Padded(kind, dim - protocell.DIM)
        bench1b.model = model
        t = time.time()
        test_x, test_p = ground_truth(bench1b.Sim(dev, 12345, 1 << 15), a.n_test, a.gt_seeds)
        hard = ((test_p >= 0.1) & (test_p <= 0.4)) | ((test_p >= 0.6) & (test_p <= 0.9))
        cr = {"kind": kind, "dim": dim, "signs": model.signs.tolist(), "hard": int(hard.sum()),
              "life": float(np.mean(test_p > 0.5)), "gt_seconds": time.time() - t, "methods": {}}
        res["configs"][cfg] = cr
        print(f"[{cfg}] ground truth {a.n_test} x {a.gt_seeds} in {time.time()-t:.0f}s; "
              f"life {cr['life']:.3f}; hard {cr['hard']}", flush=True)
        for name in a.methods.split(","):
            cr["methods"][name] = []
            for r in range(a.reps):
                t = time.time()
                x, s, n, extra = METHODS[name](bench1b.Sim(dev, 2000 + r, 1 << 15), a, r, dev)
                t_method = time.time() - t
                curve = bench1b.evaluate(x, s, n, test_x, test_p, hard, dev, a.budget, r)
                cr["methods"][name].append({"curve": curve, "seconds": time.time() - t,
                                            "method_seconds": t_method, **extra})
                print(f"[{cfg}] {name} rep {r}: {curve[-1]} ({time.time()-t:.0f}s, method {t_method:.0f}s)",
                      flush=True)
                json.dump(res, open(a.out, "w"))
    res["total_seconds"] = time.time() - t0
    json.dump(res, open(a.out, "w"))
    print(summarize(res, a.budget))
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
