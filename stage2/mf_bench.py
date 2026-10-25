"""Stage-2d benchmark: multi-fidelity COMPASS (engine5_mf) on the protocell model.

Budgets are counted in high-fidelity-equivalent cost: one reference run
(fidelity 1) costs 1, one cheap screening run (fidelity 0) costs C_LO
(measured cost ratio, fidelity_check.py). Methods:
  sobol_1seed   high-fidelity Sobol points
  engine4_soft  COMPASS on high fidelity only (stage 1d default)
  engine5_mf    multi-fidelity COMPASS: every round spends LO_SHARE of its cost
                on low-fidelity Sobol screening and the rest on high-fidelity runs
                (30% Sobol exploration, 70% soft selection on the disagreement of
                the bootstrap ensemble, evaluated at fidelity 1). The map-maker
                gets the fidelity as an extra input, so cheap runs shape the map
                and expensive runs calibrate it.
Stage-2d follow-up variants (same function, flags):
  engine5_decay   the low-fidelity share decays linearly from LO_START (first
                  round) to 0 at half the budget, so late, boundary-deciding
                  rounds are all high fidelity
  engine5_screen  low-fidelity runs only rule regions OUT: low-fi "life" results
                  are discarded (the cheap model is optimistic near the
                  boundary; its "no life" verdicts were 99% reliable in
                  fidelity_check.py), and high-fidelity exploration points that
                  the current map calls clearly lifeless (mean p < SCREEN_P and
                  ensemble sd < SCREEN_P) are skipped
  engine5_both    both
Discarded runs still count towards the cost.
All methods are scored by the same map-maker (fit_es_multi, hidden 256, inputs
x + fidelity flag), predicted at fidelity 1, against high-fidelity ground truth.
Spreading the cheap screen over every round (instead of front-loading it) keeps
every prefix of the run a valid mixed-fidelity design, so accuracy can be read at
any checkpoint.
"""
import argparse, json, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "stage1b")))
import bench1b  # noqa: E402
import protocell as model  # noqa: E402
import surrogate  # noqa: E402

bench1b.model = model
bench1b.FAST = True
C_LO = 1 / 7.8          # cost of a fidelity-0 run in fidelity-1 units (measured)
LO_SHARE = 0.3          # share of each round's cost spent on cheap screening
LO_START = 0.5          # engine5_decay: low-fi share in the first round
SCREEN_P = 0.02         # engine5_screen: "clearly lifeless" threshold
CHECKPOINTS = [2000, 4000, 8000, 16000, 32000, 64000, 128000]


class Sim:
    """Seeded simulator at a fixed fidelity (one random stream per method rep)."""

    def __init__(self, device, seed, fidelity, chunk=1 << 15):
        self.device, self.fidelity, self.chunk = device, fidelity, chunk
        self.gen = torch.Generator(device=device).manual_seed(seed)

    def __call__(self, x_np):
        out = [model.simulate(torch.tensor(x_np[i:i + self.chunk], device=self.device), generator=self.gen,
                              fidelity=self.fidelity).cpu().numpy() for i in range(0, len(x_np), self.chunk)]
        return np.concatenate(out) if out else np.zeros(0, dtype=bool)


def with_flag(x, f):
    return np.concatenate([x, np.full((len(x), 1), f, np.float32)], 1)


def m_sobol(sim_hi, budget, rep):
    x = qmc.Sobol(d=model.DIM, scramble=True, seed=1000 + rep).random(budget).astype(np.float32)
    return with_flag(x, 1.0), sim_hi(x).astype(np.float32), np.ones(budget, np.float32), np.ones(budget, bool), {}


def m_engine4(sim_hi, budget, rep, dev):
    x, s, n, extra = bench1b.m_engine3(sim_hi, budget, rep, dev, soft=True)
    return with_flag(x, 1.0), s, np.ones(len(x), np.float32), np.ones(len(x), bool), extra


def m_engine5(sim_hi, sim_lo, budget, rep, dev, explore=0.3, n_init=2048, ens=5, pool_bits=19,
              decay=False, screen=False):
    rng = np.random.default_rng(3000 + rep)
    pool = qmc.Sobol(d=model.DIM, scramble=True, seed=1000 + rep).random(1 << pool_bits).astype(np.float32)
    lo_seq = qmc.Sobol(d=model.DIM, scramble=True, seed=5000 + rep)
    used = np.zeros(len(pool), bool)
    nxt = 0
    xs, ss, cs, ks_ = [], [], [], []   # x with fidelity flag, successes, cost, keep-for-map per row
    share = lambda spent_: (LO_START * max(0.0, 1.0 - spent_ / (budget / 2)) if decay else LO_SHARE)
    skipped = 0

    def take_seq(m):
        nonlocal nxt
        idx = []
        while len(idx) < m and nxt < len(pool):
            if not used[nxt]:
                idx.append(nxt)
            nxt += 1
        return np.array(idx, int)

    def run_lo(cost):
        k = int(round(cost / C_LO))
        if k > 0:
            x = lo_seq.random(k).astype(np.float32)
            y = sim_lo(x).astype(np.float32)
            xs.append(with_flag(x, 0.0)); ss.append(y); cs.append(np.full(k, C_LO))
            ks_.append(y == 0 if screen else np.ones(k, bool))

    def run_hi(idx):
        if len(idx):
            used[idx] = True
            xs.append(with_flag(pool[idx], 1.0)); ss.append(sim_hi(pool[idx]).astype(np.float32))
            cs.append(np.ones(len(idx))); ks_.append(np.ones(len(idx), bool))

    spent = lambda: float(sum(c.sum() for c in cs))
    f0 = share(0.0)
    run_lo(f0 * n_init)
    run_hi(take_seq(int(n_init * (1 - f0))))
    rounds = 0
    pool_t = torch.tensor(with_flag(pool, 1.0), device=dev)
    while spent() < budget - 1 and nxt < len(pool):
        total = spent()
        batch = min(max(1024, total / 3), budget - total)
        f = share(total)
        run_lo(f * batch)
        m_hi = int(batch * (1 - f))
        keep = np.concatenate(ks_)
        x = torch.tensor(np.concatenate(xs)[keep], device=dev)
        s = torch.tensor(np.concatenate(ss)[keep], device=dev)
        n = torch.ones_like(s)
        preds = surrogate.predict_multi(surrogate.fit_es_multi(
            x, s, n, [len(x)] * ens, seed=rep * 100 + rounds * 10, hidden=128, max_steps=3000,
            bootstrap=True), pool_t)
        mean_p = preds.mean(0).cpu().numpy()
        score = preds.std(0).cpu().numpy()
        n_exp = int(m_hi * explore)
        if screen:
            lifeless = (mean_p < SCREEN_P) & (score < SCREEN_P)
            idx_seq = []
            while len(idx_seq) < n_exp and nxt < len(pool):
                cand = take_seq(n_exp - len(idx_seq))
                ok = ~lifeless[cand]
                skipped += int((~ok).sum())
                used[cand[~ok]] = True     # never revisit a screened-out point
                idx_seq.extend(cand[ok].tolist())
            idx_seq = np.array(idx_seq, int)
        else:
            idx_seq = take_seq(n_exp)
        used[idx_seq] = True
        score[used] = -1.0
        w = np.where(score > 0, score, 0.0) ** 2
        k = min(m_hi - len(idx_seq), int((w > 0).sum()))
        idx_act = rng.choice(len(pool), size=k, replace=False, p=w / w.sum()) if k > 0 else np.zeros(0, int)
        used[idx_seq] = False   # run_hi marks them again
        run_hi(np.concatenate([idx_seq, idx_act]).astype(int))
        rounds += 1
    x, s, c, keep = np.concatenate(xs), np.concatenate(ss), np.concatenate(cs), np.concatenate(ks_)
    n_lo = int((x[:, -1] == 0).sum())
    return x, s, c, keep, {"rounds": rounds, "n_lo": n_lo, "n_hi": int(len(x) - n_lo),
                           "lo_kept": int(keep[x[:, -1] == 0].sum()), "hi_explore_skipped": skipped}


def evaluate(x, s, cost, keep, test_x, test_p, hard, dev, budget, rep):
    cum = np.cumsum(cost)[keep]     # discarded runs still count towards the cost
    x, s = x[keep], s[keep]
    ks, cps = [], []
    for c in CHECKPOINTS:
        if c > budget + 1:
            break
        k = int(np.searchsorted(cum, c + 1e-6, side="right"))
        if k >= 10:
            ks.append(k)
            cps.append(c)
    ens = surrogate.fit_es_multi(torch.tensor(x, device=dev), torch.tensor(s, device=dev),
                                 torch.ones(len(x), device=dev), ks, seed=7 + rep, hidden=256, max_steps=10000)
    p_all = surrogate.predict_multi(ens, torch.tensor(with_flag(test_x, 1.0), device=dev)).cpu().numpy()
    curve = []
    for c, k, p in zip(cps, ks, p_all):
        ok = (p > 0.5) == (test_p > 0.5)
        curve.append((c, float(ok.mean()), float(ok[hard].mean()), int((x[:k, -1] == 0).sum()), int((x[:k, -1] == 1).sum())))
    return curve


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--budget", type=int, default=64000)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--n-test", type=int, default=4000)
    ap.add_argument("--gt-seeds", type=int, default=64)
    ap.add_argument("--methods", default="sobol_1seed,engine4_soft,engine5_mf")
    ap.add_argument("--out", default="results_mf.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    print("device:", dev, f"C_LO={C_LO:.4f} LO_SHARE={LO_SHARE}", flush=True)

    t0 = time.time()
    test_x, test_p = bench1b.ground_truth(bench1b.Sim(dev, 12345, 1 << 15), a.n_test, a.gt_seeds)
    hard = ((test_p >= 0.1) & (test_p <= 0.4)) | ((test_p >= 0.6) & (test_p <= 0.9))
    print(f"ground truth: {a.n_test} x {a.gt_seeds} runs in {time.time()-t0:.0f}s; life fraction "
          f"{np.mean(test_p > 0.5):.3f}; hard set {hard.sum()}", flush=True)

    methods = {
        "sobol_1seed": lambda r: m_sobol(Sim(dev, 2000 + r, 1), a.budget, r),
        "engine4_soft": lambda r: m_engine4(Sim(dev, 2000 + r, 1), a.budget, r, dev),
        "engine5_mf": lambda r: m_engine5(Sim(dev, 2000 + r, 1), Sim(dev, 4000 + r, 0), a.budget, r, dev),
        "engine5_decay": lambda r: m_engine5(Sim(dev, 2000 + r, 1), Sim(dev, 4000 + r, 0), a.budget, r, dev,
                                             decay=True),
        "engine5_screen": lambda r: m_engine5(Sim(dev, 2000 + r, 1), Sim(dev, 4000 + r, 0), a.budget, r, dev,
                                              screen=True),
        "engine5_both": lambda r: m_engine5(Sim(dev, 2000 + r, 1), Sim(dev, 4000 + r, 0), a.budget, r, dev,
                                            decay=True, screen=True),
    }
    res = {"budget": a.budget, "reps": a.reps, "n_test": a.n_test, "gt_seeds": a.gt_seeds, "c_lo": C_LO,
           "lo_share": LO_SHARE, "life_fraction": float(np.mean(test_p > 0.5)), "hard": int(hard.sum()),
           "methods": {}}
    for name in a.methods.split(","):
        res["methods"][name] = []
        for r in range(a.reps):
            t = time.time()
            x, s, c, keep, extra = methods[name](r)
            sim_t = time.time() - t
            curve = evaluate(x, s, c, keep, test_x, test_p, hard, dev, a.budget, r)
            res["methods"][name].append({"curve": curve, "seconds": time.time() - t, **extra})
            print(f"{name} rep {r}: {curve[-1]} ({time.time()-t:.0f}s, method {sim_t:.0f}s) {extra}", flush=True)
            json.dump(res, open(a.out, "w"))

    cps = [c for c in CHECKPOINTS if c <= a.budget]
    lines = []
    for col, label in ((1, "all test conditions"), (2, "hard (stochastic) conditions")):
        lines += [f"\nAccuracy on {label} (mean of {a.reps} reps), by high-fidelity-equivalent cost",
                  "| method | " + " | ".join(f"{c // 1000}k" for c in cps) + " |", "|---|" + "---|" * len(cps)]
        for name, reps in res["methods"].items():
            cells = []
            for c in cps:
                v = [next((pt[col] for pt in rr["curve"] if pt[0] == c), None) for rr in reps]
                cells.append("-" if any(q is None for q in v) else f"{np.mean(v):.3f} ± {np.std(v):.3f}")
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
    lines.append("\nRuns used at each checkpoint (mean over reps): low-fi / high-fi")
    for name, reps in res["methods"].items():
        cells = []
        for c in cps:
            pts = [next((pt for pt in rr["curve"] if pt[0] == c), None) for rr in reps]
            cells.append("-" if any(q is None for q in pts) else
                         f"{np.mean([q[3] for q in pts]):.0f}/{np.mean([q[4] for q in pts]):.0f}")
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    res["summary"] = "\n".join(lines)
    print(res["summary"])
    json.dump(res, open(a.out, "w"))


if __name__ == "__main__":
    main()
