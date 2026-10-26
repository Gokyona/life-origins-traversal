"""Stage-3e(a): a joint nested map for L1..L4 versus separate per-level maps (evolve.py, the 3c model).

Joint map: one network with four outputs z1..z4 and a product link
  P(>=L1) = s(z1),  P(>=L2) = P(>=L1) s(z2),  P(>=L3) = P(>=L2) s(z3),  P(>=L4) = P(>=L3) s(z4)
(s = logistic), so P(>=L4) <= P(>=L3) <= P(>=L2) <= P(>=L1) everywhere by construction.
Trained on the four labels of every simulation at once (each run reports L1..L4),
loss = sum of the four binary cross-entropies, computed in log space; early stopping
on a 10% held-out split as in surrogate.fit_es.
Separate maps: four independent fit_es_multi networks (hidden 256), one per level,
on the same runs.

Training data: the first B points of a Sobol sequence, one run each (B = 8k, 16k, 32k).
Extra comparison (--compass): the points COMPASS (engine4_soft) chooses when it
targets L4 with a 32k budget; those points are then simulated once more with
return_stats so every point carries all four labels (fresh runs), and both map
types are fitted on them.
Scored on the 3c ground truth (same 4000 test conditions x 64 runs as maps_l4.py):
accuracy and hard accuracy per level, and nesting violations of the predicted labels
(share of test points with label(L+1) = 1 and label(L) = 0).
"""
import argparse, json, os, sys, time
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "lot-scm"))
sys.path.insert(2, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "lot-run", "stage1b")))
import evolve as model  # noqa: E402
import bench1b  # noqa: E402

bench1b.model = model
bench1b.FAST = True
import surrogate  # noqa: E402

LEVELS = (1, 2, 3, 4)


class Joint(nn.Module):
    def __init__(self, dim, hidden=256, levels=4):
        super().__init__()
        self.f = nn.Sequential(nn.Linear(dim, hidden), nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU(),
                               nn.Linear(hidden, levels))

    def log_p(self, x):
        """log P(>= L) for L = 1..4, shape [B, 4]."""
        return torch.cumsum(F.logsigmoid(self.f(x)), dim=1)


def bce_log(logp, y):
    log1m = torch.log1p(-torch.exp(logp).clamp(max=1 - 1e-6))
    return -(y * logp + (1 - y) * log1m).mean()


def fit_joint(x, y, seed=0, max_steps=10000, batch=4096, lr=3e-3, val_frac=0.1, eval_every=200, patience=5):
    g = torch.Generator(device="cpu").manual_seed(seed)
    torch.manual_seed(seed)
    N = len(x)
    perm = torch.randperm(N, generator=g)
    n_val = max(1, int(N * val_frac))
    val, tr = perm[:n_val].to(x.device), perm[n_val:]
    net = Joint(x.shape[1]).to(x.device)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    best, best_state, bad = float("inf"), None, 0
    for step in range(1, max_steps + 1):
        idx = tr[torch.randint(0, len(tr), (min(batch, len(tr)),), generator=g)].to(x.device)
        loss = bce_log(net.log_p(x[idx]), y[idx])
        opt.zero_grad()
        loss.backward()
        opt.step()
        if step % eval_every == 0:
            with torch.no_grad():
                vl = bce_log(net.log_p(x[val]), y[val]).item()
            if vl < best - 1e-5:
                best, bad = vl, 0
                best_state = {k: v.clone() for k, v in net.state_dict().items()}
            else:
                bad += 1
                if bad >= patience:
                    break
    net.load_state_dict(best_state)
    net.eval()
    return net


def simulate_all(u, dev, seed, chunk=1 << 15):
    g = torch.Generator(device=dev).manual_seed(seed)
    out = {L: [] for L in LEVELS}
    for i in range(0, len(u), chunk):
        _, st = model.simulate(torch.tensor(u[i:i + chunk], device=dev), generator=g, return_stats=True)
        for L in LEVELS:
            out[L].append(st[f"l{L}"].float().cpu().numpy())
    return np.stack([np.concatenate(out[L]) for L in LEVELS], 1)          # [N, 4]


def score(pred, gt):
    """pred, gt: [N, 4] probabilities. Accuracy / hard accuracy per level, nesting violations."""
    res = {}
    for k, L in enumerate(LEVELS):
        lab, pl = gt[:, k] > 0.5, pred[:, k] > 0.5
        hard = ((gt[:, k] >= 0.1) & (gt[:, k] <= 0.4)) | ((gt[:, k] >= 0.6) & (gt[:, k] <= 0.9))
        res[f"L{L}"] = {"acc": float(np.mean(pl == lab)), "hard_acc": float(np.mean((pl == lab)[hard])),
                        "hard_n": int(hard.sum())}
    lab = pred > 0.5
    res["nesting_violations"] = float(np.mean(np.any(lab[:, 1:] & ~lab[:, :-1], axis=1)))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--budgets", default="8000,16000,32000")
    ap.add_argument("--n-test", type=int, default=4000)
    ap.add_argument("--gt-seeds", type=int, default=64)
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--compass", type=int, default=0, help="budget of the extra COMPASS-points comparison (0 = off)")
    ap.add_argument("--out", default="s3e_joint.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    budgets = [int(b) for b in a.budgets.split(",")]
    t0 = time.time()
    test_x = np.random.default_rng(12345).random((a.n_test, model.DIM)).astype(np.float32)
    gt = simulate_all(np.repeat(test_x, a.gt_seeds, 0), dev, 12345).reshape(a.n_test, a.gt_seeds, 4).mean(1)
    print(f"ground truth {a.n_test} x {a.gt_seeds} in {time.time()-t0:.0f}s; life fractions "
          + " ".join(f"L{L} {np.mean(gt[:, k] > 0.5):.3f}" for k, L in enumerate(LEVELS)), flush=True)
    res = {"budgets": budgets, "reps": a.reps, "runs": []}
    tq = torch.tensor(test_x, device=dev)
    designs = [("sobol", r) for r in range(a.reps)] + ([("compass_L4", 0)] if a.compass else [])
    for design, r in designs:
        if design == "sobol":
            x = qmc.Sobol(d=model.DIM, scramble=True, seed=1000 + r).random(max(budgets)).astype(np.float32)
            bs = budgets
        else:
            model.LEVEL = 4
            t = time.time()
            x, _, _, _ = bench1b.m_engine3(bench1b.Sim(dev, 2040, 1 << 15), a.compass, 0, dev, soft=True)
            print(f"COMPASS L4 points: {len(x)} in {time.time()-t:.0f}s", flush=True)
            bs = [a.compass]
        y = simulate_all(x, dev, 2000 + r if design == "sobol" else 2999)
        for B in bs:
            xt, yt = torch.tensor(x[:B], device=dev), torch.tensor(y[:B], device=dev)
            t = time.time()
            nets = [surrogate.fit_es_multi(xt, yt[:, k], torch.ones(B, device=dev), [B], seed=7 + r, hidden=256,
                                           max_steps=10000) for k in range(4)]
            sep = np.stack([surrogate.predict_multi(n, tq).cpu().numpy()[0] for n in nets], 1)
            t_sep = time.time() - t
            t = time.time()
            jn = fit_joint(xt, yt, seed=7 + r)
            with torch.no_grad():
                joint = torch.exp(jn.log_p(tq)).cpu().numpy()
            t_joint = time.time() - t
            row = {"design": design, "rep": r, "budget": B, "separate": score(sep, gt), "joint": score(joint, gt),
                   "seconds": {"separate": t_sep, "joint": t_joint}}
            res["runs"].append(row)
            for name in ("separate", "joint"):
                s = row[name]
                print(f"{design} rep {r} B={B // 1000}k {name:>8}: " + " | ".join(
                    f"L{L} {s[f'L{L}']['acc']:.3f}/{s[f'L{L}']['hard_acc']:.3f}" for L in LEVELS)
                      + f" | nesting violations {s['nesting_violations']:.4f} ({row['seconds'][name]:.0f}s)", flush=True)
            json.dump(res, open(a.out, "w"))
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
