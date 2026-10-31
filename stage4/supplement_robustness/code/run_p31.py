"""P3-1 batch runner (chem_p31.py): simulates a Sobol point set of one pair chemistry in batches, with the
L1 and L2 paired controls after every batch, and writes one JSON line per batch (fsync). Crash-safe as the
Stage 4 map runners: `--resume` skips finished batches (a torn last line is dropped), and each batch's seed
comes from (pair, batch), so a resumed run gives the same simulations as an uninterrupted one. Thermal guard
(thermal.py, Stage 4 GPU lock) between batches. Used for validation test 9 (kill and resume) and the cost
measurement; the P3-1 maps (COMPASS) are a later step.
    python run_p31.py --pair AS2U --out-dir runs_as2u --n 512 --batch 256 [--fid screen] [--resume]
"""
import argparse, json, math, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_p31 as chem  # noqa: E402
import thermal  # noqa: E402

PAIR_IX = {"AS2U": 1, "GC": 2, "AU": 3}


def seed_of(pair, batch):
    return (31 * 100 + PAIR_IX[pair]) * 100000 + batch


def load_lines(path):
    if not os.path.exists(path):
        return []
    good, bad = [], False
    with open(path) as f:
        for ln in f:
            try:
                good.append(json.loads(ln))
            except json.JSONDecodeError:
                bad = True
    if bad:
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            for d in good:
                f.write(json.dumps(d) + "\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    return good


def append_line(path, d):
    with open(path, "a") as f:
        f.write(json.dumps(d) + "\n")
        f.flush()
        os.fsync(f.fileno())


def fl(x):
    x = float(x)
    return None if math.isnan(x) else (("inf" if x > 0 else "-inf") if math.isinf(x) else x)


def run_summary(o, i, l1, l2):
    h = o["length_hist"][i]
    n = float(h[2:].sum())
    return {
        "E1": bool(o["E1"][i]), "E2": bool(o["E2"][i]), "liquid": bool(o["liquid"][i]),
        "L1_loose": bool(o["L1_loose"][i]), "L1_strict": bool(l1[0][i]),
        "L2_loose": bool(o["L2_loose"][i]), "L2_strict": bool(l2[0][i]), "L2_share": fl(l2[1][i]), "L2_share_ctl": fl(l2[2][i]),
        "L2_share_last5": fl(o["L2_share_last5"][i]), "L2_span": fl(o["L2_span_turnovers"][i]),
        "reached_end": bool(o["reached_end"][i]), "t_end": fl(o["t_end"][i]), "steps": int(o["steps"][i]),
        "conc": fl(o["conc_mstar"][i]), "n_chains": n,
        "mean_len": float((h[2:] * torch.arange(2, chem.L_MAX + 1)).sum() / n) if n else 0.0,
        "copies": fl(o["copies"][i]), "copies_long": fl(o["copies_long"][i]), "copy_len": fl(o["copy_mean_len"][i]),
        "copy_added": fl(o["copy_mean_added"][i]), "copy_err": fl(o["copy_error_rate"][i]),
        "copy_r": fl(o["copy_r_mean"][i]), "paired_share": fl(o["paired_template_share"][i]),
        "R": fl(o["long_copies_per_template_per_lifetime"][i]), "growth_h": fl(o["lineage_growth_rate_h"][i]),
        "t_L2_est_h": fl(o["t_L2_est_h"][i]), "growing": bool(o["growing"][i]),
        "pulse_share": fl(o["pulse_time_share"][i]), "pulses_resolved": bool(o["pulses_resolved"][i]),
        "hist": [int(v) for v in h[:33]] + [int(h[33:].sum())],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", default="AS2U", choices=list(PAIR_IX))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n", type=int, default=512)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--fid", default="screen")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    chem.PAIR = a.pair
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, f"runs_{a.pair}.jsonl")
    if os.path.exists(path) and not a.resume:
        sys.exit(f"{path} exists: use --resume")
    guard = thermal.Guard(job=f"run_p31 {a.pair} {a.out_dir}")
    dev = torch.device("xpu")
    X = qmc.Sobol(d=chem.DIM_P31, scramble=True, seed=3100 + PAIR_IX[a.pair]).random(a.n).astype(np.float32)
    done = {d["batch"] for d in load_lines(path)}
    nb = (a.n + a.batch - 1) // a.batch
    for b in range(nb):
        if b in done:
            continue
        idx = np.arange(b * a.batch, min(a.n, (b + 1) * a.batch))
        sd = seed_of(a.pair, b)
        g = torch.Generator(device=dev)
        g.manual_seed(sd)
        t = time.time()
        u = torch.tensor(X[idx], dtype=torch.float32, device=dev)
        o = chem.simulate(u, g, a.fid, record=True)
        l1 = chem.paired_control(u, sd, a.fid, o)
        l2 = chem.paired_control_l2(u, sd, a.fid, o)
        recs = [run_summary(o, i, l1, l2) for i in range(len(idx))]
        append_line(path, {"batch": b, "i0": int(idx[0]), "n": len(idx), "seed": sd, "seconds": time.time() - t,
                           "runs": recs})
        print(f"batch {b} n={len(idx)} {time.time() - t:.1f}s", flush=True)
        del o, u, l1, l2
        torch.xpu.empty_cache()
        guard.check(force=True)
    print("done", path, guard.summary(), flush=True)


if __name__ == "__main__":
    main()
