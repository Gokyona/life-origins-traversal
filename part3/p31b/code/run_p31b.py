"""P3-1b batch runner (chem_p31b.py), PAIR=None (Stage 4 chemistry): simulates a point set in batches with the L1 and
L2 paired controls after every batch and writes one JSON line per batch (fsync). Crash-safe as run_p31.py: `--resume`
skips finished batches (a torn last line is dropped), each batch's seed comes from (tag, batch), so a resumed run gives
the same simulations as an uninterrupted one. Thermal guard between batches. Every per-run output of simulate
(1-D [B] tensors, and [B, k] with k <= 64) is stored, so detector perturbations can be recomputed later.
    python run_p31b.py --set coex --tile 5x32 --fid 16,256,300 --out-dir kill_test_p31b/ka [--multi 1] [--resume]
"""
import argparse, json, math, os, sys, time, zlib
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_p31b as chem  # noqa: E402
import pointsets_p31b as PS  # noqa: E402
import thermal  # noqa: E402
from run_p31 import load_lines, append_line  # noqa: E402


def seed_of(tag, batch):
    return 3150 * 100000 + (zlib.crc32(tag.encode()) % 1000) * 100 + batch


def fl(x):
    x = float(x)
    return None if math.isnan(x) else (("inf" if x > 0 else "-inf") if math.isinf(x) else x)


def per_run(o, B):
    d = {}
    for k, v in o.items():
        if isinstance(v, torch.Tensor) and v.dim() >= 1 and v.shape[0] == B and v.dim() <= 2 and \
                (v.dim() == 1 or v.shape[1] <= 64):
            d[k] = v.detach().cpu().double()
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="coex")
    ap.add_argument("--tile", default="5x32", help="batches x batch size; points are tiled to fill them")
    ap.add_argument("--fid", default="16,256,300")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--multi", type=int, default=1)
    ap.add_argument("--tag", default="kill")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    chem.PAIR = None
    chem.MULTI_START = bool(a.multi)
    fid = tuple(int(x) for x in a.fid.split(",")) if "," in a.fid else a.fid
    nb, bs = (int(x) for x in a.tile.split("x"))
    names, X = PS.coex_points()
    X = X[:, :chem.DIM]
    idx_all = np.arange(nb * bs) % len(names)
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, "runs.jsonl")
    if os.path.exists(path) and not a.resume:
        sys.exit(f"{path} exists: use --resume")
    guard = thermal.Guard(job=f"run_p31b {a.out_dir}")
    dev = torch.device("xpu")
    done = {d["batch"] for d in load_lines(path)}
    for b in range(nb):
        if b in done:
            continue
        idx = idx_all[b * bs:(b + 1) * bs]
        sd = seed_of(a.tag, b)
        g = torch.Generator(device=dev)
        g.manual_seed(sd)
        t = time.time()
        u = torch.tensor(X[idx], dtype=torch.float32, device=dev)
        o = chem.simulate(u, g, fid, record=True)
        l1 = chem.paired_control(u, sd, fid, o)
        l2 = chem.paired_control_l2(u, sd, fid, o, both=True)
        pr = per_run(o, len(idx))
        for k, v in zip(("L1_strict", "ab_cat", "ab_ctl"), l1):
            pr["pc_" + k] = v.detach().cpu().double()
        for k, v in zip(("L2_strict", "S_main", "S_ctl"), l2):
            pr["pc_" + k] = v.detach().cpu().double()
        recs = {k: [fl(x) for x in v.flatten().tolist()] if v.dim() == 1 else [[fl(x) for x in r] for r in v.tolist()]
                for k, v in pr.items()}
        append_line(path, {"batch": b, "points": [names[i] for i in idx], "n": len(idx), "seed": sd,
                           "multi": a.multi, "seconds": time.time() - t, "runs": recs})
        print(f"batch {b} n={len(idx)} {time.time() - t:.1f}s", flush=True)
        del o, u, l1, l2
        torch.xpu.empty_cache()
        guard.check(force=True)
    print("done", path, guard.summary(), flush=True)


if __name__ == "__main__":
    main()
