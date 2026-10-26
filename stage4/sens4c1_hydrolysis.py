"""Stage 4c-1 hydrolysis-law sensitivity check (chem4c1.py itself is not changed).

The Li & Breaker (1999) law is base catalysis only and its Mg term is extrapolated far below
their 5-50 mM range. The variant floors both terms inside k_nuc:
  pH term  10^(0.983 (max(pH, 5) - 7.4))
  Mg term  (max(Mg, 0.005) / 0.005)^0.8
Baseline and variant are run on the same 2000 ground-truth points of the 4c-1 map
(rng 12345), the same number of seeds, and the same random seeds per batch (paired), at
screening fidelity. Crash-safe like run4c1_map.py (JSON lines per batch, --resume),
with the thermal guard.

  python sens4c1_hydrolysis.py --out-dir sens4c1 --n-test 2000 --seeds 8
"""
import argparse, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c1 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402
import run4c1_map as rm  # noqa: E402

K_NUC_BASE = chem.k_nuc


def k_nuc_floored(T, p):
    return (chem.K_NUC_REF * 10 ** (0.983 * (torch.clamp(p["pH"], min=5.0) - 7.4)) * 10 ** (0.07 * (T - 37))
            * (torch.clamp(p["Mg"], min=0.005) / 0.005) ** 0.8
            * (torch.clamp(p["salinity"], min=0.01) / 0.25) ** p["s_exp"])


def run_arm(a, name, fn, X, dev, guard):
    chem.k_nuc = fn
    try:
        path = os.path.join(a.out_dir, f"{name}.jsonl")
        lines = rm.load_lines(path)
        rm.run_points(path, "test", 0, 0, X, a.seeds, a.batch, dev, "screen", guard,
                      {(d["round"], d["batch"]) for d in lines})
        recs = rm.gather(rm.load_lines(path))
    finally:
        chem.k_nuc = K_NUC_BASE
    e1, e2 = rm.outcomes(recs)
    return e1.reshape(len(X), a.seeds).mean(1), e2.reshape(len(X), a.seeds).mean(1), recs


def bins(X, p2):
    d = chem.decode(torch.tensor(X))
    pH, Mg = d["pH"].numpy(), d["Mg"].numpy()
    out = {}
    for lab, m in (("pH 2-5", pH < 5), ("pH 5-8", (pH >= 5) & (pH < 8)), ("pH 8-12", pH >= 8),
                   ("Mg 1e-4..1e-3", Mg < 1e-3), ("Mg 1e-3..1e-2", (Mg >= 1e-3) & (Mg < 1e-2)),
                   ("Mg 1e-2..1e-1", Mg >= 1e-2)):
        out[lab] = {"n": int(m.sum()), "E2_rate": float((p2[m] > 0.5).mean())}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "sens4c1"))
    ap.add_argument("--n-test", type=int, default=2000)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--batch", type=int, default=512)
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    cfg = {"n_test": a.n_test, "seeds": a.seeds, "batch": a.batch}
    cfg_path = os.path.join(a.out_dir, "config.json")
    old = ckpt.load_json(cfg_path)
    if old is not None and not a.resume:
        sys.exit(f"{a.out_dir} already has a run; use --resume")
    if old is not None and old != cfg:
        sys.exit(f"settings differ from the saved run: {old} vs {cfg}")
    ckpt.atomic_json(cfg, cfg_path)
    dev = torch.device(a.device)
    guard = thermal.Guard(job=f"4c1 {os.path.basename(a.out_dir)}")
    t0 = time.time()
    X = np.random.default_rng(12345).random((a.n_test, chem.DIM)).astype(np.float32)
    b1, b2, _ = run_arm(a, "baseline", K_NUC_BASE, X, dev, guard)
    v1, v2, _ = run_arm(a, "variant", k_nuc_floored, X, dev, guard)
    res = {"config": cfg, "variant": "k_nuc with pH floored at 5 and Mg floored at 5 mM (both terms inside the law)",
           "baseline": {"E1": float((b1 > 0.5).mean()), "E2": float((b2 > 0.5).mean()), "bins": bins(X, b2)},
           "variant_floored": {"E1": float((v1 > 0.5).mean()), "E2": float((v2 > 0.5).mean()), "bins": bins(X, v2)},
           "E2_majority_flips": {"share": float(((b2 > 0.5) != (v2 > 0.5)).mean()),
                                 "lost": int(((b2 > 0.5) & ~(v2 > 0.5)).sum()),
                                 "gained": int((~(b2 > 0.5) & (v2 > 0.5)).sum())},
           "E1_majority_flips_share": float(((b1 > 0.5) != (v1 > 0.5)).mean()),
           "thermal": guard.summary(), "seconds_this_session": time.time() - t0}
    ckpt.atomic_json(res, os.path.join(a.out_dir, "result.json"))
    print(res)


if __name__ == "__main__":
    main()
