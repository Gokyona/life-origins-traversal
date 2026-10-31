"""P3-1b section 3 on the 7 confirmed L2 points (Stage 4 chemistry, PAIR=None, chem_p31b): one batch per variant
with all points x 16 seeds packed together, the published confirm fidelity (16, 256, 600) and the L2 paired
random-copy controls (paired_control_l2, both=True). The same seed for every variant (common random numbers).
One fsynced JSON line per (variant, chunk); --resume skips finished ones; thermal guard between batches.
    python run_s3_l2.py --out-dir s3_l2 --variants base,A0,A1,A2,A2a,A4,A6 [--points ...] [--resume]
Variants (module switches of chem_p31b, all off = base):
  base  published rule;  A0  POISSON_INV;  A1  G_CAP 16;  A2  MULTI_START (+ slot fix + POISSON_INV);
  A2a  SLOT_ROUND_FIX;  A4  DUPLEX_HYDRO_FACTOR 0.1;  A5  DT_PERIOD_FRAC 0.25;  A6  LONG_WASH_FACTOR 0.5
  det_LH5 / det_LH7 / det_SPAN4 / det_SPAN6: base dynamics replayed with a detector constant changed
"""
import argparse, json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_p31b as chem  # noqa: E402
import thermal  # noqa: E402
import validate_p31 as V1  # noqa: E402
from run_p31 import load_lines, append_line  # noqa: E402
from run_p31b import per_run, fl  # noqa: E402

S4 = os.path.join(HERE, "..", "..", "lot-s4")
FID, REPS, SEED = (16, 256, 600), 16, 83100
DEFAULTS = {"MULTI_START": False, "SLOT_ROUND_FIX": False, "POISSON_INV": False, "G_CAP": 4,
            "DUPLEX_HYDRO_FACTOR": 1.0, "LONG_WASH_FACTOR": 1.0, "DT_PERIOD_FRAC": None,
            "L_H": 6, "L2_SPAN": 5.0, "SLOT_CUT_RANDOM": False, "MS_MELT_FLOOR": True, "MS_ORDER_SLOT": False}
VARIANTS = {"base": {}, "A0": {"POISSON_INV": True}, "A1": {"G_CAP": 16}, "A2": {"MULTI_START": True},
            "A2a": {"SLOT_ROUND_FIX": True}, "A4": {"DUPLEX_HYDRO_FACTOR": 0.1}, "A5": {"DT_PERIOD_FRAC": 0.25},
            "A6": {"LONG_WASH_FACTOR": 0.5}, "A0_r6": {"POISSON_INV": True}, "A2_r6": {"MULTI_START": True},
            "base_rec": {}, "A2_rec": {"MULTI_START": True}, "A2b": {"SLOT_CUT_RANDOM": True}, "A2c": {"MULTI_START": True, "MS_MELT_FLOOR": False},
            "A2d": {"MULTI_START": True, "MS_ORDER_SLOT": True},
            "A2e": {"MULTI_START": True, "MS_MELT_FLOOR": False, "MS_ORDER_SLOT": True},
            "base_conc": {}, "A2b_conc": {"SLOT_CUT_RANDOM": True}, "A2c_conc": {"MULTI_START": True, "MS_MELT_FLOOR": False}, "A5_part1": {"DT_PERIOD_FRAC": 0.25},
            "A5_part2": {"DT_PERIOD_FRAC": 0.25}, "det_LH5": {"L_H": 5}, "det_LH7": {"L_H": 7},
            "det_SPAN4": {"L2_SPAN": 4.0}, "det_SPAN6": {"L2_SPAN": 6.0}}


def l2_points():
    """name -> (u_model [42], published L2_strict of 16, published seed)."""
    pts = {}
    bx = V1._box2_points()
    for i in (1558, 3709):
        pts[f"p{i}"] = (bx[i], None, None)
    for f in ("confirm4c3.json", "confirm4c3_rep1.json"):
        d = json.load(open(os.path.join(S4, f)))
        for k, v in d["points"].items():
            if k in ("engine435", "engine5704", "engine10653", "engine_rep1_1926", "engine_rep1_9823"):
                pts[k] = (np.array(v["u_model"], dtype=np.float32), v["as_is"]["L2_strict"], v["seed"])
    return pts


def set_variant(name):
    for k, v in DEFAULTS.items():
        setattr(chem, k, v)
    for k, v in VARIANTS[name].items():
        setattr(chem, k, v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--variants", default="base,A0,A1,A2,A2a,A4,A6")
    ap.add_argument("--points", default="")
    ap.add_argument("--chunk", type=int, default=0, help="runs per batch (0: all points of a variant in one batch)")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    chem.PAIR = None
    pts = l2_points()
    names = [n for n in pts if not a.points or n in a.points.split(",")]
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, "runs.jsonl")
    if os.path.exists(path) and not a.resume:
        sys.exit(f"{path} exists: use --resume")
    json.dump({"fid": FID, "reps": REPS, "seed": SEED, "points": {n: {"published_L2_strict": pts[n][1],
               "published_seed": pts[n][2], "u_model": [float(x) for x in pts[n][0]]} for n in names},
               "variants": {v: VARIANTS[v] for v in a.variants.split(",")},
               "chem_p31b_md5": os.popen(f"md5sum {os.path.join(HERE, 'chem_p31b.py')}").read().split()[0]},
              open(os.path.join(a.out_dir, f"config_{a.variants.replace(',', '_')}.json"), "w"), indent=1)
    done = {(d["variant"], d["chunk"]) for d in load_lines(path)}
    guard = thermal.Guard(job=f"p31b s3 L2 {a.out_dir}")
    dev = torch.device("xpu")
    X = np.concatenate([np.repeat(pts[n][0][None, :], REPS, 0) for n in names]).astype(np.float32)
    lab = [n for n in names for _ in range(REPS)]
    cs = a.chunk or len(X)
    for v in a.variants.split(","):
        for ci, i0 in enumerate(range(0, len(X), cs)):
            if (v, ci) in done:
                continue
            set_variant(v)
            sd = SEED + ci
            g = torch.Generator(device=dev)
            g.manual_seed(sd)
            t = time.time()
            u = torch.tensor(X[i0:i0 + cs], dtype=torch.float32, device=dev)
            torch.xpu.reset_peak_memory_stats()
            o = chem.simulate(u, g, FID, record=True)
            l2 = chem.paired_control_l2(u, sd, FID, o, both=True)
            pr = per_run(o, len(u))
            for k, x in zip(("L2_strict", "S_main", "S_ctl"), l2):
                pr["pc_" + k] = x.detach().cpu().double()
            recs = {k: [fl(x) for x in t_.flatten().tolist()] if t_.dim() == 1 else [[fl(x) for x in r] for r in t_.tolist()]
                    for k, t_ in pr.items()}
            sec = time.time() - t
            append_line(path, {"variant": v, "chunk": ci, "points": lab[i0:i0 + cs], "n": len(u), "seed": sd,
                               "seconds": sec, "peak_xpu_reserved_GiB": torch.xpu.max_memory_reserved() / 2 ** 30,
                               "runs": recs})
            st = np.array(recs["pc_L2_strict"]) > 0.5
            summ = {n: int(st[[j for j in range(len(u)) if lab[i0 + j] == n]].sum()) for n in dict.fromkeys(lab[i0:i0 + cs])}
            print(f"{v} chunk {ci} n={len(u)} {sec:.0f}s L2_strict {summ}", flush=True)
            del o, u, l2
            torch.xpu.empty_cache()
            guard.check(force=True)
    set_variant("base")
    print("done", path, guard.summary(), flush=True)


if __name__ == "__main__":
    main()
