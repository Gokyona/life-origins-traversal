"""Stage 4c-2 detector-sensitivity check at REF_L1 (validate4c2.REF_L1), same seeds per arm:
  (a) current: 4 RAF sites, refresh every 8 steps, catalysis after E1
  (b) 8 RAF sites, otherwise as (a)
  (c) 8 RAF sites, refresh every step
  (d) as (c), catalysis from t = 0
Per arm: L1_loose / L1_strict rates over 16 reps (validation fidelity 16 sites x 256 slots, as
in test 3) and the cost per run at screening fidelity (512 Sobol points, with the paired control).
Adoption: the arm with the highest L1_strict rate among those with rate >= 2 x arm (a), rate >
arm (a) and cost <= 0.3 s per run; otherwise arm (a). Writes result.json with "adopt"
(read by run4c2_map.py --detector). Crash-safe: one JSON file per finished arm, --resume.

  python sens4c2_detector.py --out-dir sens4c2
"""
import argparse, os, sys, time
import torch
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c2 as chem  # noqa: E402
import ckpt  # noqa: E402
import validate4c2 as v  # noqa: E402

ARMS = {"a": {"ns_raf": 4, "refresh": 8, "cat_after_e1": True},
        "b": {"ns_raf": 8, "refresh": 8, "cat_after_e1": True},
        "c": {"ns_raf": 8, "refresh": 1, "cat_after_e1": True},
        "d": {"ns_raf": 8, "refresh": 1, "cat_after_e1": False}}
REPS = 16
COST_LIMIT_MS = 300.0


def set_arm(s):
    chem.NS_RAF, chem.REFRESH, chem.CAT_AFTER_E1 = s["ns_raf"], s["refresh"], s["cat_after_e1"]


def cost(dev):
    u = torch.tensor(qmc.Sobol(chem.DIM, scramble=True, seed=77).random(512), dtype=torch.float32, device=dev)
    g = torch.Generator(device=dev)
    g.manual_seed(61)
    torch.xpu.synchronize()
    t = time.time()
    o = chem.simulate(u, g, "screen", record=True)
    strict, _, _ = chem.paired_control(u, 61, "screen", o)
    torch.xpu.synchronize()
    return 1000 * (time.time() - t) / 512, float(o["L1_loose"].float().mean()), float(strict.float().mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.join(HERE, "sens4c2"))
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    cfg = {"arms": ARMS, "reps": REPS, "reference": "validate4c2.REF_L1", "fidelity_L1": list(v.FID),
           "cost": "512 Sobol points, screening fidelity, with paired control"}
    ckpt.atomic_json(cfg, os.path.join(a.out_dir, "config.json"))
    dev = torch.device("xpu")
    res = {}
    for name, s in ARMS.items():
        path = os.path.join(a.out_dir, f"arm_{name}.json")
        r = ckpt.load_json(path)
        if r is None:
            set_arm(s)
            t = time.time()
            x = v.run(v.REF_L1, 301, reps=REPS)
            t_l1 = time.time() - t
            ms, loose_sobol, strict_sobol = cost(dev)
            r = {"settings": s, "L1_loose_rate": sum(x["L1_loose"]) / REPS, "L1_strict_rate": x["L1_rate"],
                 "RAF_span_turnovers": x["RAF_span_turnovers"], "RAF_size": x["RAF_size"],
                 "member_abundance_cat": x["member_abundance_cat"], "member_abundance_control": x["member_abundance_control"],
                 "seconds_REF_L1_16reps": t_l1, "ms_per_run_screening": ms,
                 "sobol512_L1_loose": loose_sobol, "sobol512_L1_strict": strict_sobol}
            ckpt.atomic_json(r, path)
            v.guard.check()
            print(name, {k: r[k] for k in ("L1_loose_rate", "L1_strict_rate", "ms_per_run_screening")}, flush=True)
        res[name] = r
    base = res["a"]["L1_strict_rate"]
    ok = [n for n, r in res.items() if n != "a" and r["L1_strict_rate"] >= 2 * base and r["L1_strict_rate"] > base
          and r["ms_per_run_screening"] <= COST_LIMIT_MS]
    pick = max(ok, key=lambda n: (res[n]["L1_strict_rate"], -res[n]["ms_per_run_screening"])) if ok else "a"
    out = {"arms": res, "adopted_arm": pick, "adopt": ARMS[pick],
           "rule": "highest L1_strict among arms with rate >= 2 x arm a, > arm a, cost <= 300 ms/run; else arm a",
           "thermal": v.guard.summary()}
    ckpt.atomic_json(out, os.path.join(a.out_dir, "result.json"))
    print("adopted", pick, ARMS[pick])


if __name__ == "__main__":
    main()
