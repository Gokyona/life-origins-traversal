"""Stage 4c-3: scan for a reference condition where L2 occurs (base of validation tests 2-5).

A grid around the 4c-2 reference (T, pH, thermo, T_m8) gave no L2: wherever copying is fast
enough (q(1) = r / (r + 1e6 /h) needs r of order 1e6 /h, i.e. high c_eff M*), initiation fills
the sites with dimers. So the reference is taken from random conditions instead:
  stage 1: 4096 Sobol points (seed 2026) of the full 42-dimensional space at screening
           fidelity with both paired controls, as in a map batch;
  stage 2: the points with L2_strict (up to 8, by information share above the control) are
           rerun at validation fidelity (16 sites x 256 slots, 600 target steps), 16 reps each;
           the point with the highest L2_strict rate (ties: larger share above control) is
           REF_L2, stored as its u vector and physical values.
Crash-safe: stage-1 batches and stage-2 points are saved as they finish (--resume by
skipping them).   python scan4c3.py [out.json] [--box]
--box: the same scan inside the sub-box where copies can be long and strands can separate
(after the full-space scan found L2_strict only at noise level): T between T_m8 - 15 and T_m8,
cycling depth >= 0.5, k_t >= 1e4, thermo >= 1e3, x_nuc >= 0.5; all other dimensions full range.
--box2: after the copy start rule changed (copies start at the template's first letter), the box
of --box plus: hot-phase temperature T + 10 depth above T_m8, and lifetime(L_H) / period >= 1
(cycles faster than an L_H-mer lives: turnover with adsorption protection, cycle-averaged
hydrolysis); points are drawn from 65536 Sobol points of the --box transform by rejection,
the first 4096 that pass.
"""
import math, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c3 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402

dev = torch.device("xpu")
guard = thermal.Guard(job="4c3 scan")
N1, SEED = 4096, 2026
BATCH = int(os.environ.get("LOT_BATCH", "512"))     # chosen by the memory pre-flight (preflight_mem.py)
FID2, REPS2, TOP = (16, 256, 600), 16, 8
LIFE_RATIO = None         # --box2: lifetime(L_H) / period of each point


def stage1_batch(X, b):
    u = torch.tensor(X[b * BATCH:(b + 1) * BATCH], dtype=torch.float32, device=dev)
    g = torch.Generator(device=dev)
    g.manual_seed(70000 + b)
    o = chem.simulate(u, g, "screen")
    s1, _, _ = chem.paired_control(u, 70000 + b, "screen", o)
    s2, sm, sc = chem.paired_control_l2(u, 70000 + b, "screen", o)
    torch.xpu.empty_cache()
    guard.check()
    pts = []
    for i in range(len(u)):
        if bool(o["L2_loose"][i]):
            pts.append({"i": b * BATCH + i, "L2_strict": bool(s2[i]), "share": float(sm[i]),
                        "share_ctl": None if math.isnan(float(sc[i])) else float(sc[i]),
                        "span": float(o["L2_span_turnovers"][i]), "copies": float(o["copies"][i]),
                        "copy_len": float(o["copy_mean_len"][i]), "E2": bool(o["E2"][i]), "L1": bool(s1[i])})
    cyc = o["cycle_averaged"] | (u[:, chem.NAMES.index("depth")].cpu() <= 0)
    long_ = o["copies_long"] > 0
    paired = o["paired_template_share"] > 0.5
    e1 = o["E1"]
    return {"E1": float(e1.float().mean()), "L2_loose": int(o["L2_loose"].sum()), "L2_strict": int(s2.sum()),
            "copying_runs": int((o["copies"] > 0).sum()),
            # which condition fails (counts over the batch): copies reaching L_H, templates paired at
            # the end, resolved cycling (depth > 0 and phases not averaged)
            "n_copies_reach_LH": int(long_.sum()), "n_templates_paired_end": int(paired.sum()),
            "n_resolved_cycling": int((~cyc).sum()), "n_long_and_unpaired": int((long_ & ~paired).sum()),
            "n_long_unpaired_resolved_E1": int((long_ & ~paired & ~cyc & e1).sum()),
            "copies_long_pct10_50_90": [float(x) for x in np.percentile(o["copies_long"].numpy(), [10, 50, 90])],
            "points": pts}


def stage2_point(x, seed):
    u = torch.tensor(np.repeat(x[None, :], REPS2, 0), dtype=torch.float32, device=dev)
    g = torch.Generator(device=dev)
    g.manual_seed(seed)
    o = chem.simulate(u, g, FID2)
    s2, sm, sc = chem.paired_control_l2(u, seed, FID2, o)
    guard.check()
    margin = (sm - torch.nan_to_num(sc, nan=0.0))[o["L2_loose"]]
    return {"L2_strict": int(s2.sum()), "L2_loose": int(o["L2_loose"].sum()),
            "share_mean": float(o["L2_share"].mean()), "margin_mean": float(margin.mean()) if margin.numel() else 0.0,
            "copies_mean": float(o["copies"].mean()), "copy_len_mean": float(o["copy_mean_len"].mean()),
            "E1": int(o["E1"].sum()), "E2": int(o["E2"].sum())}


def main():
    global LIFE_RATIO
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = args[0] if args else os.path.join(HERE, "scan4c3.json")
    res = ckpt.load_json(path) or {"stage1": {}, "stage2": {}}
    X = qmc.Sobol(d=chem.DIM, scramble=True, seed=SEED).random(N1).astype(np.float32)
    if "--box" in sys.argv or "--box2" in sys.argv:
        ix = {n: chem.NAMES.index(n) for n in ("T", "T_m8", "depth", "k_t", "thermo", "x_nuc")}
        Tm = 10 + 80 * X[:, ix["T_m8"]]
        X[:, ix["T"]] = np.clip((Tm - 15 * X[:, ix["T"]]) / 150.0, 0.0, 1.0)
        X[:, ix["depth"]] = 0.5 + 0.5 * X[:, ix["depth"]]
        X[:, ix["k_t"]] = 0.75 + 0.25 * X[:, ix["k_t"]]
        X[:, ix["thermo"]] = 0.375 + 0.625 * X[:, ix["thermo"]]
        X[:, ix["x_nuc"]] = 0.5 + 0.5 * X[:, ix["x_nuc"]]
    if "--box2" in sys.argv:
        Xb = qmc.Sobol(d=chem.DIM, scramble=True, seed=SEED).random(1 << 16).astype(np.float32)
        ix = {n: chem.NAMES.index(n) for n in ("T", "T_m8", "depth", "k_t", "thermo", "x_nuc")}
        Tm = 10 + 80 * Xb[:, ix["T_m8"]]
        T = Tm - 15 * Xb[:, ix["T"]]
        Xb[:, ix["T"]] = np.clip(T / 150.0, 0.0, 1.0)
        dmin = np.clip((Tm - T) / 10.0 + 0.02, 0.0, 1.0)                    # hot phase T + 10 d > T_m8
        Xb[:, ix["depth"]] = np.maximum(0.5 + 0.5 * Xb[:, ix["depth"]], dmin)
        Xb[:, ix["k_t"]] = 0.75 + 0.25 * Xb[:, ix["k_t"]]
        Xb[:, ix["thermo"]] = 0.375 + 0.625 * Xb[:, ix["thermo"]]
        Xb[:, ix["x_nuc"]] = 0.5 + 0.5 * Xb[:, ix["x_nuc"]]
        p = chem.decode(torch.tensor(Xb))
        tab, _, theta, _ = chem.phase_table(p)
        life = 1.0 / chem.turnover_rate(chem.L_H, tab["knuc"][:, 2], p["k_UV"] * p["UV"], 1 / p["tau"], theta, p["p_min"])
        keep = ((life / p["period"]) >= 1.0).numpy()
        sel = np.nonzero(keep)[0][:N1]
        print(f"box2: {int(keep.sum())} of {len(Xb)} pass lifetime/period >= 1; using {len(sel)}", flush=True)
        X = Xb[sel]
        LIFE_RATIO = (life / p["period"]).numpy()[sel]
    if LIFE_RATIO is not None:
        res["life_over_period_pct10_50_90"] = [float(x) for x in np.percentile(LIFE_RATIO, [10, 50, 90])]
    for b in range((len(X) + BATCH - 1) // BATCH):
        if str(b) in res["stage1"]:
            continue
        t = time.time()
        r = stage1_batch(X, b)
        r["seconds"] = time.time() - t
        res["stage1"][str(b)] = r
        ckpt.atomic_json(res, path)
        print("batch", b, {k: r[k] for k in ("E1", "copying_runs", "n_copies_reach_LH", "n_templates_paired_end",
                                             "n_resolved_cycling", "n_long_unpaired_resolved_E1", "L2_loose", "L2_strict")},
              f"({r['seconds']:.0f}s)",
              flush=True)
    cands = [p for r in res["stage1"].values() for p in r["points"] if p["L2_strict"]]
    cands.sort(key=lambda p: -(p["share"] - (p["share_ctl"] or 0.0)))
    for p in cands[:TOP]:
        if str(p["i"]) in res["stage2"]:
            continue
        t = time.time()
        r = stage2_point(X[p["i"]], 80000 + p["i"])
        r["seconds"] = time.time() - t
        res["stage2"][str(p["i"])] = r
        ckpt.atomic_json(res, path)
        print("point", p["i"], r, flush=True)
    if res["stage2"]:
        best = max(res["stage2"], key=lambda k: (res["stage2"][k]["L2_strict"], res["stage2"][k]["margin_mean"]))
        x = X[int(best)]
        phys = chem.decode(torch.tensor(x[None, :]))
        res["ref"] = {"i": int(best), "u": [float(v) for v in x], "physical": {k: float(v[0]) for k, v in phys.items()},
                      "stage2": res["stage2"][best]}
        print("REF_L2 = point", best, res["stage2"][best], flush=True)
    res["thermal"] = guard.summary()
    ckpt.atomic_json(res, path)


if __name__ == "__main__":
    main()
