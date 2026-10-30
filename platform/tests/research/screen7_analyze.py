"""Analysis of screen7 (single-genotype vesicles: A baseline, C, D, Dp x seeds 1-2) against
test_results_2026-10-29/runs/screen7_criteria.json (md5 2c13352a, written before the runs).
    python tests/screen7_analyze.py <runs dir>"""
import json, math, os, sys
import numpy as np

B = sys.argv[1]
YEAR = 8766.0
ARMS = ("A", "C", "D", "Dp")


def mean(h, k):
    v = [x[k] for x in h if x.get(k) is not None]
    return float(np.mean(v)) if v else None


def cell(name):
    d = os.path.join(B, name)
    st = json.load(open(os.path.join(d, "status.json")))
    h = [json.loads(l) for l in open(os.path.join(d, "history.jsonl"))]
    h = [x for x in h if x.get("z_all_real") is not None and x.get("z_all_control") is not None]
    dw = st["darwin"]; tr = dw["trait_response"]
    half = len(h) // 2; late = h[half:]
    yrs = np.array([x["sim_h"] / YEAR for x in h]); diff = np.array([x["z_all_real"] - x["z_all_control"] for x in h])
    dist = np.mean([x["ves_dist"] for x in late if x.get("ves_dist")], axis=0).tolist() if any(x.get("ves_dist") for x in late) else None
    sc_site = [x["ves_scored_mean"] * x["vesicles"] for x in late if x.get("ves_scored_mean") is not None]
    return {"sim_years": round(st["time"]["sim_years"], 1), "reports": len(h),
            "diff_all_late": float(diff[half:].mean()), "diff_all_slope_per_100y": float(np.polyfit(yrs[half:], diff[half:], 1)[0] * 100),
            "z_all_control_late_sd": float(np.std([x["z_all_control"] for x in late], ddof=1)),
            "within_share_late": mean(late, "ves_within_share"), "zbar_sd_late": mean(late, "ves_zbar_sd"),
            "occupied_late": mean(late, "ves_occupied_share"), "scored_per_vesicle_late": mean(late, "ves_scored_mean"),
            "scored_per_site_late": float(np.mean(sc_site)) if sc_site else None,
            "dist_0_1_2_3to4_5plus_late": dist, "share_vesicles_1_2": (dist[1] + dist[2]) if dist else None,
            "chain_share_in_1_2_late": mean(late, "ves_chain_share_1_2"),
            "div_per_vesicle_year_late": mean(late, "div_per_vesicle_year"), "div_vs_zbar_corr_late": mean(late, "div_vs_zbar_corr"),
            "mother_daughter_corr_late": mean(late, "mother_daughter_corr"), "seg_counts_final": h[-1].get("seg_counts"),
            "darwin_status": dw["status"], "shifts": tr["shifts"], "plateaus": tr["plateaus"],
            "clade_counted": dw.get("clade_counted", 0), "clade_events_real": sum(1 for e in dw.get("clade_sweeps", []) if e["half"] == "real"),
            "lineage_sweeps": len(dw["sweeps"]), "top1_lineage_real_max": max((x.get("top1_share_real") or 0) for x in h),
            "clade_top1_real_max": max((x.get("clade_top1_share_real") or 0) for x in h),
            "step_ms_mean": float(np.mean([x["step_ms"] for x in h if x.get("step_ms")])), "gpu_temp_max": max((x.get("gpu_temp_c") or 0) for x in h)}


res = {"criteria_file": "screen7_criteria.json (2c13352a)", "cells": {}}
for s in (1, 2):
    for a in ARMS:
        res["cells"][f"{a}_s{s}"] = cell(f"screen7_{a}_s{s}")
C = res["cells"]
ev = {}
for a in ("C", "D", "Dp"):
    e = {}
    for s in (1, 2):
        X, A = C[f"{a}_s{s}"], C[f"A_s{s}"]
        pooled = math.sqrt((X["z_all_control_late_sd"] ** 2 + A["z_all_control_late_sd"] ** 2) / 2)
        e[f"s{s}"] = {"1": X["clade_counted"] >= 2, "2": X["within_share_late"] is not None and X["within_share_late"] < 0.30,
                      "3": X["diff_all_late"] - A["diff_all_late"] > 2 * pooled, "3_margin": X["diff_all_late"] - A["diff_all_late"], "3_threshold": 2 * pooled,
                      "4": X["diff_all_slope_per_100y"] > 0,
                      "5": (X["occupied_late"] or 0) >= 0.5 and (X["scored_per_site_late"] or 0) >= 0.5 * (A["scored_per_site_late"] or 1)}
    e["both"] = {k: e["s1"][k] and e["s2"][k] for k in ("1", "2", "3", "4", "5")}
    ev[a] = e
res["evaluation"] = ev
json.dump(res, open(os.path.join(B, "screen7_results.json"), "w"), indent=1)
for k, c in C.items():
    print(k, {x: (round(v, 4) if isinstance(v, float) else v) for x, v in c.items()})
print(json.dumps(ev, indent=1))
