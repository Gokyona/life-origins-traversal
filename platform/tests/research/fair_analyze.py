"""Analysis of the 2x2 fair-rationing x protected-domain screen (29 Oct 2026), against the pre-registered criteria in
test_results_2026-10-29/runs/fair_criteria.json. Reads each cell's status.json and history.jsonl; writes fair_results.json.
    python tests/fair_analyze.py <runs dir>"""
import json, os, sys
import numpy as np

B = sys.argv[1]
CELLS = ["fair_slot_off", "fair_slot_prom", "fair_random_off", "fair_random_prom"]
YEAR = 8766.0


def delta(h0, h1, key, i):
    a, b = h0.get(key), h1.get(key)
    return None if a is None or b is None else b[i] - a[i]


def ratio(n, d):
    return None if n is None or not d else n / d


res = {"criteria_file": "fair_criteria.json", "cells": {}}
for c in CELLS:
    d = os.path.join(B, c)
    st = json.load(open(os.path.join(d, "status.json")))
    h = [json.loads(l) for l in open(os.path.join(d, "history.jsonl"))]
    h = [x for x in h if x.get("score_mean_real") is not None]
    dw = st["darwin"]
    yrs = np.array([x["sim_h"] / YEAR for x in h])
    diff = np.array([x["score_mean_real"] - x["score_mean_ctl"] for x in h])
    half = len(h) // 2
    slope = float(np.polyfit(yrs[half:], diff[half:], 1)[0] * 100) if len(h) - half >= 3 else None
    zr = np.array([x["score_mean_real"] for x in h]); zc = np.array([x["score_mean_ctl"] for x in h])
    h0, h1 = h[half], h[-1]
    dy = (h1["sim_h"] - h0["sim_h"]) / YEAR
    cell = {"sim_years": round(st["time"]["sim_years"], 1), "reports": len(h),
            "step_ms_mean": round(float(np.mean([x["step_ms"] for x in h if x.get("step_ms")])), 1),
            "gpu_temp_max": max((x.get("gpu_temp_c") or 0) for x in h),
            "status": dw["status"], "shifts": dw["trait_response"]["shifts"], "plateaus": dw["trait_response"]["plateaus"],
            "sd_ctl_final": dw["trait_response"]["sd_ctl"], "diff_final": dw["trait_response"]["diff"],
            "late_half": {"z_real_mean": float(zr[half:].mean()), "z_control_mean": float(zc[half:].mean()),
                          "diff_mean": float(diff[half:].mean()), "diff_slope_per_100y": slope, "years": round(dy, 1)},
            "copy_template_z_late": ratio(delta(h0, h1, "diag", 1), delta(h0, h1, "diag", 0)),
            "population_z_late": ratio(delta(h0, h1, "diag", 6), delta(h0, h1, "diag", 5)),
            "fragments_ge12_per_year_late": ratio(delta(h0, h1, "diag", 2), dy),
            "fragments_scored_per_year_late": ratio(delta(h0, h1, "diag", 4), dy),
            "fragment_z_late": ratio(delta(h0, h1, "diag", 3), delta(h0, h1, "diag", 2)),
            "copies_scored_per_year_late": ratio(delta(h0, h1, "diag", 0), dy),
            "monomer_drop_share_late": ratio(delta(h0, h1, "lim", 1), delta(h0, h1, "lim", 0)),
            "slot_drop_share_late": ratio(delta(h0, h1, "lim", 2), delta(h0, h1, "lim", 0)),
            "top1_lineage_share_real_max": max((x.get("top1_share_real") or 0) for x in h),
            "top1_lineage_share_control_max": max((x.get("top1_share_control") or 0) for x in h),
            "clade_top1_share_real_max": max((x.get("clade_top1_share_real") or 0) for x in h),
            "clade_top1_share_control_max": max((x.get("clade_top1_share_control") or 0) for x in h),
            "lineage_sweeps": len(dw["sweeps"]), "dominant_changes": dw["dominant_changes"],
            "clade_events_real": [e for e in dw.get("clade_sweeps", []) if e["half"] == "real"],
            "clade_events_control": sum(1 for e in dw.get("clade_sweeps", []) if e["half"] == "control"),
            "clade_counted": dw.get("clade_counted", 0), "clade_ctl_reports": dw.get("clade_ctl_reports"),
            "rho_above_null99_real_reports": sum(1 for x in h if x.get("rho_wz_real") is not None and x.get("null99_real") is not None and x["rho_wz_real"] > x["null99_real"]),
            "L4_final": dw["real"].get("L4"), "L4_last5": dw["real"].get("L4_last5")}
    if any("scored_share_real" in x for x in h):
        step = max(1, len(h) // 10)
        cell["scored_share_over_time"] = [{"years": round(x["sim_h"] / YEAR, 1), "real": x.get("scored_share_real"),
                                           "control": x.get("scored_share_control")} for x in h[::step] + [h[-1]]]
    res["cells"][c] = cell
base = res["cells"]["fair_slot_off"]
bp = base["plateaus"][-1] if base["plateaus"] else 0.0
for c, cell in res["cells"].items():
    p = cell["plateaus"][-1] if cell["plateaus"] else 0.0
    sd = max(cell["sd_ctl_final"] or 0, base["sd_ctl_final"] or 0)
    cell["criteria"] = {"a": c != "fair_slot_off" and (p - bp) > 2 * sd, "a_detail": {"plateau": p, "baseline_plateau": bp, "threshold": 2 * sd},
                        "b": cell["shifts"] >= 2, "c": cell["clade_counted"] >= 1}
    cell["success"] = bool(cell["criteria"]["a"] or cell["criteria"]["b"] or cell["criteria"]["c"])
json.dump(res, open(os.path.join(B, "fair_results.json"), "w"), indent=1)
for c, cell in res["cells"].items():
    print(c, {k: cell[k] for k in ("sim_years", "status", "shifts", "plateaus", "sd_ctl_final", "copy_template_z_late", "population_z_late",
                                     "fragments_scored_per_year_late", "monomer_drop_share_late", "success")})
    print("   late", cell["late_half"], "crit", cell["criteria"])
