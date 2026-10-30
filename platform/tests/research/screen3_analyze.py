"""Analysis of screen3 (random rationing x {none, domain_break_guard, promoter_required}, seeds 1-2) against
test_results_2026-10-29/runs/screen3_criteria.json (md5 db332d7b, written before the runs).
    python tests/screen3_analyze.py <runs dir>"""
import json, math, os, sys
import numpy as np

B = sys.argv[1]
YEAR = 8766.0
ARMS = ("random", "guard", "p1")


def ratio(n, d):
    return None if n is None or not d else n / d


def cell(name):
    d = os.path.join(B, name)
    st = json.load(open(os.path.join(d, "status.json")))
    h = [json.loads(l) for l in open(os.path.join(d, "history.jsonl"))]
    h = [x for x in h if x.get("score_mean_real") is not None and x.get("z_all_real") is not None]
    dw = st["darwin"]; tr = dw["trait_response"]
    yrs = np.array([x["sim_h"] / YEAR for x in h])
    prim = np.array([x["z_all_real"] - x["z_all_control"] for x in h])
    zc_all = np.array([x["z_all_control"] for x in h])
    half = len(h) // 2
    h0, h1 = h[half], h[-1]
    dy = (h1["sim_h"] - h0["sim_h"]) / YEAR
    dl = lambda k, i: None if h0.get(k) is None or h1.get(k) is None else h1[k][i] - h0[k][i]
    after10 = [x for x in h if x["sim_h"] / YEAR >= 10]
    sh_r = [x.get("scored_share_real") for x in after10]; sh_c = [x.get("scored_share_control") for x in after10]
    return {"sim_years": round(st["time"]["sim_years"], 1), "reports": len(h), "state": st["state"],
            "primary_late_mean": float(prim[half:].mean()), "primary_late_slope_per_100y": float(np.polyfit(yrs[half:], prim[half:], 1)[0] * 100),
            "control_all_z_late_sd": float(zc_all[half:].std(ddof=1)),
            "z_all_real_late": float(np.mean([x["z_all_real"] for x in h[half:]])), "z_all_control_late": float(np.mean([x["z_all_control"] for x in h[half:]])),
            "z_scored_real_late": float(np.mean([x["z_scored_real"] for x in h[half:] if x.get("z_scored_real") is not None])),
            "z_scored_control_late": float(np.mean([x["z_scored_control"] for x in h[half:] if x.get("z_scored_control") is not None])),
            "scored_share_real_min_after10": min(sh_r) if sh_r else None, "scored_share_control_min_after10": min(sh_c) if sh_c else None,
            "scored_share_real_final": h1.get("scored_share_real"), "scored_share_control_final": h1.get("scored_share_control"),
            "reports_scored_share_below_0p8": sum(1 for a, b in zip(sh_r, sh_c) if (a is not None and a < 0.8) or (b is not None and b < 0.8)),
            "darwin_status": dw["status"], "shifts": tr["shifts"], "plateaus": tr["plateaus"], "darwin_diff_final": tr["diff"], "sd_ctl_final": tr["sd_ctl"],
            "copy_template_z_late": ratio(dl("diag", 1), dl("diag", 0)), "population_scored_z_late": ratio(dl("diag", 6), dl("diag", 5)),
            "copies_scored_per_year_late": ratio(dl("diag", 0), dy),
            "fragments_ge12_per_year_late": ratio(dl("diag", 2), dy), "fragments_scored_per_year_late": ratio(dl("diag", 4), dy),
            "breaks_per_year_late": ratio(dl("breaks", 0), dy), "suppressed_breaks_per_year_late": ratio(dl("breaks", 1), dy),
            "suppressed_share_of_breaks": ratio(dl("breaks", 1), (dl("breaks", 0) or 0) + (dl("breaks", 1) or 0)),
            "monomer_drop_share_late": ratio(dl("lim", 1), dl("lim", 0)), "slot_drop_share_late": ratio(dl("lim", 2), dl("lim", 0)),
            "top1_lineage_real_max": max(x.get("top1_share_real") or 0 for x in h), "clade_top1_real_max": max(x.get("clade_top1_share_real") or 0 for x in h),
            "clade_top1_control_max": max(x.get("clade_top1_share_control") or 0 for x in h),
            "lineage_sweeps": len(dw["sweeps"]), "clade_events_real": sum(1 for e in dw.get("clade_sweeps", []) if e["half"] == "real"),
            "clade_counted": dw.get("clade_counted", 0),
            "rho_above_null99_real": sum(1 for x in h if x.get("rho_wz_real") is not None and x.get("null99_real") is not None and x["rho_wz_real"] > x["null99_real"]),
            "L4_last5": dw["real"].get("L4_last5"),
            "step_ms_mean": float(np.mean([x["step_ms"] for x in h if x.get("step_ms")])), "gpu_temp_max": max(x.get("gpu_temp_c") or 0 for x in h)}


res = {"criteria_file": "screen3_criteria.json (db332d7b)", "cells": {}}
for s in (1, 2):
    for a in ARMS:
        res["cells"][f"{a}_s{s}"] = cell(f"screen3_{a}_s{s}")
C = res["cells"]
ev = {}
for s in (1, 2):
    A, G, P = C[f"random_s{s}"], C[f"guard_s{s}"], C[f"p1_s{s}"]
    pooled = lambda X, Y: math.sqrt((X["control_all_z_late_sd"] ** 2 + Y["control_all_z_late_sd"] ** 2) / 2)
    ev[f"s{s}"] = {
        "1_P_minus_A": P["primary_late_mean"] - A["primary_late_mean"], "1_threshold": 2 * pooled(P, A),
        "1_pass": P["primary_late_mean"] - A["primary_late_mean"] > 2 * pooled(P, A),
        "2_pass": P["primary_late_slope_per_100y"] > 0 and P["primary_late_slope_per_100y"] > A["primary_late_slope_per_100y"],
        "2_slopes": [P["primary_late_slope_per_100y"], A["primary_late_slope_per_100y"]],
        "3_pass": P["shifts"] >= 2 and A["shifts"] < 2, "3_shifts_P_A": [P["shifts"], A["shifts"]],
        "4_pass": P["clade_counted"] >= 1,
        "validity_P": (P["scored_share_real_min_after10"] or 0) >= 0.7 and (P["scored_share_control_min_after10"] or 0) >= 0.7,
        "G_minus_A": G["primary_late_mean"] - A["primary_late_mean"], "G_threshold": 2 * pooled(G, A),
        "G_negative_control_pass": abs(G["primary_late_mean"] - A["primary_late_mean"]) < 2 * pooled(G, A)}
res["evaluation"] = ev
res["P_positive"] = all(ev[k]["1_pass"] for k in ev) and all(ev[k]["validity_P"] for k in ev)
res["G_negative_control_pass_both"] = all(ev[k]["G_negative_control_pass"] for k in ev)
json.dump(res, open(os.path.join(B, "screen3_results.json"), "w"), indent=1)
for k, c in C.items():
    print(k, {x: (round(v, 4) if isinstance(v, float) else v) for x, v in c.items()})
print(json.dumps(ev, indent=1))
print("P_positive", res["P_positive"], "G_negative_control_pass_both", res["G_negative_control_pass_both"])
