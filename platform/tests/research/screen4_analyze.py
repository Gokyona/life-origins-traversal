"""Analysis of screen4 (environmental change: constant / season / trend / catastrophe x seeds 1-2) against
test_results_2026-10-29/runs/screen4_criteria.json (md5 ddbb1f2f, written before the runs).
    python tests/screen4_analyze.py <runs dir>"""
import json, math, os, sys
import numpy as np

B = sys.argv[1]
YEAR = 8766.0
ARMS = ("constant", "season", "trend", "catastrophe")


def win(h, key_fn, lo, hi):
    v = [key_fn(x) for x in h if lo <= x["sim_h"] / YEAR < hi]
    v = [a for a in v if a is not None]
    return float(np.mean(v)) if v else None


def cell(name):
    d = os.path.join(B, name)
    st = json.load(open(os.path.join(d, "status.json")))
    h = [json.loads(l) for l in open(os.path.join(d, "history.jsonl"))]
    h = [x for x in h if x.get("score_mean_real") is not None and x.get("zA_real") is not None and x.get("zA_control") is not None]
    dw = st["darwin"]; tr = dw["trait_response"]
    yrs = np.array([x["sim_h"] / YEAR for x in h])
    m1 = np.array([x["score_mean_real"] - x["score_mean_ctl"] for x in h])
    dA = np.array([x["zA_real"] - x["zA_control"] for x in h]); dB = np.array([x["zB_real"] - x["zB_control"] for x in h])
    gen = np.array([x["generalist_real"] - x["generalist_control"] for x in h])
    genc = np.array([x["generalist_control"] for x in h])
    prefd = np.array([x["pref_real"] - x["pref_control"] for x in h]); w = np.array([x.get("env_w", 0.0) for x in h])
    half = len(h) // 2
    late = slice(half, None)
    # M4 tracking lag: cross-correlation of the pref difference with w(t) after year 10, lags 0..1/3 of the series
    m = yrs >= 10
    lag = None
    if m.sum() > 10 and np.std(w[m]) > 1e-6:
        a, b = prefd[m] - prefd[m].mean(), w[m] - w[m].mean()
        best = max(((float(np.dot(a[k:], b[:len(b) - k]) / (np.linalg.norm(a[k:]) * np.linalg.norm(b[:len(b) - k]) + 1e-12)), k)
                    for k in range(0, max(1, int(m.sum()) // 3))), key=lambda z: z[0])
        dyr = float(np.mean(np.diff(yrs[m])))
        lag = {"corr": best[0], "lag_reports": best[1], "lag_years": best[1] * dyr, "report_interval_years": dyr}
    win10 = []
    for y0 in range(10, int(yrs.max()) - 9, 10):
        v = m1[(yrs >= y0) & (yrs < y0 + 10)]
        if len(v):
            win10.append({"years": [y0, y0 + 10], "m1_mean": float(v.mean())})
    return {"sim_years": round(st["time"]["sim_years"], 1), "reports": len(h), "state": st["state"],
            "darwin_status": dw["status"], "shifts": tr["shifts"], "plateaus": tr["plateaus"], "sd_ctl_final": tr["sd_ctl"],
            "M1_late_mean": float(m1[late].mean()), "M1_final": float(m1[-1]), "M1_10y_windows": win10,
            "M2_late": {"zA_real": float(np.mean([x["zA_real"] for x in h[half:]])), "zB_real": float(np.mean([x["zB_real"] for x in h[half:]])),
                        "zA_control": float(np.mean([x["zA_control"] for x in h[half:]])), "zB_control": float(np.mean([x["zB_control"] for x in h[half:]])),
                        "pref_real_minus_control": float(prefd[late].mean())},
            "M3_generalist_late_mean": float(gen[late].mean()), "M3_generalist_control_late_sd": float(genc[late].std(ddof=1)),
            "M3_generalist_real_late": float(np.mean([x["generalist_real"] for x in h[half:]])),
            "M4_tracking": lag,
            "trend_phase": {"dB_5_10": win(h, lambda x: x["zB_real"] - x["zB_control"], 5, 10), "dB_45_50": win(h, lambda x: x["zB_real"] - x["zB_control"], 45, 50),
                            "dA_50_55": win(h, lambda x: x["zA_real"] - x["zA_control"], 50, 55), "dA_95_100": win(h, lambda x: x["zA_real"] - x["zA_control"], 95, 100),
                            "dA_5_10": win(h, lambda x: x["zA_real"] - x["zA_control"], 5, 10), "dB_95_100": win(h, lambda x: x["zB_real"] - x["zB_control"], 95, 100)},
            "M5": {"top1_lineage_real_max": max(x.get("top1_share_real") or 0 for x in h), "clade_top1_real_max": max(x.get("clade_top1_share_real") or 0 for x in h),
                   "clade_top1_control_max": max(x.get("clade_top1_share_control") or 0 for x in h), "lineage_sweeps": len(dw["sweeps"]),
                   "clade_events_real": sum(1 for e in dw.get("clade_sweeps", []) if e["half"] == "real"),
                   "clade_events_control": sum(1 for e in dw.get("clade_sweeps", []) if e["half"] == "control"),
                   "clade_counted": dw.get("clade_counted", 0), "clade_ctl_p95": dw.get("clade_ctl_p95"),
                   "rho_above_null99_real": sum(1 for x in h if x.get("rho_wz_real") is not None and x.get("null99_real") is not None and x["rho_wz_real"] > x["null99_real"]),
                   "L4_last5": dw["real"].get("L4_last5")},
            "M6": {"cat_events": h[-1].get("cat_events"), "cat_loss_lineage_chains": h[-1].get("cat_loss")},
            "step_ms_mean": float(np.mean([x["step_ms"] for x in h if x.get("step_ms")])), "gpu_temp_max": max(x.get("gpu_temp_c") or 0 for x in h),
            "series": {"years": [round(float(y), 2) for y in yrs], "w": [round(float(v), 3) for v in w], "M1": [round(float(v), 3) for v in m1],
                       "dA": [round(float(v), 3) for v in dA], "dB": [round(float(v), 3) for v in dB], "gen": [round(float(v), 3) for v in gen]}}


res = {"criteria_file": "screen4_criteria.json (ddbb1f2f)", "cells": {}}
for s in (1, 2):
    for a in ARMS:
        res["cells"][f"{a}_s{s}"] = cell(f"screen4_{a}_s{s}")
C = res["cells"]
ev = {}
for s in (1, 2):
    e = {}
    for a in ARMS:
        e[f"1_{a}"] = C[f"{a}_s{s}"]["M5"]["clade_counted"] >= 2
    T = C[f"trend_s{s}"]; sd = T["sd_ctl_final"] or 0
    p = T["trend_phase"]
    e["2a_trend_windows_all_above"] = bool(T["M1_10y_windows"]) and all(wd["m1_mean"] > 2 * sd for wd in T["M1_10y_windows"])
    e["2b_trend_B_rise_0to1"] = None if None in (p["dB_5_10"], p["dB_45_50"]) else p["dB_45_50"] - p["dB_5_10"] > 2 * sd
    e["2b_trend_A_rise_1to0"] = None if None in (p["dA_50_55"], p["dA_95_100"]) else p["dA_95_100"] - p["dA_50_55"] > 2 * sd
    e["2b_note"] = "None = the run did not reach the window (phase not observed)"
    S, K = C[f"season_s{s}"], C[f"constant_s{s}"]
    pooled = math.sqrt((S["M3_generalist_control_late_sd"] ** 2 + K["M3_generalist_control_late_sd"] ** 2) / 2)
    e["3_season_generalist_minus_constant"] = S["M3_generalist_late_mean"] - K["M3_generalist_late_mean"]
    e["3_threshold"] = 2 * pooled
    e["3_season_pass"] = e["3_season_generalist_minus_constant"] > 2 * pooled
    ev[f"s{s}"] = e
res["evaluation"] = ev
res["summary"] = {
    "1_clade_any_arm_both_seeds": {a: ev["s1"][f"1_{a}"] and ev["s2"][f"1_{a}"] for a in ARMS},
    "2_trend_both_seeds": all(ev[k]["2a_trend_windows_all_above"] and ev[k]["2b_trend_B_rise_0to1"] is True and ev[k]["2b_trend_A_rise_1to0"] is True for k in ev),
    "3_season_both_seeds": all(ev[k]["3_season_pass"] for k in ev)}
json.dump(res, open(os.path.join(B, "screen4_results.json"), "w"), indent=1)
for k, c in C.items():
    print(k, {x: c[x] for x in ("sim_years", "shifts", "plateaus", "sd_ctl_final", "M1_late_mean", "M3_generalist_late_mean", "M3_generalist_real_late")},
          c["M2_late"], "lag", c["M4_tracking"], c["M5"], c["M6"], "T", c["gpu_temp_max"], "ms", round(c["step_ms_mean"]))
    if k.startswith("trend"):
        print("   phase", c["trend_phase"], "windows", [round(w_["m1_mean"], 3) for w_ in c["M1_10y_windows"]])
print(json.dumps(ev, indent=1)); print(json.dumps(res["summary"], indent=1))
