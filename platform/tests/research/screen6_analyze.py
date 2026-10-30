"""Analysis of screen6 (spatial gradient: G1 moving bands + line migration 1e-3 + local pools; G2 same with m = 0;
seeds 1-2) against test_results_2026-10-29/runs/screen6_criteria.json (md5 96b709c1, written before the runs).
    python tests/screen6_analyze.py <runs dir>"""
import json, os, sys
import numpy as np

B = sys.argv[1]
YEAR = 8766.0


def tri(yr, period=50.0):
    ph = (yr / period) % 1.0
    return 2 * ph if ph < 0.5 else 2 - 2 * ph


def corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and a.std() > 0 and b.std() > 0 else None


def cell(name):
    d = os.path.join(B, name)
    st = json.load(open(os.path.join(d, "status.json"))); cfg = json.load(open(os.path.join(d, "config.json")))
    h = [json.loads(l) for l in open(os.path.join(d, "history.jsonl"))]
    h = [x for x in h if x.get("cline_corr_real") is not None and x.get("cline_corr_control") is not None]
    dw = st["darwin"]; rg = dw.get("regional") or {}
    N = int(cfg["sites"]); hN = N // 2
    yrs = np.array([x["sim_h"] / YEAR for x in h]); half = len(h) // 2
    cr = np.array([x["cline_corr_real"] for x in h]); cc = np.array([x["cline_corr_control"] for x in h])
    m10 = [i for i, y in enumerate(yrs) if y >= 10 and h[i].get("top_clade_centroid") is not None]
    cen = [h[i]["top_clade_centroid"] for i in m10]; tr_ = [tri(yrs[i]) for i in m10]
    # per-site monomers vs local fitness and local z (late half, pooled over reports), real and control
    ms_r, zl_r, ms_c, zl_c = [], [], [], []
    for x in h[half:]:
        for s_ in range(N):
            if x["site_zlocal"][s_] is None:
                continue
            (ms_r if s_ < hN else ms_c).append(x["site_Ms"][s_]); (zl_r if s_ < hN else zl_c).append(x["site_zlocal"][s_])
    m1 = [np.mean([v for v in x["site_zlocal"][:hN] if v is not None]) - np.mean([v for v in x["site_zlocal"][hN:] if v is not None]) for x in h]
    return {"sim_years": round(st["time"]["sim_years"], 1), "slots": cfg["slots"], "reports": len(h), "state": st["state"],
            "migration_m": cfg.get("migration_m"), "mig_moves": h[-1].get("mig_moves"),
            "cline_corr_real_late": float(cr[half:].mean()), "cline_corr_control_late": float(cc[half:].mean()),
            "cline_corr_control_late_sd": float(cc[half:].std(ddof=1)),
            "cline_slope_real_late": float(np.mean([x["cline_slope_real"] for x in h[half:] if x.get("cline_slope_real") is not None])),
            "M1_local_late": float(np.mean(m1[half:])), "darwin_status": dw["status"], "shifts": dw["trait_response"]["shifts"],
            "regional": {"counted": rg.get("counted"), "real_events": rg.get("real_events"), "control_events": rg.get("control_events"),
                         "max_window_share_real": max((x.get("regional_max_share_real") or 0) for x in h),
                         "max_window_share_control": max((x.get("regional_max_share_control") or 0) for x in h),
                         "ctl_p95": rg.get("ctl_p95"), "events_real": [e for e in rg.get("events", []) if e["half"] == "real"]},
            "top_clade": {"share_max": max((x.get("top_clade_share") or 0) for x in h), "sites_max": max((x.get("top_clade_sites") or 0) for x in h),
                          "centroid_vs_tri_corr": corr(cen, tr_), "n": len(cen)},
            "global_clade_counted": dw.get("clade_counted", 0), "lineage_sweeps": len(dw["sweeps"]),
            "top1_lineage_real_max": max((x.get("top1_share_real") or 0) for x in h),
            "monomers": {"Ms_real_late_mean": float(np.mean(ms_r)) if ms_r else None, "Ms_control_late_mean": float(np.mean(ms_c)) if ms_c else None,
                         "corr_Ms_zlocal_real": corr(ms_r, zl_r), "corr_Ms_zlocal_control": corr(ms_c, zl_c)},
            "step_ms_mean": float(np.mean([x["step_ms"] for x in h if x.get("step_ms")])), "gpu_temp_max": max((x.get("gpu_temp_c") or 0) for x in h)}


res = {"criteria_file": "screen6_criteria.json (96b709c1)", "cells": {}}
for s in (1, 2):
    for a in ("G1", "G2"):
        res["cells"][f"{a}_s{s}"] = cell(f"screen6_{a}_s{s}")
C = res["cells"]
ev = {}
for s in (1, 2):
    e = {}
    for a in ("G1", "G2"):
        c = C[f"{a}_s{s}"]
        e[f"i_{a}"] = c["cline_corr_real_late"] > c["cline_corr_control_late"] + 2 * c["cline_corr_control_late_sd"]
    e["ii_G1_counted"] = C[f"G1_s{s}"]["regional"]["counted"]; e["ii_G2_counted"] = C[f"G2_s{s}"]["regional"]["counted"]
    e["ii_pass"] = (e["ii_G1_counted"] or 0) >= 2 and (e["ii_G2_counted"] or 0) < (e["ii_G1_counted"] or 0)
    cv = C[f"G1_s{s}"]["top_clade"]["centroid_vs_tri_corr"]
    e["iii_corr"] = cv; e["iii_pass"] = cv is not None and cv < -0.5
    ev[f"s{s}"] = e
res["evaluation"] = ev
res["summary"] = {"i_G1_both": all(ev[k]["i_G1"] for k in ev), "ii_both": all(ev[k]["ii_pass"] for k in ev), "iii_both": all(ev[k]["iii_pass"] for k in ev)}
json.dump(res, open(os.path.join(B, "screen6_results.json"), "w"), indent=1)
for k, c in C.items():
    print(k, {x: c[x] for x in ("sim_years", "slots", "mig_moves", "cline_corr_real_late", "cline_corr_control_late", "cline_corr_control_late_sd", "cline_slope_real_late", "M1_local_late", "shifts")})
    print("   regional", {x: v for x, v in c["regional"].items() if x != "events_real"}, "top", c["top_clade"], "glob", c["global_clade_counted"], "top1", c["top1_lineage_real_max"])
    print("   monomers", c["monomers"], "ms", round(c["step_ms_mean"]), "T", c["gpu_temp_max"])
print(json.dumps(ev, indent=1)); print(res["summary"])
