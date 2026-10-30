"""Final analysis of a fair_long run against test_results_2026-10-29/runs/fair_long_criteria.json (written before
the runs). Reference: darwin_long1 (same settings, slot-order rationing). CPU only, reads status/history.
    python tests/fair_long_analyze.py <run name> <out.json>"""
import json, math, os, sys
import numpy as np

YEAR = 8766.0
G = os.path.expanduser("~/.gokyona-env/runs")
name, out = sys.argv[1], sys.argv[2]


def load(n):
    h = [json.loads(l) for l in open(f"{G}/{n}/history.jsonl")]
    return json.load(open(f"{G}/{n}/status.json")), [x for x in h if x.get("score_mean_real") is not None]


def window(h, Y):
    w = [x for x in h if abs(x["sim_h"] / YEAR - Y) <= 25]
    if len(w) < 2:
        return None
    return {"n": len(w), "diff_mean": float(np.mean([x["score_mean_real"] - x["score_mean_ctl"] for x in w])),
            "ctl_sd": float(np.std([x["score_mean_ctl"] for x in w], ddof=1))}


def delta(h0, h1, key, i):
    a, b = h0.get(key), h1.get(key)
    return None if a is None or b is None else b[i] - a[i]


def ratio(n, d):
    return None if n is None or not d else n / d


st, h = load(name)
_, href = load("darwin_long1")
dw = st["darwin"]; tr = dw["trait_response"]
yrs = np.array([x["sim_h"] / YEAR for x in h]); diff = np.array([x["score_mean_real"] - x["score_mean_ctl"] for x in h])
half = len(h) // 2
slope = float(np.polyfit(yrs[half:], diff[half:], 1)[0] * 100)
h0, h1 = h[half], h[-1]
dy = (h1["sim_h"] - h0["sim_h"]) / YEAR
blocks = {}
for y, d in zip(yrs, diff):
    blocks.setdefault(int(y // 250) * 250, []).append(d)
res = {"run": name, "sim_years": round(st["time"]["sim_years"], 1), "wall_s": round(st["time"]["wall_s"]), "reports": len(h),
       "state": st["state"], "status": dw["status"], "shifts": tr["shifts"], "plateaus": tr["plateaus"], "sd_ctl_final": tr["sd_ctl"],
       "diff_final": tr["diff"], "block_means_250y": {k: round(float(np.mean(v)), 3) for k, v in sorted(blocks.items())},
       "late_half_slope_per_100y": slope, "late_half_diff_mean": float(diff[half:].mean()),
       "copy_template_z_late": ratio(delta(h0, h1, "diag", 1), delta(h0, h1, "diag", 0)),
       "population_z_late": ratio(delta(h0, h1, "diag", 6), delta(h0, h1, "diag", 5)),
       "fragments_ge12_per_year_late": ratio(delta(h0, h1, "diag", 2), dy), "fragment_z_late": ratio(delta(h0, h1, "diag", 3), delta(h0, h1, "diag", 2)),
       "monomer_drop_share_late": ratio(delta(h0, h1, "lim", 1), delta(h0, h1, "lim", 0)),
       "slot_drop_share_late": ratio(delta(h0, h1, "lim", 2), delta(h0, h1, "lim", 0)),
       "top1_lineage_share_real_max": max(x.get("top1_share_real") or 0 for x in h),
       "top1_lineage_sites_real_max": max(x.get("top1_sites_real") or 0 for x in h),
       "clade_top1_share_real_max": max(x.get("clade_top1_share_real") or 0 for x in h),
       "clade_top1_share_control_max": max(x.get("clade_top1_share_control") or 0 for x in h),
       "lineage_sweeps": len(dw["sweeps"]), "dominant_changes": dw["dominant_changes"],
       "clade_events_real": [e for e in dw.get("clade_sweeps", []) if e["half"] == "real"],
       "clade_events_control": sum(1 for e in dw.get("clade_sweeps", []) if e["half"] == "control"),
       "clade_counted": dw.get("clade_counted", 0),
       "rho_above_null99_real_reports": sum(1 for x in h if x.get("rho_wz_real") is not None and x.get("null99_real") is not None and x["rho_wz_real"] > x["null99_real"]),
       "L4_final": dw["real"].get("L4"), "gpu_temp_max": max(x.get("gpu_temp_c") or 0 for x in h),
       "step_ms_mean": float(np.mean([x["step_ms"] for x in h if x.get("step_ms")]))}
crit = {}
p = tr["plateaus"][-1] if tr["plateaus"] else 0.0
ref_sd = json.load(open(f"{G}/darwin_long1/status.json"))["darwin"]["trait_response"]["sd_ctl"]
thr = 2 * max(tr["sd_ctl"] or 0, ref_sd or 0)
crit["a"] = {"pass": p - 0.569 > thr, "plateau": p, "ref_plateau": 0.569, "threshold": thr, "sd_ctl_run": tr["sd_ctl"], "sd_ctl_ref": ref_sd}
crit["b"] = {"pass": tr["shifts"] >= 2, "shifts": tr["shifts"]}
crit["c"] = {"pass": res["clade_counted"] >= 1, "counted": res["clade_counted"]}
crit["d"] = {}
for Y in (1000, 2000):
    w, r = window(h, Y), window(href, Y)
    if w is None or r is None:
        crit["d"][str(Y)] = {"evaluated": False}
        continue
    pooled = math.sqrt((r["ctl_sd"] ** 2 + w["ctl_sd"] ** 2) / 2)
    crit["d"][str(Y)] = {"evaluated": True, "run": w, "ref": r, "pooled_sd": pooled, "margin": w["diff_mean"] - r["diff_mean"],
                         "pass": w["diff_mean"] - r["diff_mean"] > 2 * pooled}
crit["d"]["pass"] = all(v["pass"] for k, v in crit["d"].items() if isinstance(v, dict) and v.get("evaluated"))
crit["e"] = {"pass": slope > 0, "slope_per_100y": slope}
res["criteria"] = crit
json.dump(res, open(out, "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("clade_events_real",)}, indent=1))
