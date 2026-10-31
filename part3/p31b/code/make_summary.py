"""Builds results/l2_summary.csv (per variant x point means over the 16 runs) from results/l2/runs.jsonl."""
import csv, json, os, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, "..", "results", "l2", "runs.jsonl")
L = [json.loads(l) for l in open(src)]
keys = ["pc_L2_strict", "L2_loose", "copies", "copy_tries", "copy_mean_len", "L2_share_last5", "copy_conc_gini", "copy_conc_top10",
        "copy_conc_bins", "break_events", "paired_break_events", "round_cut_copies_lost", "ms_melt_wait_share",
        "ms_copies_drawn_per_template_step", "ms_limited_slots_share", "reached_end", "steps"]
rows = []
for d in L:
    pts = d["points"]
    for n in dict.fromkeys(pts):
        ix = [i for i, p in enumerate(pts) if p == n]
        row = {"variant": d["variant"], "point": n, "runs": len(ix), "seconds_batch": round(d["seconds"], 1)}
        for k in keys:
            if k in d["runs"]:
                x = np.array([np.nan if d["runs"][k][i] is None else d["runs"][k][i] for i in ix], float)
                row[k] = int(np.nansum(x > 0.5)) if k == "pc_L2_strict" else (round(float(np.nanmean(x)), 5) if np.isfinite(x).any() else "")
            else:
                row[k] = ""
        rows.append(row)
out = os.path.join(here, "..", "results", "l2_summary.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["variant", "point", "runs", "seconds_batch"] + keys)
    w.writeheader()
    w.writerows(rows)
print("wrote", out, len(rows), "rows")
