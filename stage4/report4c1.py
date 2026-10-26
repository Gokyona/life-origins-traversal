"""Text report of a run4c1_map.py result: python report4c1.py map4c1 > map4c1/report.md"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c1 as chem  # noqa: E402


def reach(curve, target, col):
    cps = [c[0] for c in curve]
    v = [c[col] for c in curve]
    for i, (c, x) in enumerate(zip(cps, v)):
        if x >= target - 1e-12:
            if i == 0:
                return f"<= {c / 1000:.1f}k"
            f = (target - v[i - 1]) / (x - v[i - 1])
            return f"{np.exp(np.log(cps[i - 1]) + f * (np.log(c) - np.log(cps[i - 1]))) / 1000:.1f}k"
    return "not reached"


def main(d):
    r = json.load(open(os.path.join(d, "result.json")))
    out = []
    p = out.append
    p(f"# Stage 4c-1 map ({os.path.basename(d)})\n")
    p(f"config: {r['config']}\n")
    p(f"life fraction (ground truth, P > 0.5): E1 {r['life_fraction']['E1']:.3f}, E2 (= E1 & E2_L8) {r['life_fraction']['E2']:.3f}")
    p(f"GT run flags: {r['gt_flags']}; E2raw & not E1 in {r['gt_nesting_violation_E2raw_not_E1']:.4f} of GT runs; "
      f"slot overflow in {r['gt_overflow_run_share']:.4f} of runs; GT cost {r['gt_ms_per_run']:.1f} ms/run\n")
    for lvl in ("E1", "E2"):
        for col, lab in ((1, "all"), (2, "hard")):
            p(f"\n{lvl} accuracy ({lab}; hard n = {r['methods']['sobol'][0]['curves'][lvl + '_hard_n']})")
            curves = {m: r["methods"][m][0]["curves"][lvl] for m in r["methods"]}
            cps = [c[0] for c in curves["sobol"]]
            p("| method | " + " | ".join(f"{c // 1000}k" for c in cps) + " | reaches sobol@max |")
            p("|---|" + "---|" * (len(cps) + 1))
            target = curves["sobol"][-1][col]
            for m, cv in curves.items():
                p(f"| {m} | " + " | ".join(f"{c[col]:.3f}" for c in cv) + f" | {reach(cv, target, col)} |")
    for m in r["methods"]:
        x = r["methods"][m][0]
        p(f"\n{m}: {x['seconds_this_session']:.0f} s this session" + (f", {x['rounds']} rounds" if "rounds" in x else ""))
    if "maps" in r:
        for lvl in ("E1", "E2"):
            mp = r["maps"][lvl]
            p(f"\n## {lvl} final map (engine data), accuracy {mp['final_map_accuracy']:.3f}")
            imp = sorted(mp["importance"].items(), key=lambda kv: -kv[1])
            p("permutation importance: " + ", ".join(f"{k} {v:.3f}" for k, v in imp))
            p("partial dependence (P at u = 0, 0.25, 0.5, 0.75, 1):")
            for k, v in sorted(mp["partial_dependence"].items(), key=lambda kv: -(max(kv[1]) - min(kv[1]))):
                p(f"  {k:12s} " + " ".join(f"{v[i]:.3f}" for i in (0, 3, 6, 9, 12)) + f"   range {max(v) - min(v):.3f}")
            if lvl == "E2":
                for name, s in mp["slices"].items():
                    P = np.array(s["p"])
                    p(f"slice {name} (rows {s['y']} u = 0..1 step 0.25, cols {s['x']} u = 0..1 step 0.25):")
                    for yi in (0, 5, 10, 15, 20):
                        p("  " + " ".join(f"{P[yi, xi]:.2f}" for xi in (0, 5, 10, 15, 20)))
        c = r.get("composition_E2_positive", {})
        p("\n## composition of E2-positive GT runs")
        p(json.dumps(c, indent=1))
    p(f"\nthermal: {r.get('thermal')}")
    print("\n".join(out))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "map4c1"))
