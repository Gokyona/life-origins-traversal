"""Box2 stage-2 points: long copies (>= L_H) per template-capable chain per L_H-mer lifetime
(a reproduction number for lineages), from reruns of each stage-2 point (16 reps, validation
fidelity, main run only; same seeds as stage 2). Also lists the layer-A values of points with
>= 8/16 L2_strict.   python s2stats4c3.py [scan json] [out json]
"""
import json, os, sys
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c3 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402
import scan4c3  # noqa: E402

dev = torch.device("xpu")
guard = thermal.Guard(job="4c3 box2 stage-2 stats")


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "scan4c3_box2.json")
    path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "s2stats4c3.json")
    r = json.load(open(src))
    # the box2 point set, as in scan4c3.py --box2
    sys.argv = [sys.argv[0], "--box2"]
    from scipy.stats import qmc
    Xb = qmc.Sobol(d=chem.DIM, scramble=True, seed=2026).random(1 << 16).astype(np.float32)
    ix = {n: chem.NAMES.index(n) for n in ("T", "T_m8", "depth", "k_t", "thermo", "x_nuc")}
    Tm = 10 + 80 * Xb[:, ix["T_m8"]]
    T = Tm - 15 * Xb[:, ix["T"]]
    Xb[:, ix["T"]] = np.clip(T / 150.0, 0.0, 1.0)
    Xb[:, ix["depth"]] = np.maximum(0.5 + 0.5 * Xb[:, ix["depth"]], np.clip((Tm - T) / 10.0 + 0.02, 0.0, 1.0))
    Xb[:, ix["k_t"]] = 0.75 + 0.25 * Xb[:, ix["k_t"]]
    Xb[:, ix["thermo"]] = 0.375 + 0.625 * Xb[:, ix["thermo"]]
    Xb[:, ix["x_nuc"]] = 0.5 + 0.5 * Xb[:, ix["x_nuc"]]
    p = chem.decode(torch.tensor(Xb))
    tab, _, th, _ = chem.phase_table(p)
    life = 1.0 / chem.turnover_rate(chem.L_H, tab["knuc"][:, 2], p["k_UV"] * p["UV"], 1 / p["tau"], th, p["p_min"])
    X = Xb[np.nonzero(((life / p["period"]) >= 1.0).numpy())[0][:4096]]
    out = {"points": {}}
    for k, s2 in r.get("stage2", {}).items():
        i = int(k)
        u = torch.tensor(np.repeat(X[i][None, :], 16, 0), dtype=torch.float32, device=dev)
        g = torch.Generator(device=dev)
        g.manual_seed(80000 + i)
        o = chem.simulate(u, g, scan4c3.FID2)
        guard.check()
        R0 = o["long_copies_per_template_per_lifetime"].numpy()
        ph = chem.decode(torch.tensor(X[i][None, :]))
        d = {"stage2": s2, "R_long_per_template_per_lifetime_pct10_50_90": [float(x) for x in np.percentile(R0, [10, 50, 90])],
             "templates_per_site_mean": float(o["templates_avg_per_site"].mean()),
             "life_over_period": float(life[np.nonzero(((life / p["period"]) >= 1.0).numpy())[0][i]] /
                                       p["period"][np.nonzero(((life / p["period"]) >= 1.0).numpy())[0][i]])}
        if s2["L2_strict"] >= 8:
            d["layer_A"] = {n: float(ph[n][0]) for n in ("T", "T_m8", "depth", "period", "pH", "tau", "k_t", "eps", "L_t")}
        out["points"][k] = d
        ckpt.atomic_json(out, path)
        print(k, {kk: d[kk] for kk in d if kk != "stage2"}, "strict", s2["L2_strict"], "/16 margin", round(s2["margin_mean"], 3), flush=True)
    out["thermal"] = guard.summary()
    ckpt.atomic_json(out, path)


if __name__ == "__main__":
    main()
