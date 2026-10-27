"""Write a REF_L2 file for validate4c3.py (LOT_REF_SCAN) at a chosen point of the box2 scan.

scan4c3_box2.json stores the u vector of its own pick only (#3709); this rebuilds the --box2
point set exactly as scan4c3.py does (65536 Sobol points, box transform, rejection by
lifetime(L_H) / period >= 1, first 4096) and writes {"ref": {"i", "u", "physical", "stage2"}}
for point i. Check: the rebuilt u of the scan's own pick must equal the stored one.
    python make_ref4c3.py 1558 [scan4c3_box2.json] -> scan4c3_box2_ref1558.json
"""
import json, os, sys
import numpy as np
import torch
from scipy.stats import qmc
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c3 as chem  # noqa: E402
import ckpt  # noqa: E402


def box2_points(n1=4096, seed=2026):
    Xb = qmc.Sobol(d=chem.DIM, scramble=True, seed=seed).random(1 << 16).astype(np.float32)
    ix = {n: chem.NAMES.index(n) for n in ("T", "T_m8", "depth", "k_t", "thermo", "x_nuc")}
    Tm = 10 + 80 * Xb[:, ix["T_m8"]]
    T = Tm - 15 * Xb[:, ix["T"]]
    Xb[:, ix["T"]] = np.clip(T / 150.0, 0.0, 1.0)
    dmin = np.clip((Tm - T) / 10.0 + 0.02, 0.0, 1.0)
    Xb[:, ix["depth"]] = np.maximum(0.5 + 0.5 * Xb[:, ix["depth"]], dmin)
    Xb[:, ix["k_t"]] = 0.75 + 0.25 * Xb[:, ix["k_t"]]
    Xb[:, ix["thermo"]] = 0.375 + 0.625 * Xb[:, ix["thermo"]]
    Xb[:, ix["x_nuc"]] = 0.5 + 0.5 * Xb[:, ix["x_nuc"]]
    p = chem.decode(torch.tensor(Xb))
    tab, _, theta, _ = chem.phase_table(p)
    life = 1.0 / chem.turnover_rate(chem.L_H, tab["knuc"][:, 2], p["k_UV"] * p["UV"], 1 / p["tau"], theta, p["p_min"])
    sel = np.nonzero(((life / p["period"]) >= 1.0).numpy())[0][:n1]
    return Xb[sel]


def main():
    i = int(sys.argv[1])
    src = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "scan4c3_box2.json")
    scan = json.load(open(src))
    X = box2_points()
    own = scan["ref"]["i"]
    if not np.allclose(X[own], np.array(scan["ref"]["u"], np.float32), atol=0, rtol=0):
        sys.exit(f"rebuilt u of #{own} differs from the stored one: point set not reproduced")
    x = X[i]
    phys = chem.decode(torch.tensor(x[None, :]))
    ref = {"i": i, "u": [float(v) for v in x], "physical": {k: float(v[0]) for k, v in phys.items()},
           "stage2": scan["stage2"].get(str(i))}
    out = os.path.join(HERE, os.path.basename(src).replace(".json", f"_ref{i}.json"))
    ckpt.atomic_json({"ref": ref, "source": os.path.basename(src)}, out)
    print("wrote", out, {k: round(ref["physical"][k], 3) for k in ("T", "pH", "T_m8", "period", "eps")}, flush=True)


if __name__ == "__main__":
    main()
