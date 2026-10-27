"""Stage 4c-3: is "no L2" built in by the layer-C bounds? (coordinator request, 26 Oct 2026)

(1) q(l) table: for the box reference #3493 (scan4c3_box.json) and the 61 L2_loose points of
    the box scan, the mean extension rate r at copies (copy_r_mean) and the duplex exponent of
    the wet and dry phase give q(l) = r / (r + K0 exp((l - 1) e)) for l = 1, 2, 3, 4, 6, and the
    expected number of letters added after a primer of length 2 or 3 on a 24-letter template.
(2) Sensitivity at #3493 (validation fidelity, 16 reps, paired copy_random control): k_t raised
    to 1e6 and 1e7 (layer-C bound 1e5), primers of length 2-6 (PRIMER_LENS), and both.
    python sens4c3.py [out.json]
"""
import json, math, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c3 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402

dev = torch.device("xpu")
guard = thermal.Guard(job="4c3 sensitivity")
FID, REPS = (16, 256, 600), 16


def qtable(r, e):
    q = lambda l: r / (r + chem.K0_OFF * math.exp(min((l - 1) * e, 60.0)))
    out = {f"q{l}": q(l) for l in (1, 2, 3, 4, 6)}
    for l0 in (2, 3):
        surv, add = 1.0, 0.0
        for j in range(l0, 24):
            surv *= q(j)
            add += surv
        out[f"added_after_primer{l0}"] = add
    return out


def stats(o, strict=None, s_m=None, s_c=None):
    cp = o["copies"]
    has = cp > 0
    r = float(o["copy_r_mean"][has].median()) if bool(has.any()) else 0.0
    e = o["dup_e"].median(0).values
    d = {"copying_runs": int(has.sum()), "r_median": r, "e_wet": float(e[0]), "e_dry": float(e[1]),
         "q_wet": qtable(r, float(e[0])) if r > 0 else None, "q_dry": qtable(r, float(e[1])) if r > 0 else None,
         "copy_len_mean": float((o["copy_mean_len"] * cp).sum() / cp.sum().clamp(min=1)),
         "copies_long_share": float(o["copies_long"].sum() / cp.sum().clamp(min=1)),
         "primed_share": float((o["copy_primed_share"] * cp).sum() / cp.sum().clamp(min=1)),
         "L2_loose": int(o["L2_loose"].sum())}
    if strict is not None:
        mk = o["L2_loose"]
        d["L2_strict"] = int(strict.sum())
        d["margin_mean"] = float((s_m - s_c)[mk].mean()) if bool(mk.any()) else None
    return d


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "sens4c3.json")
    res = ckpt.load_json(path) or {}
    box = json.load(open(os.path.join(HERE, "scan4c3_box.json")))
    ref = box["ref"]
    res["ref_point"] = ref["i"]
    # (1) the 61 loose points of the box scan, one run each at screening fidelity
    if "loose_points" not in res:
        X = qmc.Sobol(d=chem.DIM, scramble=True, seed=2026).random(4096).astype(np.float32)
        ix = {n: chem.NAMES.index(n) for n in ("T", "T_m8", "depth", "k_t", "thermo", "x_nuc")}
        Tm = 10 + 80 * X[:, ix["T_m8"]]
        X[:, ix["T"]] = np.clip((Tm - 15 * X[:, ix["T"]]) / 150.0, 0.0, 1.0)
        X[:, ix["depth"]] = 0.5 + 0.5 * X[:, ix["depth"]]
        X[:, ix["k_t"]] = 0.75 + 0.25 * X[:, ix["k_t"]]
        X[:, ix["thermo"]] = 0.375 + 0.625 * X[:, ix["thermo"]]
        X[:, ix["x_nuc"]] = 0.5 + 0.5 * X[:, ix["x_nuc"]]
        pts = sorted({p["i"] for r in box["stage1"].values() for p in r["points"]})
        u = torch.tensor(X[pts], dtype=torch.float32, device=dev)
        g = torch.Generator(device=dev)
        g.manual_seed(901)
        o = chem.simulate(u, g, "screen")
        guard.check()
        per = []
        for k in range(len(pts)):
            r = float(o["copy_r_mean"][k])
            e = o["dup_e"][k]
            per.append({"i": pts[k], "r": r, "e_wet": float(e[0]), "e_dry": float(e[1]),
                        "q_wet": qtable(r, float(e[0])) if r > 0 else None})
        qs = [p["q_wet"] for p in per if p["q_wet"]]
        med = {k: float(np.median([q[k] for q in qs])) for k in qs[0]} if qs else None
        res["loose_points"] = {"n": len(pts), "n_copying": len(qs), "q_wet_median": med, "per_point": per}
        ckpt.atomic_json(res, path)
        print("loose points", len(pts), "q_wet median", med, flush=True)
    # (2) sensitivity at the reference
    u = torch.tensor([ref["u"]] * REPS, dtype=torch.float32, device=dev)
    cfgs = {"base": ({}, (2, 3)), "k_t_1e6": ({"k_t": 1e6}, (2, 3)), "k_t_1e7": ({"k_t": 1e7}, (2, 3)),
            "primers_2_6": ({}, (2, 3, 4, 5, 6)), "k_t_1e7_primers_2_6": ({"k_t": 1e7}, (2, 3, 4, 5, 6))}
    res.setdefault("configs", {})
    for name, (ov, pl) in cfgs.items():
        if name in res["configs"]:
            continue
        chem.PRIMER_LENS = pl
        t = time.time()
        g = torch.Generator(device=dev)
        g.manual_seed(911)
        o = chem.simulate(u, g, FID, override=dict(ov))
        strict, s_m, s_c = chem.paired_control_l2(u, 911, FID, o, override=ov)
        guard.check()
        d = stats(o, strict, s_m, s_c)
        d["seconds"] = time.time() - t
        res["configs"][name] = d
        ckpt.atomic_json(res, path)
        print(name, {k: d[k] for k in ("L2_loose", "L2_strict", "margin_mean", "copy_len_mean", "copies_long_share",
                                       "r_median")}, "q_wet", d["q_wet"], flush=True)
    chem.PRIMER_LENS = (2, 3)
    res["thermal"] = guard.summary()
    ckpt.atomic_json(res, path)


if __name__ == "__main__":
    main()
