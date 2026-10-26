"""Stage 4c-2 validation test 1: the known RAF threshold with the fixed-hash catalysis rule.

Well-mixed, constant concentrations: every binary string of length 1..n is present, the food
set is every string up to length 2, every molecule may catalyse (L_cat = 1, p_cat = p0 for
all lengths), and ligations are reversible (the stage-1e convention). f = p0 * R is the mean
number of reactions a molecule catalyses. Different hash salts play the role of independent
random catalysis assignments. P(RAF) must rise steeply at f about 1-2 (Hordijk, Kauffman &
Steel 2011); the stage-1e reference crossings of P = 0.5 are f* = 1.06 (n = 6), 1.25 (n = 8)
and 1.31 (n = 10), from lot-raf/threshold_check.log.

  python validate4c2_t1.py [out.json]
"""
import json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import raf4c2 as raf  # noqa: E402
import ckpt  # noqa: E402

FGRID = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0]
REF = {6: 1.06, 8: 1.25, 10: 1.31}


def crossing(fs, p):
    for i in range(1, len(fs)):
        if p[i - 1] < 0.5 <= p[i]:
            w = (0.5 - p[i - 1]) / (p[i] - p[i - 1])
            return math.exp(math.log(fs[i - 1]) + w * (math.log(fs[i]) - math.log(fs[i - 1])))
    return float("nan")


def run(n, salts, dev):
    chem, L, seq, ra, rb, rab = (x.to(dev) for x in raf.full_polymer_network(n))
    M, R = L.numel(), ra.numel()
    food = (L <= 2)
    # hash of every (catalyst, reaction) pair: catalyst species c, reaction (a, b)
    c_idx = torch.arange(M, device=dev).repeat_interleave(R)
    r_idx = torch.arange(R, device=dev).repeat(M)
    out = {f: [] for f in FGRID}
    for s in range(salts):
        h = raf.hash_u01([chem[c_idx], L[c_idx], seq[c_idx], chem[ra[r_idx]], L[ra[r_idx]], seq[ra[r_idx]],
                          chem[rb[r_idx]], L[rb[r_idx]], seq[rb[r_idx]]], salt=1000 + s)
        for f in FGRID:
            p0 = f / R
            m = h < p0
            in_raf = raf.max_raf(M, food, ra, rb, rab, r_idx[m], c_idx[m], reversible=True)
            out[f].append(bool(in_raf.any()))
    return {"n": n, "molecules": M, "reactions": R, "salts": salts,
            "P_RAF": [float(np.mean(out[f])) for f in FGRID]}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "validate4c2_t1.json")
    dev = torch.device("xpu")
    res = {"f_grid": FGRID, "rows": []}
    for n, salts in ((6, 400), (8, 200)):
        t = time.time()
        r = run(n, salts, dev)
        r["f_star"] = crossing(FGRID, r["P_RAF"])
        r["f_star_stage1e"] = REF[n]
        r["seconds"] = time.time() - t
        res["rows"].append(r)
        print(n, [round(p, 3) for p in r["P_RAF"]], "f* =", round(r["f_star"], 3), "stage 1e", REF[n],
              f"({r['seconds']:.0f}s)", flush=True)
        ckpt.atomic_json(res, path)
    ok = all(1.0 <= r["f_star"] <= 2.0 or abs(r["f_star"] - r["f_star_stage1e"]) / r["f_star_stage1e"] < 0.15
             for r in res["rows"])
    res["pass"] = bool(ok)
    ckpt.atomic_json(res, path)
    print("pass", ok)


if __name__ == "__main__":
    main()
