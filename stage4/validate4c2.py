"""Stage 4c-2 validation tests 2-4 (spec section 3); test 1 is validate4c2_t1.py.

L1 (strict) = RAF >= 3 reactions, persisting >= 5 member turnover times up to the end, with E1
(chem4c2; not nested in E2), and
self-sustaining: the products of the RAF reactions in the best site are more abundant at the end
than in a paired run of the same condition and seed with catalysis switched off (p0 = 0).
(Comparing all members would count the reactants the RAF consumes, which catalysis depletes.)

  (2) catalysis off: p0 = 0 with strong mineral catalysis -> L1 never; minerals still speed up
      chain formation (E2 statistics with versus without mineral catalysis)
  (3) shuffled control: catalysis re-hashed every step -> persistent RAF sets disappear
  (4) hydrolysis link: at the reference catalysis settings, raising T or pH shortens member
      lifetimes and lowers the L1 rate
Reference condition REF: a hot, short-chain regime where persistent RAF sets form (found by a
scan). At REF catalysis depletes the RAF products (they are ligated into longer chains), so the
self-sustaining test fails; L1 does occur where hydrolysis is fast and uncatalysed growth makes
few short chains: REF_L1 = REF at pH 9.5 is the base for tests 3 and 4 (L1 needs E2 in the map's
nested ladder, but these tests look at L1 itself).

  python validate4c2.py [out.json] [tests, e.g. 234]
"""
import json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c2 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402

dev = torch.device("xpu")
guard = thermal.Guard(job="4c2 validation")
FID = (16, 256, 600)          # sites, slots, target steps: small sites for the validation runs
REPS = 8
# k_c is now the uncatalysed rate at 25 C with Ea_unc = 60 kJ/mol: 0.03 / 43.4 gives the same
# uncatalysed rate at 80 C as the k_c = 0.03 (no T dependence) used before; Ea_cat = 0.5 * 60.
# kappa is relative to the uncatalysed rate at 25 C, so kappa = 1e6 keeps the catalysed rate at
# 80 C close to the earlier setting (k_c 0.03, kappa 1e4, Ea_cat 50: 6.9e3 vs 4.6e3 now).
ARR80 = math.exp(60000 / 8.314 * (1 / 298.15 - 1 / 353.15))
REF = dict(T=80.0, pH=7.5, Mg=0.005, salinity=0.25, s_exp=0.0, P=100.0, feed=10.0, tau=1e4, k_c=0.03 / ARR80,
           p0=0.3, kappa=1e6, L_cat=3.0, lam=0.5, Ea_unc=60.0, r_Ea=0.5, depth=0.0, thermo=1.0, mineral_type=0.0,
           x_nuc=1.0, Y_UV=0.1, UV=1.0, Y_H2=0.0, H2=1e-6, k_UV=0.0, avail_time=1e8)
REF_L1 = dict(REF, pH=9.5)


def run(override, seed, reps=REPS, control=True):
    u = torch.full((reps, chem.DIM), 0.5, device=dev)
    g = torch.Generator(device=dev)
    g.manual_seed(seed)
    o = chem.simulate(u, g, FID, record=True, override=dict(override))
    guard.check()
    cand = o["L1_loose"]
    sustained = torch.zeros(reps, dtype=torch.bool)
    ab_cat = torch.zeros(reps)
    ab_ctl = torch.zeros(reps)
    if control:
        sustained, ab_cat, ab_ctl = chem.paired_control(u, seed, FID, o, override=override)
        guard.check()
    member_life = []
    return {"E2": o["E2"].tolist(), "RAF_final": o["RAF_final"].tolist(), "RAF_size": o["RAF_size"].tolist(),
            "RAF_span_turnovers": o["RAF_span_turnovers"].tolist(), "L1_loose": cand.tolist(),
            "L1": sustained.tolist(), "L1_rate": float(sustained.float().mean()), "E1": o["E1"].tolist(),
            "RAF_site_catalysed_share": o["RAF_site_catalysed_share"].tolist(),
            "RAF_product_mass": o["RAF_product_mass"].tolist(),
            "member_abundance_cat": ab_cat.tolist(), "member_abundance_control": ab_ctl.tolist(),
            "cat_events": o["cat_events"].tolist(), "polymer_catalysed_share": o["polymer_catalysed_share"].tolist(),
            "mean_len": _mean_len(o), "frac_sites_L8": o["frac_sites_L8"].tolist(),
            "catalyst_L": o["RAF_catalyst_L"], "t_end_h": o["t_end"].tolist()}


def _mean_len(o):
    h = o["length_hist"][:, 2:].double()
    Ls = torch.arange(2, chem.L_MAX + 1, dtype=torch.float64)
    return ((h * Ls).sum(1) / h.sum(1).clamp(min=1)).tolist()


def member_lifetime(override):
    """Mean lifetime of a 4-mer (3 bonds) against hydrolysis and washout at this condition."""
    u = torch.full((1, chem.DIM), 0.5)
    p = chem.decode(u)
    for k, v in override.items():
        p[k] = torch.as_tensor(v, dtype=torch.float64).expand(1).clone()
    tab, *_ = chem.phase_table(p)
    k = 3 * float(tab["knuc"][0, 2] if override.get("x_nuc", 1.0) >= 0.5 else tab["kpep"][0, 2]) + 1 / float(p["tau"][0])
    return 1 / k


def test_2():
    rows = {}
    rows["REF_p0_0_mineral_clay"] = run(dict(REF, p0=0.0, mineral_type=1.0, mineral=1.0, m_clay=100.0), 201)
    rows["REF_p0_0_no_mineral"] = run(dict(REF, p0=0.0), 202)
    # 64 random conditions with p0 = 0 and strong mineral catalysis
    u = torch.rand(64, chem.DIM, generator=torch.Generator().manual_seed(5)).to(dev)
    g = torch.Generator(device=dev)
    g.manual_seed(203)
    o = chem.simulate(u, g, FID, record=False, override={"p0": 0.0, "mineral_type": 1.0, "mineral": 1.0, "m_clay": 100.0})
    guard.check()
    rnd = {"n": 64, "RAF_final": int(o["RAF_final"].sum()), "L1_loose": int(o["L1_loose"].sum())}
    ok = (not any(rows["REF_p0_0_mineral_clay"]["RAF_final"])) and rnd["RAF_final"] == 0
    speed = np.mean(rows["REF_p0_0_mineral_clay"]["mean_len"]) > np.mean(rows["REF_p0_0_no_mineral"]["mean_len"])
    return {"rows": rows, "random_64_p0_0": rnd, "mineral_lengthens_chains": bool(speed), "pass": bool(ok)}


def test_3():
    fixed = run(REF_L1, 301)
    shuffled = run(dict(REF_L1, salt_per_step=True), 301)
    ok = fixed["L1_rate"] > 0 and shuffled["L1_rate"] < fixed["L1_rate"] and \
        np.mean(shuffled["RAF_span_turnovers"]) < np.mean(fixed["RAF_span_turnovers"])
    return {"fixed": fixed, "shuffled": shuffled, "pass": bool(ok)}


def test_4():
    out = {}
    for name, ov in (("REF", {}), ("T+10", {"T": REF_L1["T"] + 10}), ("T+20", {"T": REF_L1["T"] + 20}),
                     ("pH+1", {"pH": REF_L1["pH"] + 1}), ("pH+2", {"pH": REF_L1["pH"] + 2})):
        cond = dict(REF_L1, **ov)
        r = run(cond, 401)
        r["tetramer_lifetime_h"] = member_lifetime(cond)
        out[name] = r
    life = [out[k]["tetramer_lifetime_h"] for k in ("REF", "T+10", "T+20")]
    lifep = [out[k]["tetramer_lifetime_h"] for k in ("REF", "pH+1", "pH+2")]
    rate_T = [out[k]["L1_rate"] for k in ("REF", "T+10", "T+20")]
    rate_p = [out[k]["L1_rate"] for k in ("REF", "pH+1", "pH+2")]
    ok = (life[0] > life[1] > life[2]) and (lifep[0] > lifep[1] > lifep[2]) and \
        rate_T[2] <= rate_T[0] and rate_p[2] <= rate_p[0] and (rate_T[2] < rate_T[0] or rate_p[2] < rate_p[0])
    out["summary"] = {"lifetime_T": life, "lifetime_pH": lifep, "L1_rate_T": rate_T, "L1_rate_pH": rate_p}
    out["pass"] = bool(ok)
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "validate4c2.json")
    which = sys.argv[2] if len(sys.argv) > 2 else "234"
    res = ckpt.load_json(path) or {}
    res["REF"] = REF
    res["REF_L1"] = REF_L1
    res["fidelity"] = FID
    for k, fn in (("2", test_2), ("3", test_3), ("4", test_4)):
        if k not in which:
            continue
        t = time.time()
        res["test_" + k] = fn()
        res["test_" + k]["seconds"] = time.time() - t
        print("test", k, "pass" if res["test_" + k]["pass"] else "FAIL", f"({time.time() - t:.0f}s)", flush=True)
        ckpt.atomic_json(res, path)
    res["thermal"] = guard.summary()
    ckpt.atomic_json(res, path)


if __name__ == "__main__":
    main()
