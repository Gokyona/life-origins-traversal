"""P3-1 validation tests 7 and 8 (simulation positive control), called by validate_p31.py.

No natural L2 reference exists for the P3-1 chemistry before the maps, so, as in Stage 4 test 4, an explicit
positive control is used. G/C pair (fast copying, duplexes that melt only when hot), base 5 C, hot pulses of
5 min at 100 C every 0.5 h (16-mer duplexes melt in each pulse, are stable between pulses), a master template
of 16 letters (8 per site), primers of 5 letters (PRIMER_LENS = (5,), outside the map setting (2, 3): with
measured stacks a 2-3-mer primer adds one letter only in ~1/150 starts at 5 C, see the P3-1 report) that
match the first 5 letters of the master and of its complement, topped up to 100 each per site every step
(clamp_primers, validation fixture), M* clamped, no growth or initiation, no hydrolysis, turnover by washout (tau 5 h; with hydrolysis
k_h 0.02 /h, tried first, fragments of copies that keep the primer prefix dominated the L >= 6 chain mass and
kept S near 1 with random letters), run 80 h, 16 sites x 1024 slots, 1200 target steps. A copy adds L_copy = 11 letters (a 7-letter copied region, tried first, has only 128
variants, so random-letter controls kept S near 1).
  (7) eps_0 = eps_1 = x / L_copy, x = eps * L_copy in XS; 16 reps, seed 701, paired copy_random controls
      (validate4c3.run machinery: c2 where c1 alone would pass). Criterion (Stage 4 test 9b, pre-registered,
      written to the output before any run): (a) S (mean L2_share_last5) non-increasing in eps, one upward
      step <= 0.03 tolerated; (b) L2_strict >= 12/16 at x <= 0.3; (c) L2_strict <= 2/16 at x >= 1.5.
  (8) the same condition with random letters (eps_0 = eps_1 = 0.5), 64 reps, seed 801, both copy_random
      controls for every candidate; pass: L2_strict <= 2 of 64.
"""
import math, time
import torch
import chem_p31 as chem
import validate_p31 as V

MASTER = 0b1011001110001011
L_COPY = 11
XS = (0.15, 0.3, 0.5, 0.7, 1.0, 1.5, 3.0)
FID7 = (16, 1024, 1200)
CRITERION = ("(a) S non-increasing in eps (one upward step <= 0.03 tolerated); (b) L2_strict >= 12/16 at eps*L_copy "
             "<= 0.3; (c) L2_strict <= 2/16 at eps*L_copy >= 1.5; L_copy = 11, eps_0 = eps_1 = x / 7, x in " + str(list(XS)))


def fixture(eps):
    return dict(V.BASE, T=5.0, T_p=100.0, tau_p=5 / 60, f_p=(5 / 60) / 0.5, thermo=1.0, eps_0=eps, eps_1=eps,
                k_h_fixed=0.0, tau=5.0, clamp_mstar=1000.0, t_end=80.0, L_t=16.0,
                init_templates={"n": [8], "L": [16], "seq": [MASTER]},
                clamp_primers={"n": [100, 100], "L": [5, 5], "seq": [(~MASTER) & 31, MASTER & 31]})


def run(ov, seed, reps, both=False):
    chem.PAIR = "GC"
    old = chem.PRIMER_LENS
    chem.PRIMER_LENS = (5,)
    try:
        u = torch.full((reps, chem.DIM_P31), 0.5, device=V.dev)
        o = chem.simulate(u, V.gen(seed), FID7, override=dict(ov))
        V.G().check()
        strict, s_m, s_c = chem.paired_control_l2(u, seed, FID7, o, override=dict(ov), both=both)
        V.G().check()
    finally:
        chem.PRIMER_LENS = old
        chem.PAIR = None
    f = lambda k: [round(float(x), 5) for x in o[k]]
    nn = lambda t: [None if math.isnan(float(x)) else round(float(x), 5) for x in t]
    return {"L2_strict": int(strict.sum()), "L2_loose": int(o["L2_loose"].sum()), "S": f("L2_share_last5"),
            "S_main_rerun": nn(s_m), "S_ctl": nn(s_c), "copies": f("copies"), "copy_len": f("copy_mean_len"),
            "copy_err": f("copy_error_rate"), "R": f("long_copies_per_template_per_lifetime"),
            "templates_per_site": f("templates_avg_per_site"), "tau_LH_h": float(o["L2_tau_h"][0]),
            "steps": int(o["steps"][0]), "reached_end": int(o["reached_end"].sum())}


def judge(rows):
    S = [r["S_mean"] for r in rows]
    ups = [S[k + 1] - S[k] for k in range(len(S) - 1) if S[k + 1] > S[k]]
    a = len(ups) == 0 or (len(ups) == 1 and ups[0] <= 0.03)
    b = all(r["L2_strict"] >= 12 for r in rows if r["x"] <= 0.3)
    c = all(r["L2_strict"] <= 2 for r in rows if r["x"] >= 1.5)
    return {"a_S_nonincreasing": a, "b_strict_ge12_low": b, "c_strict_le2_high": c, "pass": bool(a and b and c)}


def test_7():
    out = {"criterion_preregistered": CRITERION, "L_copy": L_COPY, "fidelity": list(FID7), "rows": []}
    for x in XS:
        t = time.time()
        r = run(fixture(x / L_COPY), 701, V.REPS)
        row = {"x": x, "eps": x / L_COPY, "L2_strict": r["L2_strict"], "L2_loose": r["L2_loose"],
               "S_mean": sum(r["S"]) / len(r["S"]), "copy_len_mean": sum(r["copy_len"]) / len(r["copy_len"]),
               "copy_err_mean": sum(r["copy_err"]) / len(r["copy_err"]), "R_mean": sum(r["R"]) / len(r["R"]),
               "seconds": round(time.time() - t, 1), "detail": r}
        out["rows"].append(row)
        print("test 7", {k: v for k, v in row.items() if k != "detail"}, flush=True)
    out["verdict"] = judge(out["rows"])
    out["pass"] = out["verdict"]["pass"]
    return out


def test_8():
    r = run(fixture(0.5), 801, 64, both=True)
    return {"eps": 0.5, "reps": 64, "result": r, "pass": r["L2_strict"] <= 2}
