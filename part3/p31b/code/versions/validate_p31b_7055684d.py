"""Part III step P3-1b, section 1 validation tests 1-4 (spec docs/part3_p31b_spec.md, e02a9fe) for chem_p31b.py.
Fixtures and criteria: IMPLEMENTATION_P31B.md (written before any result). Test 5 is kill_resume_p31b.sh, test 6
preflight_p31b.py.
    python validate_p31b.py out.json [tests, e.g. 1234]
"""
import json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_p31b as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402
import validate_p31 as V1  # noqa: E402  (P3-1 helpers: _box2_points, _flat, _same, BASE)

dev = torch.device("xpu")
guard = None


def G():
    global guard
    if guard is None:
        guard = thermal.Guard(job="P3-1b validation")
    return guard


def gen(seed):
    g = torch.Generator(device=dev)
    g.manual_seed(seed)
    return g


def md5(f):
    return os.popen(f"md5sum {os.path.join(HERE, f)}").read().split()[0]


# ---------- test 1: off switch ----------
def test_1():
    import pickle
    import chem_p31
    res = {"chem_p31_md5": md5("chem_p31.py"), "chem_p31b_md5": md5("chem_p31b.py")}
    chem.MULTI_START = False
    # (a) the 16 stored box2 runs (P3-1 test 4)
    ref = pickle.load(open(os.path.join(HERE, "regress_4c3_stored.pkl"), "rb"))
    idx = [1558, 3709] + list(range(14))
    u = torch.tensor(V1._box2_points()[idx], dtype=torch.float32, device=dev)
    fid, seed = (16, 256, 600), 4401

    def run(mod, uu, fid_, sd):
        o = mod.simulate(uu, gen(sd), fid_, record=True)
        G().check()
        s1 = mod.paired_control(uu, sd, fid_, o)
        s2 = mod.paired_control_l2(uu, sd, fid_, o, both=True)
        G().check()
        d = V1._flat(o)
        for k, v in zip(("L1_strict", "ab_cat", "ab_ctl"), s1):
            d["pc_" + k] = v.cpu()
        for k, v in zip(("L2_strict", "S_main", "S_ctl"), s2):
            d["pc_" + k] = v.cpu()
        return d
    chem.PAIR = None
    t0 = time.time()
    new = run(chem, u, fid, seed)
    bad_a = V1._same(ref, new)
    res["a"] = {"runs": 16, "fields": len(ref), "differing": bad_a, "copies": float(ref["copies"].sum()),
                "L2_strict": int(ref["pc_L2_strict"].sum()), "seconds": round(time.time() - t0, 1)}
    # (b) 8 AS2U points that copy, chem_p31 vs chem_p31b
    from scipy.stats import qmc
    pts = [287, 16, 324, 552, 31, 184, 328, 303]
    X = qmc.Sobol(d=chem.DIM_P31, scramble=True, seed=3101).random(640).astype(np.float32)[pts]
    ub = torch.tensor(X, dtype=torch.float32, device=dev)
    chem.PAIR = chem_p31.PAIR = "AS2U"
    t0 = time.time()
    r_old = run(chem_p31, ub, "screen", 3110)
    r_new = run(chem, ub, "screen", 3110)
    bad_b = V1._same(r_old, r_new)
    res["b"] = {"points": pts, "fields": len(r_old), "differing": bad_b, "copies": r_old["copies"].tolist(),
                "seconds": round(time.time() - t0, 1)}
    res["pass"] = not bad_a and not bad_b
    chem.MULTI_START = True
    return res


# ---------- fixtures ----------
def p31_fixture(extra):
    ov = dict(V1.BASE)
    ov.update(extra)
    return ov


def gc_rates(c):
    """GC: (r, required letter) for template letter 0 (G, requires C) and letter 1 (C, requires G), /h, at c mM."""
    km = chem.PAIRS["GC"]["km_kmax"]
    rC = km[1][1] * c / (km[1][0] + c)        # C*C on G (template letter 0)
    rG = km[0][1] * c / (km[0][0] + c)        # G*G on C (template letter 1)
    return rC, rG


def one_step(pair, u_dim, B, fid, ov, seed, multi, copy_log=False, state=False):
    chem.PAIR = pair
    chem.MULTI_START = multi
    ov = dict(ov)
    if copy_log:
        ov["copy_log"] = True
    if state:
        ov["return_state"] = True
    u = torch.full((B, u_dim), 0.5, dtype=torch.float32, device=dev)
    o = chem.simulate(u, gen(seed), fid, override=ov)
    G().check()
    return o


def _cmp(on, off, per):
    """mean copies per template-step, on/off ratio and z from per-run counts."""
    a, b = on.double() / per, off.double() / per
    ma, mb = float(a.mean()), float(b.mean())
    se = math.sqrt(float(a.var()) / a.numel() + float(b.var()) / b.numel())
    return {"mean_on": ma, "mean_off": mb, "ratio": ma / mb if mb else float("nan"),
            "z": (ma - mb) / se if se else float("nan")}


# ---------- test 2: low-rate limit ----------
def test_2():
    """Review 2: statistic = copy tries (after the start, primer, >= 2-letter and block filters, before the slot step);
    post-slot copies reported, not gated. 2b dropped (R1.1, not a low-rate case). 2c low-rate (R2.c). 2d (R2.d)."""
    res = {"2b_dropped_record": {"53711ca2": {"ratio": 0.8762, "z": -8.85}, "21085a25": {"ratio": 0.1989, "z": -78.3},
                                 "note": "dropped by Review 1: lambda(t_copy+t_melt) = 1.5 is not a low-rate case"}}
    fid = (64, 64, 1)
    nt = 32
    per = fid[0] * nt
    Ms = 780.0                                   # c = 0.39 mM = K_m of C*C: r = 19.5 /h (template letter 0)
    c = Ms / 2000
    r = gc_rates(c)[0]
    kf = 100 * r
    q = r / (r + kf)
    lam = r * q
    dt = 0.02 / lam
    ov2a = p31_fixture({"clamp_mstar": Ms, "koff_fixed": kf, "t_end": dt,
                        "init_templates": {"n": nt, "L": 12, "seq": 0}})
    E = 1 + sum(q ** m for m in range(1, 11))
    lt = lam * (E / r + 1 / kf)
    B = 512
    on = one_step("GC", chem.DIM_P31, B, fid, ov2a, 7001, True)
    off = one_step("GC", chem.DIM_P31, B, fid, ov2a, 7002, False)
    d = _cmp(on["copy_tries"], off["copy_tries"], per)
    d.update({"post_slot": _cmp(on["copies"], off["copies"], per), "lam_dt": 0.02, "lam_tcopy_plus_tmelt": lt, "B": B,
              "round_cut_copies_lost_off": float(off["round_cut_copies_lost"].sum()),
              "expected_ratio_renewal_vs_old": (0.02 / (1 + lt)) / (1 - math.exp(-0.02))})
    d["pass"] = abs(d["ratio"] - 1) <= 0.03 and abs(d["z"]) < 3
    res["2a_P31_start_limited"] = d
    print("2a", d, flush=True)
    # 2c: 4c-3 rule (PAIR = None), start rate = extension rate r = k_c k_t M*/V. Low rate (R2.c): T = T_m8 + 10 C
    # (k_off rises with length above T_m8) and r = 1e4 /h, so q(1) = r / (r + K0) ~ 0.01 and lambda(t_copy+t_melt) ~ 0.01
    # (raising T alone cannot go below q(1) = 1e5 / (1e5 + K0) = 0.09 because k_off(1) = K0 at any T). Copy tries are
    # ~2e-4 per template-step, so 32 batches of 2048 runs x 64 sites are pooled (memory: [nT, 64] per-template draws).
    r4 = 1e4
    ov = {"T": 35.0, "T_m8": 25.0, "k_c": 1.0, "k_t": r4, "thermo": 1.0, "depth": 0.0, "mineral_type": 0.0,
          "x_nuc": 1.0, "init_factor": 0.0, "k_h_fixed": 0.0, "k_UV": 0.0, "tau": 1e12, "eps": 1e-3, "p0": 0.0,
          "avail_time": 1e8, "P": 100.0, "feed": 10.0, "L_t": 12.0, "pH": 7.5, "clamp_mstar": 1000.0,
          "t_end": 0.02 / r4, "init_templates": {"n": nt, "L": 12, "seq": 0}}
    fid_c = (64, 64, 1)
    per_c = fid_c[0] * nt
    ons, offs, ons_p, offs_p = [], [], [], []
    for k in range(32):
        o1 = one_step(None, chem.DIM, 2048, fid_c, ov, 7003 + 10 * k, True)
        o0 = one_step(None, chem.DIM, 2048, fid_c, ov, 7004 + 10 * k, False)
        ons.append(o1["copy_tries"]); offs.append(o0["copy_tries"])
        ons_p.append(o1["copies"]); offs_p.append(o0["copies"])
    d = _cmp(torch.cat(ons), torch.cat(offs), per_c)
    d.update({"post_slot": _cmp(torch.cat(ons_p), torch.cat(offs_p), per_c), "B": 32 * 2048, "sites": fid_c[0],
              "tries_on": float(torch.cat(ons).sum()), "tries_off": float(torch.cat(offs).sum())})
    # lambda (t_copy + t_melt) from the model's own k_off (dH_bp, T_m8, T as above)
    e = (25 + 20 * 0.5) * 1000 / 8.314 * (1 / 298.15 - 1 / 308.15)          # dH_bp at u = 0.5: 35 kJ/mol
    ko = [chem.K0_OFF * math.exp((l - 1) * e) for l in range(1, 13)]
    a_, E4 = 1.0, 0.0
    for l in range(1, 12):
        a_ *= r4 / (r4 + ko[l - 1])
        E4 += a_
    lt4 = r4 * (E4 / r4 + 1 / ko[11])
    d.update({"lam_dt": 0.02, "lam_tcopy_plus_tmelt": lt4,
              "expected_ratio_renewal_vs_old": (0.02 / (1 + lt4)) / (1 - math.exp(-0.02))})
    d["pass"] = abs(d["ratio"] - 1) <= 0.03 and abs(d["z"]) < 3
    res["2c_4c3_rule"] = d
    print("2c", d, flush=True)
    # 2d (R2.d): fixture 2a, SLOT_ROUND_FIX on in both paths, post-slot copies ON vs OFF
    chem.SLOT_ROUND_FIX = True
    try:
        on = one_step("GC", chem.DIM_P31, B, fid, ov2a, 7011, True)
        off = one_step("GC", chem.DIM_P31, B, fid, ov2a, 7012, False)
    finally:
        chem.SLOT_ROUND_FIX = False
    d = _cmp(on["copies"], off["copies"], per)
    d.update({"tries": _cmp(on["copy_tries"], off["copy_tries"], per), "B": B,
              "round_cut_copies_lost_off": float(off["round_cut_copies_lost"].sum())})
    d["pass"] = abs(d["ratio"] - 1) <= 0.03 and abs(d["z"]) < 3
    res["2d_slot_round_fix"] = d
    print("2d", d, flush=True)
    res["pass"] = all(v["pass"] for k, v in res.items() if isinstance(v, dict) and "pass" in v)
    chem.MULTI_START = True
    return res


# ---------- tests 3 and 4: rate sensitivity, resource limits ----------
T3_KOFF = 100.0
T3_DT = 10.0


def t3_fixture(Ms, thermo=1.0, unlimited=True):
    ov = p31_fixture({"clamp_mstar": Ms, "koff_fixed": T3_KOFF, "t_end": T3_DT, "thermo": thermo,
                      "init_templates": {"n": [16, 16], "L": [12, 12], "seq": [0, (1 << 12) - 1]}})
    if unlimited:
        ov["copy_unlimited"] = True
    return ov


def t3_expected(c):
    rC, rG = gc_rates(c)
    out = {}
    for name, r in (("slow_tmpl0", rC), ("fast_tmpl1", rG)):
        q = r / (r + T3_KOFF)
        lam = r * q
        E = 1 + sum(q ** m for m in range(1, 11))
        cyc = 1 / lam + E / r + 1 / T3_KOFF
        out[name] = {"r": r, "q": q, "lam_dt": lam * T3_DT, "t_copy": E / r, "cycle": cyc, "mu": T3_DT / cyc,
                     "cap": math.ceil(T3_DT / (E / r + 1 / T3_KOFF))}
    out["ratio_fast_slow"] = out["slow_tmpl0"]["cycle"] / out["fast_tmpl1"]["cycle"]
    return out


def counts_by_template(o):
    n0 = n1 = 0
    for rec in o["copy_log"]:
        ts = rec["tseq"]
        n0 += int((ts == 0).sum())
        n1 += int((ts == (1 << 12) - 1).sum())
    return n0, n1


def t3_Ms():
    # c with r(G required) = 2 r(C required): 39 c / (0.39 + c) * 2 = 24 c / (0.006 + c)
    km = chem.PAIRS["GC"]["km_kmax"]
    a, b = km[1], km[0]
    # 2 a1 (b0 + c) = b1 (a0 + c)  ->  c = (b1 a0 - 2 a1 b0) / (2 a1 - b1)
    c = (b[1] * a[0] - 2 * a[1] * b[0]) / (2 * a[1] - b[1])
    return round(2000 * c)


def test_3():
    Ms = t3_Ms()
    c = Ms / 2000
    exp_ = t3_expected(c)
    fid = (16, 2048, 1)
    res = {"Ms": Ms, "c_mM": c, "r_ratio": exp_["fast_tmpl1"]["r"] / exp_["slow_tmpl0"]["r"], "expected": exp_}
    on = one_step("GC", chem.DIM_P31, 64, fid, t3_fixture(Ms), 7101, True, copy_log=True)
    n0, n1 = counts_by_template(on)
    res["multi"] = {"copies_slow": n0, "copies_fast": n1, "ratio": n1 / n0,
                    "rel_dev": n1 / n0 / exp_["ratio_fast_slow"] - 1,
                    "per_template_slow": n0 / (64 * 16 * 16), "per_template_fast": n1 / (64 * 16 * 16),
                    "occ_capped_share": float(on["ms_occupancy_capped_share"].mean()),
                    "limited_any_share": float(on["ms_limited_any_share"].mean())}
    off = one_step("GC", chem.DIM_P31, 64, fid, t3_fixture(Ms), 7102, False, copy_log=True)
    m0, m1 = counts_by_template(off)
    res["old"] = {"copies_slow": m0, "copies_fast": m1, "ratio": m1 / m0}
    res["pass"] = abs(res["multi"]["rel_dev"]) <= 0.05 and abs(res["old"]["ratio"] - 1) <= 0.02
    chem.MULTI_START = True
    print("t3", res, flush=True)
    return res


def test_4():
    Ms = t3_Ms()
    c = Ms / 2000
    exp_ = t3_expected(c)
    res = {"expected_ratio": exp_["ratio_fast_slow"]}
    nt = 32
    # (a) blocks low: thermo x10 and 1/10 the blocks (same c)
    Ms_a = round(Ms / 10)
    fid = (16, 2048, 1)
    o = one_step("GC", chem.DIM_P31, 256, fid, t3_fixture(Ms_a, thermo=10.0, unlimited=False), 7201, True,
                 copy_log=True, state=True)
    L = o["state"]["L"]
    used = (L[:, nt:].clamp(min=1) - 1).sum(1)            # letters added by the copies (unprimed: L - 1) per site
    n0, n1 = counts_by_template(o)
    res["a_blocks"] = {"supply_per_site": Ms_a, "max_used_per_site": int(used.max()), "mean_used": float(used.double().mean()),
                       "copies_slow": n0, "copies_fast": n1, "ratio": n1 / n0, "rel_dev": n1 / n0 / exp_["ratio_fast_slow"] - 1,
                       "limited_blocks_share": float(o["ms_limited_blocks_share"].mean())}
    res["a_blocks"]["pass"] = res["a_blocks"]["max_used_per_site"] <= Ms_a and abs(res["a_blocks"]["rel_dev"]) <= 0.05
    # (b) slots low: 64 slots per site (32 templates + 32 free)
    fid = (16, 64, 1)
    o = one_step("GC", chem.DIM_P31, 256, fid, t3_fixture(Ms), 7202, True, copy_log=True, state=True)
    L = o["state"]["L"]
    ncp = (L[:, nt:] > 0).sum(1)
    n0, n1 = counts_by_template(o)
    res["b_slots"] = {"free_slots_per_site": 64 - nt, "max_copies_per_site": int(ncp.max()), "mean_copies": float(ncp.double().mean()),
                      "copies_slow": n0, "copies_fast": n1, "ratio": n1 / n0, "rel_dev": n1 / n0 / exp_["ratio_fast_slow"] - 1,
                      "limited_slots_share": float(o["ms_limited_slots_share"].mean()),
                      "copies_logged": n0 + n1, "copies_in_state": int(ncp.sum())}
    res["b_slots"]["pass"] = (res["b_slots"]["max_copies_per_site"] <= 64 - nt and n0 + n1 == int(ncp.sum())
                              and abs(res["b_slots"]["rel_dev"]) <= 0.05)
    res["pass"] = res["a_blocks"]["pass"] and res["b_slots"]["pass"]
    print("t4", res, flush=True)
    return res


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "validate_p31b.json")
    tests = sys.argv[2] if len(sys.argv) > 2 else "1234"
    res = ckpt.load_json(out) or {}
    res["chem_p31b_md5"] = md5("chem_p31b.py")
    for t in tests:
        t0 = time.time()
        r = globals()[f"test_{t}"]()
        r["seconds"] = round(time.time() - t0, 1)
        res[f"test_{t}"] = r
        ckpt.atomic_json(res, out)
        print(f"TEST {t}", "PASS" if r.get("pass") else "FAIL", flush=True)
    res["thermal"] = G().summary()
    ckpt.atomic_json(res, out)


if __name__ == "__main__":
    main()
