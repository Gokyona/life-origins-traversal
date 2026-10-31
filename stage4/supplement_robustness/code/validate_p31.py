"""Part III step P3-1 validation tests 1-9 (spec docs/part3_p31_spec.md, section 5) for chem_p31.py.

  (1) duplex thermodynamics: chem_p31's Delta G(T) for the Fang et al. 2025 9-pair duplexes
      5'-CUGA X GUAG-3' (X:Y = A:U, G:C centre) at 25 C, from its reference function nn_seq_dG, against an
      independent NNDB calculation written here from the alphabet_tm.py table values (pass: within 0.3
      kcal/mol; the measured Table 1 values are reported beside them). The GPU path (_nn_hs, used in
      simulate) is checked against nn_seq_dG on all 2-letter 8-mers of each pair (pass: max |diff| < 0.01).
      AS2U: the offset gives AS2U - AU = -delta per pair at 25 C (-1.1 at delta = 1.1), by construction;
  (2) extension kinetics: simulate with fixed homopolymer templates (40 letters), 3-mer primers, M*
      clamped (so the per-letter block concentration c is fixed), k_off fixed to the expected r (q ~ 0.5),
      no hydrolysis, growth or initiation, 296 K and 100 mM Mg (f_T = f_Mg = 1); c = 0.1, 0.3, 1, 3, 10 K_m
      for both letters of every pair. From the logged copies, q = added / (added + stops) and r = k_off q /
      (1 - q); pass: |r - r_MM| <= 3 standard errors in every case (r_MM = k_max c / (K_m + c) at the c
      actually set, after rounding of the clamped count);
  (3) pulse integration: closed form (chem_p31.pulse_on_time + integrated hazard, as in simulate) against
      an explicit stochastic simulation of 1e6 duplexes and 1e6 bonds with every pulse resolved in 20
      sub-steps (and the cool part of each period in 20 sub-steps), tau_p = 1 s, 1 min, 1 h; pass: survival
      within 2% (relative) in all six cases;
  (4) Stage 4 regression: chem_p31 with PAIR = None against stored chem4c3 results on 16 runs (box2 points of
      the 4c-3 scan incl. #1558 and #3709; main run, L1 and L2 paired controls): bit-identical;
  (5) activation calibration (cold-gap change, 30 Oct 2026): simulator reagent network = calibration ODE on
      Aitken Fig. 5 (<= 0.15 mM), and the k_act range brackets the observed 2-h A*A (see test_5);
  (6) cold-gap acceptance (replaces the old cold-benefit check): see test_6 / cold_accept.py;
  (7) error threshold: explicit positive control (see test_7), eps * L_copy = 0.15 ... 3, criteria of
      Stage 4 test 9b, written to the output before any run;
  (8) null false positives: the test 7 condition with random letters (eps_0 = eps_1 = 0.5), 64 reps, both
      copy_random controls for every candidate; pass: L2_strict <= 2 of 64;
  (9) kill and resume: kill_resume_p31.sh with run_p31.py (bit-identical results).
    python validate_p31.py [out.json] [tests, e.g. 123456789]
"""
import json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_p31 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402

dev = torch.device("xpu")
guard = None
FID = (16, 256, 600)
REPS = 16
# a neutral explicit condition (P3-1 dimensions; overrides in each test)
BASE = dict(T=22.85, pH=7.5, Mg=0.1, salinity=0.25, s_exp=0.0, P=100.0, feed=10.0, tau=1e12, k_c=1e-12, p0=0.0,
            depth=0.0, thermo=1.0, mineral_type=0.0, x_nuc=1.0, k_UV=0.0, avail_time=1e8, L_t=12.0, act_a=100.0,
            T_p=22.85, tau_p=1 / 60, f_p=1e-3, delta=1.1, k_on=3.16e6, m_Mg=1.0, E_ext=50.0, eps_0=1e-3, eps_1=1e-3,
            t_half_act=3.0, E_act_h=75.0, k_h_fixed=0.0, init_factor=0.0)
BASE["E_ext"] = 18.4                 # E_ext band 15-20 kJ/mol (f_T = 1 at 296 K in tests 2 and 3)
BASE["act_a"] = 0.0                  # no reagent in the clamped-M* fixtures (tests 2, 7, 8): no side products, I = 0
BASE["K_I"] = 4.0


def G():
    global guard
    if guard is None:
        guard = thermal.Guard(job="P3-1 validation")
    return guard


def gen(seed):
    g = torch.Generator(device=dev)
    g.manual_seed(seed)
    return g


# ---------- test 1 ----------
# independent NNDB numbers (Turner 2004, as tabulated in part3/alphabet_tm.py)
_DG = {'AA': -0.9, 'UU': -0.9, 'AU': -1.1, 'UA': -1.3, 'CU': -2.1, 'AG': -2.1, 'CA': -2.1, 'UG': -2.1, 'GU': -2.2,
       'AC': -2.2, 'GA': -2.4, 'UC': -2.4, 'CG': -2.4, 'GG': -3.3, 'CC': -3.3, 'GC': -3.4}
_DH = {'AA': -6.8, 'UU': -6.8, 'AU': -9.4, 'UA': -7.7, 'CU': -10.5, 'AG': -10.5, 'CA': -10.4, 'UG': -10.4, 'GU': -11.4,
       'AC': -11.4, 'GA': -12.4, 'UC': -12.4, 'CG': -10.6, 'GG': -13.4, 'CC': -13.4, 'GC': -14.9}


def nndb_dG(s, T_C):
    g, h = 4.09, 3.61
    for i in range(len(s) - 1):
        g += _DG[s[i:i + 2]]
        h += _DH[s[i:i + 2]]
    for e in (s[0], s[-1]):
        if e in "AU":
            g += 0.45
            h += 3.72
    S = (h - g) / 310.15
    return h - (T_C + 273.15) * S


def test_1():
    out = {"fang": {}}
    fang = {"A:U": ("CUGAAGUAG", -13.6, 46.8), "G:C": ("CUGAGGUAG", -15.2, 55.6)}
    ok = True
    for k, (s, meas, tm) in fang.items():
        code = chem.nn_seq_dG(s, 25.0)
        ref = nndb_dG(s, 25.0)
        out["fang"][k] = {"seq": s, "code_dG25": round(code, 4), "nndb_dG25": round(ref, 4), "diff": round(code - ref, 4),
                          "measured_dG25_Fang_Table1": meas, "measured_Tm_C": tm,
                          "code_minus_measured": round(code - meas, 3)}
        ok = ok and abs(code - ref) <= 0.3
    # GPU path (_nn_hs) vs the string reference on all 2-letter 8-mers of each pair, several temperatures
    gpu = {}
    for pair in ("AS2U", "GC", "AU"):
        pt = chem._pair_tables(pair, "cpu")
        nm = chem.PAIRS[pair]["letters"]
        delta = 1.1 if chem.PAIRS[pair]["offset"] else 0.0
        cs = torch.arange(256, dtype=torch.int64)
        l = torch.full((256,), 8, dtype=torch.int64)
        H, S = chem._nn_hs(cs, l, pt, torch.full((256,), delta))
        md = 0.0
        for T in (0.0, 25.0, 60.0, 95.0):
            g_code = (H - (T + 273.15) * S).double()
            for i in range(256):
                st = "".join(nm[(i >> j) & 1] for j in range(8))
                md = max(md, abs(float(g_code[i]) - chem.nn_seq_dG(st, T, delta=delta)))
        gpu[pair] = round(md, 6)
        ok = ok and md < 0.01
    out["gpu_vs_reference_max_abs_diff"] = gpu
    # AS2U offset: one 9-pair A/s2U duplex vs the same A/U duplex at 25 C, delta = 1.1 -> -1.1 per pair
    s = "AUUAAUAUU"
    d_as2u = chem.nn_seq_dG(s, 25.0, delta=1.1) - chem.nn_seq_dG(s, 25.0)
    out["AS2U_minus_AU_9mer_dG25"] = round(d_as2u, 4)
    out["AS2U_minus_AU_per_pair"] = round(d_as2u / 9, 4)
    out["Fang_A_s2U_minus_A_U_single_pair"] = -1.1
    ok = ok and abs(d_as2u / 9 + 1.1) < 1e-9
    out["pass"] = bool(ok)
    return out


# ---------- test 2 ----------

def test_2():
    out = {}
    allok = True
    S_, n_t, n_p, Lt = 4096, 8, 2000, 40
    for pair in ("AS2U", "GC", "AU"):
        chem.PAIR = pair
        tab = chem.PAIRS[pair]["km_kmax"]
        for X in (0, 1):                         # template letter X requires letter 1 - X
            req = 1 - X
            Km, kmax = tab[req]
            tseq = 0 if X == 0 else (1 << Lt) - 1
            pseq = 7 if req == 1 else 0          # 3-mer primer of the required letter
            cs = [0.1, 0.3, 1.0, 3.0, 10.0]
            B = len(cs)
            thermo, ms, c_set, r_exp = [], [], [], []
            for c in cs:
                target = c * Km
                th = min(1.0, 2 * target * chem.V / 1000.0)      # >= 1000 molecules in the clamped pool
                m = round(2 * target * chem.V / th)
                thermo.append(th); ms.append(float(m))
                c_act = 0.5 * th * m / chem.V
                c_set.append(c_act)
                r_exp.append(kmax * c_act / (Km + c_act))
            ov = dict(BASE, thermo=torch.tensor(thermo, dtype=torch.float64), L_t=float(Lt),
                      clamp_mstar=torch.tensor(ms, dtype=torch.float64),
                      koff_fixed=torch.tensor(r_exp, dtype=torch.float64), copy_unlimited=True, copy_log=True,
                      t_end=torch.tensor([400.0 / r for r in r_exp], dtype=torch.float64),
                      init_templates={"n": [n_t, n_p], "L": [Lt, 3], "seq": [tseq, pseq]})
            u = torch.full((B, chem.DIM_P31), 0.5, device=dev)
            o = chem.simulate(u, gen(200 + X), (16, S_, 300), override=ov)
            G().check()
            log = {k: torch.cat([d[k] for d in o["copy_log"]]) for k in o["copy_log"][0]}
            rows = []
            for b in range(B):
                m = log["run"] == b
                l, l0, mx = log["l"][m], log["l0"][m], log["maxl"][m]
                succ = int((l - l0).sum()) - int(m.sum())        # first letter is forced (start = >= 1 letter)
                stops = int((l < mx).sum())
                n_tr = succ + stops
                q = succ / max(n_tr, 1)
                kf = r_exp[b]
                r_hat = kf * q / (1 - q) if q < 1 else float("inf")
                se = kf / (1 - q) ** 2 * math.sqrt(q * (1 - q) / max(n_tr, 1))
                ok = abs(r_hat - r_exp[b]) <= 3 * se
                allok = allok and ok
                rows.append({"c_over_Km": cs[b], "c_mM": c_set[b], "r_MM": r_exp[b], "r_hat": r_hat, "se": se,
                             "z": (r_hat - r_exp[b]) / se if se > 0 else None, "copies": int(m.sum()),
                             "trials": n_tr, "pass": bool(ok)})
            out[f"{pair}_template{X}_needs{req}"] = {"K_m": Km, "k_max": kmax, "rows": rows}
    chem.PAIR = None
    out["pass"] = bool(allok)
    return out


# ---------- test 3 ----------

def test_3():
    """Duplex: the G/C 8-mer 01101001 (k_on 3.16e6 /M/s), base 20 C; bonds: k_nuc at pH 9,
    5 mM Mg, 0.25 M salt, base 20 C. For each tau_p (f_p = 0.01, period 100 tau_p) T_p is chosen so that the
    pulses of one step give survival ~0.5 (duplex: 10 pulses per step; bonds: n pulses per step as needed),
    then the step is integrated in closed form and by explicit sub-stepping (fixed seed)."""
    rng = np.random.default_rng(31)
    pt = chem._pair_tables("GC", "cpu")
    H, S = chem._nn_hs(torch.tensor([~0b01101001]), torch.tensor([8]), pt, torch.tensor([0.0]))
    H, S = float(H), float(S)
    kon = 3.16e6 * 3600
    koff = lambda T: kon * math.exp(min((H - (T + 273.15) * S) / (chem.R_KCAL * (T + 273.15)), 60.0))
    p = {"pH": torch.tensor([9.0]), "Mg": torch.tensor([0.005]), "salinity": torch.tensor([0.25]), "s_exp": torch.tensor([0.0])}
    kh = lambda T: float(chem.k_nuc(torch.tensor([float(T)], dtype=torch.float64), p))
    Tb = 20.0
    N = 1_000_000
    out, ok = {}, True
    for name, tau in (("1s", 1 / 3600), ("1min", 1 / 60), ("1h", 1.0)):
        f_p = 0.01
        P = tau / f_p
        for kind in ("duplex", "bond"):
            kf = koff if kind == "duplex" else kh
            if kind == "duplex":
                n = 10
            else:
                n = max(1, int(round(0.69 / (kh(100.0) * tau))))     # pulses per step at <= 100 C
            # T_p such that n tau k(T_p) ~ ln 2 (bisection on the code's rate)
            lo, hi = Tb, 150.0
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if n * tau * kf(mid) < math.log(2):
                    lo = mid
                else:
                    hi = mid
            Tp = 0.5 * (lo + hi)
            dt = n * P
            h = float(chem.pulse_on_time(torch.tensor([0.0], dtype=torch.float64), torch.tensor([dt], dtype=torch.float64),
                                         torch.tensor([P], dtype=torch.float64), torch.tensor([tau], dtype=torch.float64)))
            closed = math.exp(-(kf(Tb) * (dt - h) + kf(Tp) * h))
            # explicit: each period = cool part (20 sub-steps) + pulse (20 sub-steps); binomial losses
            alive = N
            pc = 1 - math.exp(-kf(Tb) * (P - tau) / 20)
            ph = 1 - math.exp(-kf(Tp) * tau / 20)
            for _ in range(n):
                for _ in range(20):
                    alive -= rng.binomial(alive, pc)
                for _ in range(20):
                    alive -= rng.binomial(alive, ph)
            expl = alive / N
            rel = abs(closed - expl) / expl
            ok = ok and rel <= 0.02
            out[f"{kind}_{name}"] = {"tau_p_h": tau, "f_p": f_p, "pulses": n, "T_base": Tb, "T_p": round(Tp, 3),
                                     "pulse_on_time_h": h, "survival_closed": closed, "survival_explicit": expl,
                                     "rel_diff": rel}
    out["pass"] = bool(ok)
    return out


# ---------- test 4 ----------

def _box2_points(n1=4096, seed=2026):
    """The 4c-3 box2 scan point set (lot-s4/make_ref4c3.box2_points, re-implemented with chem4c3)."""
    import chem4c3
    from scipy.stats import qmc
    Xb = qmc.Sobol(d=chem4c3.DIM, scramble=True, seed=seed).random(1 << 16).astype(np.float32)
    ix = {n: chem4c3.NAMES.index(n) for n in ("T", "T_m8", "depth", "k_t", "thermo", "x_nuc")}
    Tm = 10 + 80 * Xb[:, ix["T_m8"]]
    T = Tm - 15 * Xb[:, ix["T"]]
    Xb[:, ix["T"]] = np.clip(T / 150.0, 0.0, 1.0)
    dmin = np.clip((Tm - T) / 10.0 + 0.02, 0.0, 1.0)
    Xb[:, ix["depth"]] = np.maximum(0.5 + 0.5 * Xb[:, ix["depth"]], dmin)
    Xb[:, ix["k_t"]] = 0.75 + 0.25 * Xb[:, ix["k_t"]]
    Xb[:, ix["thermo"]] = 0.375 + 0.625 * Xb[:, ix["thermo"]]
    Xb[:, ix["x_nuc"]] = 0.5 + 0.5 * Xb[:, ix["x_nuc"]]
    p = chem4c3.decode(torch.tensor(Xb))
    tab, _, theta, _ = chem4c3.phase_table(p)
    life = 1.0 / chem4c3.turnover_rate(chem4c3.L_H, tab["knuc"][:, 2], p["k_UV"] * p["UV"], 1 / p["tau"], theta, p["p_min"])
    sel = np.nonzero(((life / p["period"]) >= 1.0).numpy())[0][:n1]
    return Xb[sel]


def _flat(o):
    """Comparable view of a simulate() result: tensors as CPU tensors, lists / dicts as repr."""
    d = {}
    for k, v in o.items():
        if isinstance(v, torch.Tensor):
            d[k] = v.detach().cpu()
        elif v is None:
            d[k] = None
        else:
            d[k] = repr(v)
    return d


def _same(a, b):
    bad = []
    for k in sorted(set(a) | set(b)):
        x, y = a.get(k, "<missing>"), b.get(k, "<missing>")
        if isinstance(x, torch.Tensor) and isinstance(y, torch.Tensor):
            if x.shape != y.shape or x.dtype != y.dtype:
                bad.append(k)
            elif x.dtype.is_floating_point:
                if not torch.equal(torch.nan_to_num(x, nan=-12345.0), torch.nan_to_num(y, nan=-12345.0)) or \
                        not torch.equal(torch.isnan(x), torch.isnan(y)):
                    bad.append(k)
            elif not torch.equal(x, y):
                bad.append(k)
        elif x != y:
            bad.append(k)
    return bad


def test_4():
    import pickle
    import chem4c3
    store = os.path.join(HERE, "regress_4c3_stored.pkl")
    idx = [1558, 3709] + list(range(14))
    X = _box2_points()[idx]
    u = torch.tensor(X, dtype=torch.float32, device=dev)
    fid = FID
    seed = 4401

    def run(mod):
        o = mod.simulate(u, gen(seed), fid, record=True)
        G().check()
        s1 = mod.paired_control(u, seed, fid, o)
        s2 = mod.paired_control_l2(u, seed, fid, o, both=True)
        G().check()
        d = _flat(o)
        for k, v in zip(("L1_strict", "ab_cat", "ab_ctl"), s1):
            d["pc_" + k] = v.cpu()
        for k, v in zip(("L2_strict", "S_main", "S_ctl"), s2):
            d["pc_" + k] = v.cpu()
        return d
    t0 = time.time()
    if os.path.exists(store):
        ref = pickle.load(open(store, "rb"))
    else:
        ref = run(chem4c3)
        tmp = store + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(ref, f)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, store)
    t_ref = time.time() - t0
    chem.PAIR = None
    new = run(chem)
    bad = _same(ref, new)
    md5 = lambda f: os.popen(f"md5sum {os.path.join(HERE, f)}").read().split()[0]
    return {"points_box2": idx, "fidelity": list(fid), "seed": seed, "fields_compared": len(ref),
            "differing_fields": bad, "chem4c3_md5": md5("chem4c3.py"), "stored": os.path.basename(store),
            "copies_total": float(ref["copies"].sum()), "L2_loose": int(ref["L2_loose"].sum()),
            "L2_strict": int(ref["pc_L2_strict"].sum()), "L1_loose": int(ref["L1_loose"].sum()),
            "seconds_ref": t_ref, "pass": len(bad) == 0}


# ---------- test 5 (cold-gap version) ----------

def test_5():
    """Activation calibration (cold-gap change (a)): (i) the simulator's reagent + monomer network (CPU, Fig. 5
    protocol: 20 mM monomer, 200 mM reagent once, 22 C and -3 C, 400 min) equals the calibration ODE (reagent_fit.py,
    fit B on Aitken et al. 2023 Fig. 5) within 0.15 mM at all 33 record points; (ii) at 22 C, the A*A after 2 h at the
    low and high ends of the k_act range brackets the observed 2.0-3.3 mM (2.0 mM at 2 h with 1MeI; 2.7-3.3 mM peak
    with the other catalysts)."""
    import subprocess
    r = subprocess.run([sys.executable, os.path.join(HERE, "reagent_simcheck.py"), "chem_p31"], capture_output=True,
                       text=True, cwd=HERE)
    chk = json.load(open(os.path.join(HERE, "reagent_simcheck.json")))
    ok1 = all(v["max_abs_diff_mM"] <= 0.15 and v["liquid"] for v in chk.values())
    import reagent_fit as RF
    FITB = json.load(open(os.path.join(HERE, "reagent_fit.json")))["B_with_Ea"]["RT"]
    ends = []
    for f in (0.1, 1.0, 10.0):
        x = RF.sim((FITB["k_act"] * f, FITB["k_r"], FITB["y"], FITB["k_h"]), [2.0])[2][0]
        ends.append(float(x))
    ok2 = ends[0] <= 2.0 and ends[2] >= 3.3
    return {"sim_vs_ode": {k: {kk: v[kk] for kk in ("max_abs_diff_mM", "liquid", "reagent_end_mM")} for k, v in chk.items()},
            "AA_2h_RT_at_k_act_low_fit_high_mM": ends, "observed_2h_mM": [2.0, 3.3],
            "pass_sim_vs_ode": bool(ok1), "pass_bracket": bool(ok2), "pass": bool(ok1 and ok2)}


# ---------- test 6 (cold-gap acceptance) ----------

def test_6():
    """Cold-gap acceptance test (replaces the old test 6): cold_accept.py. One parameter set (reagent-pool fit B,
    E_ext scanned over 0-20 kJ/mol) against (A) >= +1 at 24 h 65/87%, (B) >= +7 at 5 d 7/60%, (C) pre-activated
    >= +7 at 24 h 42/8% (+3 C); tolerance +-5 points. Fails honestly: the residual rate factor is reported."""
    import subprocess
    r = subprocess.run([sys.executable, os.path.join(HERE, "cold_accept.py")], capture_output=True, text=True, cwd=HERE)
    acc = json.load(open(os.path.join(HERE, "cold_accept.json")))
    return {"best": acc["best"], "at_best": acc["at_best"], "accepted_sets": acc["accepted_sets"],
            "A_residual_unexplained": acc["A_residual_unexplained"], "phi": acc["phi"], "scan": acc["scan"],
            "no_inhibition_scan": acc["no_inhibition_scan"], "stdout": r.stdout[-3000:], "pass": bool(acc["pass"])}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "validate_p31.json")
    which = sys.argv[2] if len(sys.argv) > 2 else "123456"
    res = ckpt.load_json(path) or {}
    res["fidelity"] = FID
    res["chem_p31_md5"] = os.popen(f"md5sum {os.path.join(HERE, 'chem_p31.py')}").read().split()[0]
    tests = {"1": test_1, "2": test_2, "3": test_3, "4": test_4, "5": test_5, "6": test_6}
    try:
        import validate_p31_sim as VS          # tests 7 and 8 (simulation positive control)
        tests.update({"7": VS.test_7, "8": VS.test_8})
    except ImportError:
        pass
    for k, fn in tests.items():
        if k not in which:
            continue
        t = time.time()
        if k == "7":                                   # criterion on disk before any run
            res["test_7_criterion_preregistered"] = VS.CRITERION
            ckpt.atomic_json(res, path)
        res["test_" + k] = fn()
        res["test_" + k]["seconds"] = time.time() - t
        print("test", k, "pass" if res["test_" + k]["pass"] else "FAIL", f"({time.time() - t:.0f}s)", flush=True)
        ckpt.atomic_json(res, path)
    if guard is not None:
        res["thermal"] = guard.summary()
    ckpt.atomic_json(res, path)


if __name__ == "__main__":
    main()
