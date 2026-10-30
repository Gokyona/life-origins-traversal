"""P3-1 test 6 follow-up: the cold-benefit gap (Aitken, Wright, Radakovic & Szostak 2023, JACS 145:16142).
CPU only. Matched protocols and metrics:
  (A) Table 2, standard vs entry 14: CMP + GMP 10 mM each, 2AI 10 mM, catalysts 50 mM, MBA/MeNC 200 mM once,
      Mg 30 mM, pH 8, primer 2 uM on template overhang 5'-CCG-3' (first added letter C on G); yield >= +1 at 24 h:
      65% at room temperature (RT) vs 87% at -3 C (liquid, no eutectic);
  (B) Fig. 7C / S4: NMPs 10 mM each, helper trimers 500 uM each, 2AI 30 mM, catalysts 60 mM, MBA/MeNC 200 mM added
      every 24 h, template 5'-GCGCCUCAU overhang; yield >= +7 after 5 days: 60% at -3 C vs 7% at RT;
  (C) Fig. 7A/B pre-activated control: *N 10 mM, *NNN 500 uM, no reagent, 24 h: >= +7 42% at RT vs 8% at 3 C.
Activation kinetics from Fig. 5 (31P NMR, AMP 20 mM, 2AI 10 mM, catalysts 50 mM, MBA/MeNC 200 mM once; RT and
-3 C), digitised by eye from a 300-dpi render (+-0.2 mM): AMP and A*A.
Minimal activation model per temperature: reagent R (relative, +1 per addition, loss k_r), AMP M consumed at
k_m R M, A*A formed at y k_m R M / 2 (two nucleotides per A*A; y = share of consumed AMP that ends as A*A rather
than *A, AppA, impurities), A*A lost at k_l (hydrolysis + exchange with 2AI). Copying: homogeneous per-letter
Poisson rate r(t) = s k_max c / (K_m + c) f_T, c = A*A-equivalent of the required letter (scaled with that
nucleotide's concentration), the scale s fitted once at RT per protocol, then -3 C predicted. Room temperature
taken as 22 C.
    python cold_gap.py [out.json]
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
from scipy.stats import poisson

R = 8.314
T_RT, T_C = 22.0, -3.0
# Fig. 5 digitised (t min, AMP mM, A*A mM)
FIG5 = {
    T_RT: [(0, 20.0, 0.0), (10, 14.5, 0.6), (20, 12.5, 1.8), (30, 11.6, 2.3), (45, 10.6, 2.6), (60, 9.8, 2.8),
           (80, 9.1, 2.9), (100, 8.8, 2.9), (120, 8.6, 2.8), (150, 8.3, 2.6), (200, 7.9, 2.3), (250, 7.7, 2.0),
           (300, 7.6, 1.8), (350, 7.5, 1.6), (400, 7.5, 1.5)],
    T_C: [(0, 20.0, 0.0), (10, 15.5, 0.1), (30, 14.6, 0.4), (60, 13.8, 1.1), (80, 13.4, 1.5), (100, 13.0, 1.9),
          (130, 12.6, 2.3), (160, 12.3, 2.6), (200, 11.6, 2.9), (250, 11.1, 3.2), (300, 10.6, 3.35), (350, 10.3, 3.4),
          (400, 10.0, 3.4)],
}


def act_rhs(k, M0):
    k_m, k_r, y, k_l = k

    def f(t, z):
        Rg, M, x = z
        v = k_m * Rg * M
        return [-k_r * Rg, -v + 2 * k_l * x, y * v / 2 - k_l * x]
    return f


def act_sim(k, M0, t_h, adds=(0.0,)):
    """A*A x(t) (mM) and AMP for additions of reagent at times adds (h); A*A loss returns AMP."""
    ts = np.asarray(t_h, float)
    z = [0.0, M0, 0.0]
    out_x, out_M = np.zeros_like(ts), np.zeros_like(ts)
    edges = list(adds) + [max(ts.max(), adds[-1]) + 1e-9]
    for i in range(len(adds)):
        z[0] += 1.0
        a, b = edges[i], edges[i + 1]
        m = (ts >= a) & (ts <= b) if i == len(adds) - 1 else (ts >= a) & (ts < b)
        sol = solve_ivp(act_rhs(k, M0), (a, b), z, t_eval=ts[m] if m.any() else None, rtol=1e-8, atol=1e-10,
                        dense_output=True)
        if m.any():
            out_M[m], out_x[m] = sol.y[1], sol.y[2]
        z = list(sol.sol(b))
    return out_x, out_M


def fit_fig5(T):
    d = np.array(FIG5[T])
    t, M, x = d[:, 0] / 60, d[:, 1], d[:, 2]

    def res(lk):
        k = np.exp(lk)
        k[2] = min(k[2], 1.0)
        xs, Ms = act_sim(k, 20.0, t)
        return np.concatenate([xs - x, 0.5 * (Ms - M)])
    best = None
    for g in ([1, 1, 0.5, 0.1], [0.3, 0.1, 0.5, 0.01], [3, 3, 0.5, 0.3], [0.1, 0.02, 0.5, 0.01]):
        r = least_squares(res, np.log(g), bounds=(np.log([1e-4, 1e-4, 1e-3, 1e-5]), np.log([100, 100, 1.0, 10])))
        if best is None or r.cost < best.cost:
            best = r
    k = np.exp(best.x)
    return {"k_m_per_h": k[0], "k_r_per_h": k[1], "y": k[2], "k_l_per_h": k[3],
            "rms_mM": float(np.sqrt(np.mean(best.fun ** 2)))}, k


def arr_ratio(E, T1=T_C, T0=T_RT):
    return math.exp(-E * 1000 / R * (1 / (T1 + 273.15) - 1 / (T0 + 273.15)))


def E_from_ratio(ratio, T1=T_C, T0=T_RT):
    return -R * math.log(ratio) / (1 / (T1 + 273.15) - 1 / (T0 + 273.15)) / 1000


def copy_lambda(k, fT, s, Km, kmax, M0_letter, t_end, adds, Km_T=None):
    ts = np.linspace(0, t_end, 6001)
    x, _ = act_sim(k, 2 * M0_letter, ts, adds)          # fit made at 20 mM AMP; x scales with M0
    c = x / 2                                           # A*A-equivalent of this letter (letter's own M0)
    Kd = Km if Km_T is None else Km_T
    r = s * kmax * c / (Kd + c) * fT
    return float(np.sum(0.5 * (r[1:] + r[:-1]) * np.diff(ts)))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "cold_gap.json")
    out = {"fig5_fit": {}}
    kk = {}
    for T in (T_RT, T_C):
        out["fig5_fit"][str(T)], kk[T] = fit_fig5(T)
    fr, fc = out["fig5_fit"][str(T_RT)], out["fig5_fit"][str(T_C)]
    out["apparent_E_kJ"] = {n: E_from_ratio(fc[n] / fr[n]) for n in ("k_m_per_h", "k_r_per_h", "k_l_per_h")}
    out["apparent_E_kJ"]["note"] = "E from the -3 C / 22 C ratio of each fitted rate; y (selectivity) is a ratio"
    out["y_ratio_cold_over_RT"] = fc["y"] / fr["y"]
    # (C) pre-activated control: constant c (slow monomer hydrolysis, 4e-3 /h, Kaddour/Walton kinetic model): Poisson
    # rates giving 42% (RT) and 8% (3 C) >= +7 at 24 h -> apparent extension E (3 C, not -3 C)
    lam = lambda p: float(next(l for l in np.linspace(0.01, 30, 300000) if poisson.sf(6, l) >= p))
    l42, l8 = lam(0.42), lam(0.08)
    E_ext_app = E_from_ratio(l8 / l42, T1=3.0)
    out["preactivated_control"] = {"lambda_RT": l42, "lambda_3C": l8, "per_step_ratio_3C_over_RT": l8 / l42,
                                   "E_ext_apparent_kJ": E_ext_app}
    # matched copying, protocols A and B; letter-specific K_m, k_max for the first added letter of A (C on G,
    # Fang 2025 C*C on G: 0.39 mM, 39 /h); B uses the same values (4-letter template: homogeneous step rate)
    Km, kmax = 0.39, 39.0
    protos = {"A_plus1_24h": dict(M0=10.0, t=24.0, adds=(0.0,), n=1, p_RT=0.65, p_C=0.87),
              "B_plus7_5d": dict(M0=10.0, t=120.0, adds=(0.0, 24.0, 48.0, 72.0, 96.0), n=7, p_RT=0.07, p_C=0.60)}
    yld = lambda lam_, n: float(poisson.sf(n - 1, lam_))

    def scenario(name, k_cold, fT, Km_cold=None):
        res = {}
        for pn, P in protos.items():
            unit = copy_lambda(kk[T_RT], 1.0, 1.0, Km, kmax, P["M0"], P["t"], P["adds"])
            # s: RT yield = observed
            s = None
            lo, hi = 1e-9, 1e3
            for _ in range(200):
                mid = math.sqrt(lo * hi)
                if yld(mid * unit, P["n"]) < P["p_RT"]:
                    lo = mid
                else:
                    hi = mid
            s = math.sqrt(lo * hi)
            lc = copy_lambda(k_cold, fT, s, Km, kmax, P["M0"], P["t"], P["adds"], Km_T=Km_cold)
            res[pn] = {"s": s, "yield_RT": P["p_RT"], "yield_cold_model": yld(lc, P["n"]), "yield_cold_obs": P["p_C"],
                       "factor_model": yld(lc, P["n"]) / P["p_RT"], "factor_obs": P["p_C"] / P["p_RT"]}
        return res

    kC, kR = kk[T_C], kk[T_RT]
    mix = lambda use_cold: np.array([kC[i] if use_cold[i] else kR[i] for i in range(4)])
    sc = {}
    # current chem_p31 (test 6): constant reagent (k_r = 0), activation T-independent, activated-block loss with
    # E_act_h, extension E_ext; here on the fitted RT activation with k_r -> 0
    for Eh in (50, 75, 100):
        for Ee in (20, 50, 80):
            k0 = kR.copy(); k0[1] = 1e-6
            k1 = k0.copy(); k1[3] = kR[3] * arr_ratio(Eh)
            r_ = {}
            for pn, P in protos.items():
                unit = copy_lambda(k0, 1.0, 1.0, Km, kmax, P["M0"], P["t"], P["adds"])
                lo, hi = 1e-9, 1e3
                for _ in range(200):
                    mid = math.sqrt(lo * hi)
                    lo, hi = (mid, hi) if yld(mid * unit, P["n"]) < P["p_RT"] else (lo, mid)
                lc = copy_lambda(k1, arr_ratio(Ee), math.sqrt(lo * hi), Km, kmax, P["M0"], P["t"], P["adds"])
                r_[pn] = yld(lc, P["n"]) / P["p_RT"]
            sc[f"current_model_Eh{Eh}_Eext{Ee}"] = r_
    out["current_model_factors"] = sc
    # decomposition with the fitted activation (RT -> -3 C one term at a time), E_ext from the pre-activated control
    fT_c = arr_ratio(E_ext_app)
    terms = ["k_m (activation rate)", "k_r (reagent loss)", "y (selectivity)", "k_l (A*A loss)"]
    dec = {"all_RT_terms": scenario("rt", kR, 1.0), "all_cold_terms_Eext_app": scenario("all", kC, fT_c),
           "all_cold_terms_Eext0": scenario("all0", kC, 1.0)}
    for i, n in enumerate(terms):
        u = [False] * 4; u[i] = True
        dec["only_" + n] = scenario(n, mix(u), 1.0)
        v = [True] * 4; v[i] = False
        dec["all_but_" + n] = scenario(n, mix(v), fT_c)
    dec["only_extension_Eext_app"] = scenario("ext", kR, fT_c)
    # K_m(T) with binding enthalpy dH_b (tighter in the cold), E_ext 0 and the apparent value
    for dHb in (-20, -50, -100):
        Kc = Km * math.exp(dHb * 1000 / R * (1 / (T_C + 273.15) - 1 / (T_RT + 273.15)))
        dec[f"all_cold_terms_Eext_app_Km_dH{dHb}"] = scenario("km", kC, fT_c, Km_cold=Kc)
    out["decomposition"] = dec
    out["per_step_obs_A"] = 0.87 / 0.65
    out["per_step_obs_A_pow7"] = (0.87 / 0.65) ** 7
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"fig5_fit": out["fig5_fit"], "apparent_E_kJ": out["apparent_E_kJ"],
                      "y_ratio": out["y_ratio_cold_over_RT"], "preactivated": out["preactivated_control"]}, indent=1))
    for k, v in dec.items():
        print(f"{k:55s} A: {v['A_plus1_24h']['yield_cold_model']:.3f} (obs 0.87)  B: {v['B_plus7_5d']['yield_cold_model']:.3f} (obs 0.60)")
    print("current model factors (A, B):", {k: (round(v['A_plus1_24h'], 3), round(v['B_plus7_5d'], 3)) for k, v in sc.items()})


if __name__ == "__main__":
    main()


# ---------- empirical version (robust to the poorly identified ODE fit) ----------
# A*A after one reagent addition at 20 mM AMP: Fig. 5 points (0-400 min), then exponential decay with k_tail:
# RT: from the Fig. 5 tail (2.9 mM at 80 min -> 1.5 mM at 400 min): ln(2.9/1.5)/5.33 h = 0.124 /h;
# -3 C: Fig. 5 flat to 400 min; Fig. 6 (-3 C, 40 mM AMP, 30 mM 2AI) first cycle 4.0 mM at ~3 h -> 1.2 mM at 24 h
# with activation still running: k <= ln(4/1.2)/21 = 0.057 /h (upper bound); 0.02 /h used as a lower case.
K_TAIL = {T_RT: math.log(2.9 / 1.5) / (320 / 60), T_C: math.log(4.0 / 1.2) / 21.0}


def aa_single(T, t_h, k_tail=None):
    d = np.array(FIG5[T])
    tt, xx = d[:, 0] / 60, d[:, 2]
    k = K_TAIL[T] if k_tail is None else k_tail
    t_h = np.asarray(t_h, float)
    return np.where(t_h <= tt[-1], np.interp(t_h, tt, xx), xx[-1] * np.exp(-k * (t_h - tt[-1])))


def aa_curve(T, t_h, adds, k_tail=None, head_from=None):
    """A*A (mM, per 20 mM nucleotide) for reagent additions at times adds: superposed single-addition curves
    (Fig. 6 shows similar cycles, peaks 4-6 mM). head_from: take the 0-400 min part from that temperature."""
    t_h = np.asarray(t_h, float)
    x = np.zeros_like(t_h)
    for a in adds:
        m = t_h >= a
        tt = t_h[m] - a
        if head_from is None:
            x[m] += aa_single(T, tt, k_tail)
        else:
            d = np.array(FIG5[head_from])
            h = np.interp(np.minimum(tt, d[-1, 0] / 60), d[:, 0] / 60, d[:, 2])
            k = K_TAIL[T] if k_tail is None else k_tail
            x[m] += np.where(tt <= d[-1, 0] / 60, h, d[-1, 2] * np.exp(-k * (tt - d[-1, 0] / 60)))
    return x


def empirical():
    Km, kmax = 0.39, 39.0
    protos = {"A_plus1_24h": dict(t=24.0, adds=(0.0,), n=1, p_RT=0.65, p_C=0.87),
              "B_plus7_5d": dict(t=120.0, adds=(0.0, 24.0, 48.0, 72.0, 96.0), n=7, p_RT=0.07, p_C=0.60)}
    yld = lambda l, n: float(poisson.sf(n - 1, l))
    ts = np.linspace(0, 120, 12001)
    E_app = 18.4

    def lam(xc, fT, s, Kd):
        c = xc / 2                                   # 10 mM of the letter vs 20 mM AMP in Fig. 5
        r = s * kmax * c / (Kd + c) * fT
        return np.concatenate([[0], np.cumsum(0.5 * (r[1:] + r[:-1]) * np.diff(ts))])

    def run(xc_cold, fT, Kd_cold=Km):
        res = {}
        for pn, P in protos.items():
            xr = aa_curve(T_RT, ts, P["adds"])
            i = int(np.searchsorted(ts, P["t"]))
            unit = lam(xr, 1.0, 1.0, Km)[i]
            lo, hi = 1e-9, 1e4
            for _ in range(200):
                mid = math.sqrt(lo * hi)
                lo, hi = (mid, hi) if yld(mid * unit, P["n"]) < P["p_RT"] else (lo, mid)
            s = math.sqrt(lo * hi)
            xc = xc_cold(P["adds"])
            res[pn] = round(yld(lam(xc, fT, s, Kd_cold)[i], P["n"]), 3)
        return res
    sc = {}
    cold = lambda adds, **kw: aa_curve(T_C, ts, adds, **kw)
    for kt_name, kt in (("tail0.057", K_TAIL[T_C]), ("tail0.02", 0.02)):
        sc[f"cold A*A ({kt_name}), E_ext 18.4 (preactivated control)"] = run(lambda a: cold(a, k_tail=kt), arr_ratio(E_app))
        sc[f"cold A*A ({kt_name}), E_ext 0"] = run(lambda a: cold(a, k_tail=kt), 1.0)
        sc[f"cold A*A ({kt_name}), E_ext 50"] = run(lambda a: cold(a, k_tail=kt), arr_ratio(50))
        for dHb in (-50, -100):
            Kc = Km * math.exp(dHb * 1000 / R * (1 / (T_C + 273.15) - 1 / (T_RT + 273.15)))
            sc[f"cold A*A ({kt_name}), E_ext 18.4, K_m dH_bind {dHb}"] = run(lambda a: cold(a, k_tail=kt), arr_ratio(E_app), Kc)
    sc["RT A*A, E_ext 18.4 (extension term only)"] = run(lambda a: aa_curve(T_RT, ts, a), arr_ratio(E_app))
    sc["RT first 400 min + cold tail 0.057, E_ext 0 (persistence only)"] = run(lambda a: aa_curve(T_C, ts, a, head_from=T_RT), 1.0)
    sc["cold first 400 min + RT tail, E_ext 0 (formation only)"] = run(lambda a: aa_curve(T_RT, ts, a, head_from=T_C), 1.0)
    return sc


if __name__ == "__main__" and os.environ.get("COLD_EMPIRICAL", "1") == "1":
    em = empirical()
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cold_gap.json")
    d = json.load(open(p))
    d["empirical"] = em
    json.dump(d, open(p, "w"), indent=1)
    print("EMPIRICAL (obs cold: A 0.87, B 0.60; RT: A 0.65, B 0.07)")
    for k, v in em.items():
        print(f"  {k:70s} {v}")
