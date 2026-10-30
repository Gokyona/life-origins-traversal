"""P3-1 cold-gap acceptance test (replaces validation test 6): one parameter set must reproduce
  (A) >= +1 at 24 h: 65% (22 C) / 87% (-3 C), reagent 200 mM once, CMP + GMP 10 mM each;
  (B) >= +7 at 5 d, reagent 200 mM every 24 h, NMPs 10 mM each (4 letters): 7% / 60%;
  (C) pre-activated control, >= +7 at 24 h, *N 10 mM each, no reagent: 42% (22 C) / 8% (+3 C);
within +-5 percentage points (stated uncertainty: gel quantification and figure reading).
Chemistry: the reagent-pool model calibrated on Aitken Fig. 5 (reagent_fit.json, fit B; one-experiment calibration),
equations as in reagent_fit.py / chem_p31 (cold-gap version). Copying: per-letter Poisson rate r = s k_max c / (K_m + c)
f_T(E_ext), C*C on G values (K_m 0.39 mM, k_max 39 /h, Fang 2025), c = blocks of the required letter = X / (number of
letters in the pool); the rate scale s is fitted at 22 C separately for each protocol (different templates, helpers
and primers: system-specific), the cold value is then predicted. E_ext is scanned over its approved band 0-20 kJ/mol.
Residual: the per-letter rate factor that would be needed at -3 C (or +3 C) to hit the observed cold yield.
    python cold_accept.py   -> cold_accept.json
"""
import json, math, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
R = 8.314
T_RT = 22.0
FIT = json.load(open(os.path.join(HERE, "reagent_fit.json")))["B_with_Ea"]
Km, kmax = 0.39, 39.0
TOL = 0.05


def par(T):
    """Rates at T from the fit (Arrhenius through the two calibration temperatures; y logistic in 1/T)."""
    a = lambda E: math.exp(-E * 1000 / R * (1 / (T + 273.15) - 1 / (T_RT + 273.15)))
    rt = FIT["RT"]
    y0 = rt["y"]
    y = 1 / (1 + (1 - y0) / y0 * math.exp(FIT["dE_y_kJ"] * 1000 / R * (1 / (T + 273.15) - 1 / (T_RT + 273.15))))
    return rt["k_act"] * a(FIT["E_activation_kJ"]), rt["k_r"] * a(FIT["E_r_kJ"]), y, rt["k_h"] * a(FIT["E_h_kJ"])


def blocks(T, M0, adds, t_end, X0=0.0, n=4801):
    k_act, k_r, y, k_h = par(T)
    ts = np.linspace(0, t_end, n)
    z = [0.0, M0, X0, 0.0]
    X = np.zeros(n)
    edges = list(adds) + [t_end]
    if not adds:
        edges = [0.0, t_end]

    def f(t, z):
        Rg, M, Xb, W = z
        v = k_act * Rg * M
        return [-k_r * Rg - v, -v + k_h * Xb, (1 - y) * v - k_h * Xb, y * v]
    for i in range(len(edges) - 1):
        if adds:
            z[0] += 200.0
        a, b = edges[i], edges[i + 1]
        m = (ts >= a) & ((ts < b) if i < len(edges) - 2 else (ts <= b))
        sol = solve_ivp(f, (a, b), z, rtol=1e-9, atol=1e-12, method="LSODA", dense_output=True)
        X[m] = sol.sol(ts[m])[2]
        z = list(sol.sol(b))
    return ts, X


def Lam(ts, c, fT):
    r = kmax * c / (Km + c) * fT
    return float(np.sum(0.5 * (r[1:] + r[:-1]) * np.diff(ts)))


PROT = {
    "A_plus1_24h": dict(M0=20.0, letters=2, adds=[0.0], t=24.0, n=1, X0=0.0, T_cold=-3.0, p_RT=0.65, p_cold=0.87),
    "B_plus7_5d": dict(M0=40.0, letters=4, adds=[0.0, 24.0, 48.0, 72.0, 96.0], t=120.0, n=7, X0=0.0, T_cold=-3.0,
                       p_RT=0.07, p_cold=0.60),
    "C_preactivated_plus7_24h": dict(M0=0.0, letters=1, adds=[], t=24.0, n=7, X0=10.0, T_cold=3.0, p_RT=0.42, p_cold=0.08),
}


def lam_for(p, n):
    lo, hi = 1e-9, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if poisson.sf(n - 1, mid) < p else (lo, mid)
    return math.sqrt(lo * hi)


def run(E_ext):
    out = {}
    for name, P in PROT.items():
        fT = math.exp(-E_ext * 1000 / R * (1 / (P["T_cold"] + 273.15) - 1 / (T_RT + 273.15)))
        ts, X = blocks(T_RT, P["M0"], P["adds"], P["t"], P["X0"])
        u = Lam(ts, X / P["letters"], 1.0)
        s = lam_for(P["p_RT"], P["n"]) / u
        tc, Xc = blocks(P["T_cold"], P["M0"], P["adds"], P["t"], P["X0"])
        lc = s * Lam(tc, Xc / P["letters"], fT)
        yc = float(poisson.sf(P["n"] - 1, lc))
        need = lam_for(P["p_cold"], P["n"]) / lc
        out[name] = {"s": s, "yield_RT": P["p_RT"], "yield_cold_model": round(yc, 4), "yield_cold_obs": P["p_cold"],
                     "within_tol": abs(yc - P["p_cold"]) <= TOL, "residual_rate_factor_needed": round(need, 3),
                     "mean_blocks_RT_mM": float(X.mean()), "mean_blocks_cold_mM": float(Xc.mean())}
    return out


def main():
    res = {"fit_used": {k: FIT[k] for k in ("RT", "m3C", "E_activation_kJ", "E_r_kJ", "E_h_kJ", "dE_y_kJ", "rms_mM")},
           "tolerance_abs": TOL, "scan": {}}
    for E in (0.0, 5.0, 10.0, 15.0, 18.4, 20.0):
        res["scan"][str(E)] = run(E)
    best = min(res["scan"], key=lambda e: sum(abs(v["yield_cold_model"] - v["yield_cold_obs"]) for v in res["scan"][e].values()))
    res["best_E_ext"] = float(best)
    res["pass"] = all(v["within_tol"] for v in res["scan"][best].values())
    json.dump(res, open(os.path.join(HERE, "cold_accept.json"), "w"), indent=1)
    print("rates at 22 / -3 / +3 C (k_act, k_r, y, k_h):", [tuple(round(x, 5) for x in par(T)) for T in (22.0, -3.0, 3.0)])
    for E, v in res["scan"].items():
        print(f"E_ext {E:>5}: " + "  ".join(f"{k.split('_')[0]} {x['yield_cold_model']:.3f} (obs {x['yield_cold_obs']}, need x{x['residual_rate_factor_needed']})" for k, x in v.items()))
    print("best E_ext", best, "PASS" if res["pass"] else "FAIL")


if __name__ == "__main__":
    main()
