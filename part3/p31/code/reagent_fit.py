"""P3-1 cold-gap change (a): calibrate the activating-reagent pool on Aitken et al. 2023 Fig. 5 (one experiment).
Single well-mixed site, the equations used by chem_p31 (reagent pool R, monomer M, activated blocks X = A*A, waste):
  dR/dt = J - k_r(T) R - k_act R M                       (reagent loss + one reagent per reaction with a monomer)
  dM/dt = -k_act R M + k_h(T) X                          (activated-block hydrolysis returns the monomer)
  dX/dt = (1 - y(T)) k_act R M - k_h(T) X                 (share 1 - y activates; share y -> side products, the
  dW/dt = y(T) k_act R M                                  monomer is lost to them: AppA, Passerini, *A in this mapping)
Mapping (one-experiment calibration, stated): X = bridged A*A (mM; the species whose K_m is measured, one block per
consumed nucleotide as in chem_p31); everything else made from AMP (*A, AppA, impurities) is "side products". The
monomeric *A is in fact also an activated species; counting it as side product is the conservative reading.
Data: Fig. 5 digitised by eye (+-0.2 mM), AMP + Nu-pA (unactivated pool) and A*A, 0-400 min, RT (22 C) and -3 C; R(0) = 200 mM, J = 0.
Fit: k_act shared (T-independent, as in the spec) with an optional activation energy E_a (fit B), and per temperature
k_r, y, k_h. Reported: Arrhenius energies from the two temperatures (k_r, k_h) and the logistic temperature
coefficient of y. Output reagent_fit.json.
"""
import json, math, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cold_gap import FIG5, T_RT, T_C  # noqa: E402

# Nu-pA (catalyst adduct, not activated; Fig. 5, read by eye) is added to AMP: the unactivated pool M. Fig. 5 Nu-pA:
# RT 2.5 (10 min), 1.2 (20), 0.4 (30), ~0.1 later; -3 C 3.8 (10), 3.9 (30), 2.5 (60), 1.8 (80), 1.4 (100), 1.0 (130),
# 0.6 (160), 0.3 (200), ~0 later
NUPA = {T_RT: {10: 2.5, 20: 1.2, 30: 0.4}, T_C: {10: 3.8, 30: 3.9, 60: 2.5, 80: 1.8, 100: 1.4, 130: 1.0, 160: 0.6, 200: 0.3}}
FIG5 = {T: [(t, m + NUPA[T].get(t, 0.0), x) for (t, m, x) in v] for T, v in FIG5.items()}

R = 8.314


def rhs(k_act, k_r, y, k_h, J=0.0):
    def f(t, z):
        Rg, M, X, W = z
        v = k_act * Rg * M
        return [J - k_r * Rg - v, -v + k_h * X, (1 - y) * v - k_h * X, y * v]
    return f


def sim(par, t_h, R0=200.0, M0=20.0):
    sol = solve_ivp(rhs(*par), (0, max(t_h)), [R0, M0, 0.0, 0.0], t_eval=t_h, rtol=1e-9, atol=1e-12, method="LSODA")
    return sol.y


def unpack(v, with_Ea):
    k_act = math.exp(v[0])
    Ea = v[1] if with_Ea else 0.0
    o = 2 if with_Ea else 1
    per = {}
    for i, T in enumerate((T_RT, T_C)):
        k_r, y, k_h = math.exp(v[o + 3 * i]), 1 / (1 + math.exp(-v[o + 3 * i + 1])), math.exp(v[o + 3 * i + 2])
        ka = k_act * math.exp(-Ea * 1000 / R * (1 / (T + 273.15) - 1 / (T_RT + 273.15)))
        per[T] = (ka, k_r, y, k_h)
    return per


def fit(with_Ea):
    def res(v):
        per = unpack(v, with_Ea)
        out = []
        for T in (T_RT, T_C):
            d = np.array(FIG5[T])
            Rg, M, X, W = sim(per[T], d[:, 0] / 60)
            out += list(X - d[:, 2]) + list(M - d[:, 1])
        return np.array(out)
    best = None
    for g in ([-6, 0.5, -2, -3, 1, 0.5, -3], [-4, 1, -1, -2, 3, 1, -5], [-8, -1, 0, -1, -1, -2, -2]):
        v0 = [g[0]] + ([0.0] if with_Ea else []) + [g[1], g[2], g[3], g[4], g[5], g[6]]
        lo = [-20] + ([-100] if with_Ea else []) + [-12, -8, -12] * 2
        hi = [5] + ([150] if with_Ea else []) + [5, 8, 3] * 2
        try:
            r = least_squares(res, v0, bounds=(lo, hi))
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    per = unpack(best.x, with_Ea)
    rms = float(np.sqrt(np.mean(best.fun ** 2)))
    E = lambda a, b: -R * math.log(b / a) / (1 / (T_C + 273.15) - 1 / (T_RT + 273.15)) / 1000
    yR, yC = per[T_RT][2], per[T_C][2]
    # y(T) = 1 / (1 + ((1 - y_ref) / y_ref) exp(dE_y/R (1/T - 1/T_ref))): dE_y from the two temperatures
    dEy = R * math.log(((1 - yC) / yC) / ((1 - yR) / yR)) / (1 / (T_C + 273.15) - 1 / (T_RT + 273.15)) / 1000
    return {"with_Ea_activation": with_Ea, "rms_mM": rms,
            "k_act_per_mM_h_at_RT": per[T_RT][0], "E_activation_kJ": (best.x[1] if with_Ea else 0.0),
            "RT": dict(zip(("k_act", "k_r", "y", "k_h"), per[T_RT])), "m3C": dict(zip(("k_act", "k_r", "y", "k_h"), per[T_C])),
            "E_r_kJ": E(per[T_RT][1], per[T_C][1]), "E_h_kJ": E(per[T_RT][3], per[T_C][3]), "dE_y_kJ": dEy,
            "residuals_X_then_M_RT_then_m3C": [round(float(x), 3) for x in best.fun]}


def main():
    out = {"A_no_Ea": fit(False), "B_with_Ea": fit(True)}
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reagent_fit.json")
    json.dump(out, open(p, "w"), indent=1)
    for k, v in out.items():
        print(k, {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "residuals_X_then_M_RT_then_m3C"})


if __name__ == "__main__":
    main()
