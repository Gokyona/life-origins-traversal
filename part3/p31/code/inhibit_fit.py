"""Cold-gap follow-up (CPU, analysis only, no model change): do side products inhibit copying at room temperature?
Evidence used (Aitken et al. 2023 main text + SI; curves digitised by eye from Fig. 5 and Fig. S3, +-0.2 mM):
  Table 2 >= +1 at 24 h (CMP + GMP, standard conditions): MBA 65%, acetaldehyde 47%, glycolaldehyde 41%;
  Fig. S3 (AMP 20 mM, RT): A*A and "impurities" (Passerini products, AppA, others) for the three aldehydes.
Hypothesis H: the per-letter extension rate is divided by (1 + I(t) / K_I), I = accumulated side products (mM).
Without H, the copying integral int c/(K_m + c) dt is LARGER for acetaldehyde and glycolaldehyde (more A*A) than for
MBA, while the yields are lower. Fit: rate scale s and K_I on the three RT yields (A*A and I curves per aldehyde, C*C on G
K_m 0.39 mM; c = A*A / 2 for CMP + GMP). Then, with the same K_I:
  A (MBA, -3 C): cold A*A curve (Fig. 5), impurities ~0.1 mM (Fig. 5: near zero);
  B (5 days, daily re-addition, 4 NMPs 10 mM): impurities accumulate per addition: RT +1.5 mM per 20 mM nucleotide
    per addition (MBA, Fig. S3/5 at 400 min, ~1.2-1.5), -3 C +0.16 (Fig. 6: AppA 1.3 mM after 96 h with 40 mM AMP);
    scaled x2 for 40 mM NMP; A*A per addition as in cold_gap.empirical (superposed curves);
  C (pre-activated): no reagent, no side products: unchanged.
E_ext scanned in 0-20 kJ/mol. Output inhibit_fit.json.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ["COLD_EMPIRICAL"] = "0"
import cold_gap as CG  # noqa: E402

Km, kmax = 0.39, 39.0
R = 8.314
pts = lambda L: (np.array([a for a, _ in L]) / 60.0, np.array([b for _, b in L]))
AA = {"MBA": pts([(t, x) for t, _, x in CG.FIG5[22.0]]),
      "acetaldehyde": pts([(0, 0), (10, 2.0), (20, 3.0), (40, 3.6), (60, 3.9), (90, 4.0), (150, 3.6), (200, 3.3),
                           (300, 2.7), (400, 2.2)]),
      "glycolaldehyde": pts([(0, 0), (10, 1.8), (20, 2.6), (40, 3.2), (60, 3.3), (100, 3.2), (150, 2.9), (200, 2.4),
                             (300, 1.7), (400, 1.3)])}
IMP = {"MBA": pts([(0, 0), (50, 0.4), (100, 0.6), (200, 0.9), (300, 1.1), (400, 1.2)]),
       "acetaldehyde": pts([(0, 0), (10, 1.0), (30, 2.2), (60, 3.0), (100, 3.4), (200, 3.9), (300, 4.3), (400, 4.6)]),
       "glycolaldehyde": pts([(0, 0), (10, 1.5), (30, 3.0), (60, 4.4), (100, 5.4), (200, 6.5), (300, 7.3), (400, 8.0)])}
Y_RT = {"MBA": 0.65, "acetaldehyde": 0.47, "glycolaldehyde": 0.41}
ts = np.linspace(0, 120, 12001)


def curve(tp, xp, t, tail=None):
    """interpolate; after the last point: exponential tail (A*A) or constant (impurities)"""
    k = tail if tail is not None else 0.0
    return np.where(t <= tp[-1], np.interp(t, tp, xp), xp[-1] * np.exp(-k * (t - tp[-1])))


def tail_rate(tp, xp):
    i = int(np.argmax(xp))
    return max(math.log(max(xp[i], 1e-9) / max(xp[-1], 1e-9)) / max(tp[-1] - tp[i], 1e-9), 0.0)


def lam(c, I, KI, fT=1.0, t_end=24.0):
    r = kmax * c / (Km + c) * fT / (1 + I / KI if KI else 1.0)
    m = ts <= t_end + 1e-9
    return float(np.sum(0.5 * (r[m][1:] + r[m][:-1]) * np.diff(ts[m])))


def plus1(lmb):
    return 1 - math.exp(-lmb)


def fit():
    curves = {}
    for a in AA:
        tp, xp = AA[a]
        curves[a] = (curve(tp, xp, ts, tail_rate(tp, xp)) / 2, curve(*IMP[a], ts))
    out = {"no_inhibition_relative_integral": {a: lam(*curves[a], None) / lam(*curves["MBA"], None) for a in AA}}
    best = None
    for KI in np.logspace(-2, 2, 801):
        s = -math.log(1 - Y_RT["MBA"]) / lam(*curves["MBA"], KI)
        err = sum((plus1(s * lam(*curves[a], KI)) - Y_RT[a]) ** 2 for a in AA)
        if best is None or err < best[0]:
            best = (err, KI, s)
    _, KI, s = best
    out["K_I_mM"] = KI
    out["fit_yields"] = {a: plus1(s * lam(*curves[a], KI)) for a in AA}
    out["observed"] = Y_RT
    out["rms_points"] = math.sqrt(best[0] / 3)
    return out, KI


def protocols(KI, E_ext):
    fT = math.exp(-E_ext * 1000 / R * (1 / 270.15 - 1 / 295.15))
    res = {}
    # A: one addition, MBA
    cR, IR = CG.aa_curve(22.0, ts, (0.0,)) / 2, curve(*IMP["MBA"], ts)
    cC, IC = CG.aa_curve(-3.0, ts, (0.0,), k_tail=0.02) / 2, np.full_like(ts, 0.1)
    s = -math.log(1 - 0.65) / lam(cR, IR, KI)
    res["A_plus1_24h"] = {"cold_model": plus1(s * lam(cC, IC, KI, fT)), "cold_obs": 0.87}
    # B: daily additions, 4 NMPs (c = A*A-equivalent per letter = 2 x curve / 4), impurities accumulate per addition
    adds = (0.0, 24.0, 48.0, 72.0, 96.0)
    cRB = CG.aa_curve(22.0, ts, adds) * 2 / 4
    cCB = CG.aa_curve(-3.0, ts, adds, k_tail=0.02) * 2 / 4
    IRB = sum(np.where(ts >= a, 2 * 1.5 * np.interp(ts - a, *IMP["MBA"]) / IMP["MBA"][1][-1], 0.0) for a in adds)
    ICB = sum(np.where(ts >= a, 2 * 0.16, 0.0) for a in adds)
    g = lambda s_, c, I, f: float(poisson.sf(6, s_ * lam(c, I, KI, f, 120.0)))
    sB = brentq(lambda s_: g(s_, cRB, IRB, 1.0) - 0.07, 1e-9, 1e4)
    res["B_plus7_5d"] = {"cold_model": g(sB, cCB, ICB, fT), "cold_obs": 0.60}
    # C: unchanged by side products; the cold_accept value at this E_ext (model block decay mapping)
    return res


def main():
    out, KI = fit()
    acc = json.load(open(os.path.join(HERE, "cold_accept.json")))["scan"]
    out["scan"] = {}
    for E in (0.0, 10.0, 15.0, 18.4, 20.0):
        r = protocols(KI, E)
        r["C_preactivated_from_cold_accept"] = acc.get(str(E), {}).get("C_preactivated_plus7_24h", {}).get("yield_cold_model")
        out["scan"][str(E)] = r
    out["scan_KI_x0.5_x2_at_Eext_18.4"] = {f: protocols(KI * f, 18.4) for f in (0.5, 2.0)}
    json.dump(out, open(os.path.join(HERE, "inhibit_fit.json"), "w"), indent=1)
    print("no-inhibition copying integral relative to MBA:", {k: round(v, 3) for k, v in out["no_inhibition_relative_integral"].items()})
    print("K_I =", round(KI, 4), "mM; fitted RT yields", {k: round(v, 3) for k, v in out["fit_yields"].items()}, "obs", Y_RT)
    for E, r in out["scan"].items():
        print(f"E_ext {E:>5}: A {r['A_plus1_24h']['cold_model']:.3f} (0.87)  B {r['B_plus7_5d']['cold_model']:.3f} (0.60)  C {r['C_preactivated_from_cold_accept']} (0.08)")
    for f, r in out["scan_KI_x0.5_x2_at_Eext_18.4"].items():
        print(f"K_I x{f} @18.4: A {r['A_plus1_24h']['cold_model']:.3f}  B {r['B_plus7_5d']['cold_model']:.3f}")


if __name__ == "__main__":
    main()
