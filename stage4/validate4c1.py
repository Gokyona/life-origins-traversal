"""Stage 4c-1 validation (spec section 4) and cost per run.

  (a) chain-length distribution at steady state (fixed M*, no flow, no cycling) versus the
      analytic steady state of the same kinetics (linear mean-field master equation, solved
      exactly), at 3 growth/hydrolysis ratios g/h = 0.5, 1, 1.5; also versus the best
      geometric (Flory) fit, for reference. Flory's most-probable distribution belongs to
      reversible step-growth (any two chains join, any bond breaks); here growth is monomer
      addition only, so the stationary distribution is not geometric.
  (b) Arrhenius: hydrolysis of a fixed pool of dimers, slope of log10 k vs T (0.07 per C) and
      the reference rate 1.4e-7 /min at 37 C, pH 7.4, 5 mM Mg, 0.25 M; peptide half-lives
  (c) wet-dry cycling versus constant wet at the same mean water activity
  (d) washout (tau = 1 h) -> E2 fails; no feedstock -> E1 fails
  (f) the inactive variable (amphiphile fraction) leaves results bit-identical
  (g) cost per run at screening and full fidelity
Writes validate4c1.json (atomic); usage: python validate4c1.py [out.json] [letters, e.g. abcdf].
"""
import json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c1 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402

dev = torch.device("xpu")
guard = thermal.Guard(job="4c1 validation")


def gen(seed):
    g = torch.Generator(device=dev)
    g.manual_seed(seed)
    return g


def base_u(B):
    return torch.full((B, chem.DIM), 0.5, device=dev)


def sim(B, override, seed, fid="screen", **kw):
    o = chem.simulate(base_u(B), gen(seed), fid, record=True, override=override, **kw)
    guard.check()
    return o


CLEAN = {"tau": 1e12, "depth": 0.0, "k_UV": 0.0, "mineral_type": 0.0, "thermo": 1.0, "evap": 1.0,
         "x_nuc": 1.0, "feed": 0.0, "avail_time": 1e12, "P": 100.0, "mix": 0.0}


def meanfield_steady(g, h, I=1.0):
    """Steady state n_L (L = 2..L_MAX) of growth g per chain end, scission h per bond (both
    fragments kept, fragments of length 1 leave), initiation I dimers per unit time."""
    Lm = chem.L_MAX
    n = Lm - 1
    A = np.zeros((n, n))
    for i, L in enumerate(range(2, Lm + 1)):
        if L >= 3:
            A[i, i - 1] += g
        if L < Lm:
            A[i, i] -= g
        A[i, i] -= (L - 1) * h
        for m in range(L + 1, Lm + 1):
            A[i, m - 2] += 2 * h
    b = np.zeros(n)
    b[0] = -I
    return np.linalg.solve(A, b)


def test_a():
    h = 0.01
    # with M* held fixed, a stationary state exists only for g < 2h (above that, fragments of
    # >= 2 act as new growing chains and the chain number grows without bound); see report
    gs = np.array([0.005, 0.01, 0.015])
    N_target = 200.0
    tot_per_I = np.array([meanfield_steady(g, h).sum() for g in gs])
    I0 = 1.0 * (gs * chem.V) ** 2 / (2 * chem.V)                      # k_c = 1, c_eff = 1
    init_factor = N_target / tot_per_I / I0
    ov = dict(CLEAN, k_c=1.0, clamp_mstar=torch.tensor(gs * chem.V), k_h_fixed=h,
              init_factor=torch.tensor(init_factor), t_end=3000.0)
    t = time.time()
    o = sim(3, ov, 11, fid=(256, 2048, 2000))
    secs = time.time() - t
    rows = []
    for i, g in enumerate(gs):
        hist = o["length_hist"][i].numpy()[2:]
        P = hist / hist.sum()
        an = meanfield_steady(g, h)
        Pa = an / an.sum()
        tv = 0.5 * np.abs(P - Pa).sum()
        # best geometric P(L) ~ p^(L-2) on L = 2..L_MAX by maximum likelihood of the mean
        Ls = np.arange(2, chem.L_MAX + 1)
        best = None
        for pp in np.linspace(0.01, 0.999, 999):
            Pg = pp ** (Ls - 2)
            Pg /= Pg.sum()
            d = 0.5 * np.abs(P - Pg).sum()
            if best is None or d < best[1]:
                best = (pp, d)
        rows.append({"g_per_chain_per_h": g, "h_per_bond_per_h": h, "g_over_h": g / h,
                     "chains_simulated": int(hist.sum()), "mean_len_sim": float((Ls * P).sum()),
                     "mean_len_analytic": float((Ls * Pa).sum()), "TV_vs_analytic": float(tv),
                     "best_geometric_p": float(best[0]), "TV_vs_best_geometric": float(best[1]),
                     "P_sim_L2_10": P[:9].round(4).tolist(), "P_analytic_L2_10": Pa[:9].round(4).tolist(),
                     "overflow_events": float(o["overflow_events"][i]), "stationary": bool(o["stationary"][i])})
    return {"rows": rows, "pass": all(r["TV_vs_analytic"] < 0.05 for r in rows), "seconds": secs}


def test_b():
    out = {}
    nuc_T = [25.0, 37.0, 50.0, 65.0, 80.0]
    pep_T = [25.0, 60.0, 100.0]
    rows = []
    for chem_name, Ts in (("nucleotide", nuc_T), ("peptide", pep_T)):
        B = len(Ts)
        T = torch.tensor(Ts, dtype=torch.float64)
        p = {"pH": torch.full((B,), 7.4, dtype=torch.float64), "Mg": torch.full((B,), 0.005, dtype=torch.float64),
             "salinity": torch.full((B,), 0.25, dtype=torch.float64), "s_exp": torch.zeros(B, dtype=torch.float64)}
        k_exp = chem.k_nuc(T, p) if chem_name == "nucleotide" else chem.k_pep(T)
        t_end = math.log(2) / k_exp
        ov = dict(CLEAN, x_nuc=1.0 if chem_name == "nucleotide" else 0.0, k_c=0.0, pH=7.4, Mg=0.005,
                  salinity=0.25, s_exp=0.0, T=T, t_end=t_end)
        n0 = 500
        o = sim(B, ov, 21 if chem_name == "nucleotide" else 22, fid=(64, 512, 300), init={"dimers": n0})
        for i, Tv in enumerate(Ts):
            left = float(o["length_hist"][i][2])
            S = left / (64 * n0)
            k = -math.log(S) / float(o["t_reached"][i])
            rows.append({"chem": chem_name, "T_C": Tv, "survival": S, "k_sim_per_h": k,
                         "reached_end": bool(o["reached_end"][i]),
                         "k_expected_per_h": float(k_exp[i]), "ratio": k / float(k_exp[i]),
                         "half_life_sim_h": math.log(2) / k})
    nuc = [r for r in rows if r["chem"] == "nucleotide"]
    slope = np.polyfit([r["T_C"] for r in nuc], [math.log10(r["k_sim_per_h"]) for r in nuc], 1)[0]
    k37 = next(r["k_sim_per_h"] for r in nuc if r["T_C"] == 37.0) / 60
    pep = {r["T_C"]: r for r in rows if r["chem"] == "peptide"}
    out = {"rows": rows, "log10_slope_per_C": float(slope), "slope_rel_err": float(abs(slope - 0.07) / 0.07),
           "k37_per_min": k37, "k37_ratio_to_1.4e-7": k37 / 1.4e-7,
           "peptide_half_life_25C_years": pep[25.0]["half_life_sim_h"] / 8766,
           "peptide_half_life_100C_weeks": pep[100.0]["half_life_sim_h"] / 168}
    out["pass"] = out["slope_rel_err"] < 0.10 and 1 / 1.5 <= out["k37_ratio_to_1.4e-7"] <= 1.5
    return out


CYCLE_BASE = {"T": 50.0, "P": 10.0, "pH": 7.0, "Mg": 0.01, "salinity": 0.25, "s_exp": 0.0, "x_nuc": 1.0,
              "feed": 1.0, "tau": 1e4, "mix": 0.01, "thermo": 1.0, "mineral_type": 0.0, "avail_time": 1e8,
              "k_c": 1.0, "Y_UV": 0.01, "UV": 1.0, "Y_H2": 0.0, "H2": 1e-6, "Ea_deact": 70.0, "k_UV": 0.0,
              "evap": 10.0, "f_dry": 10.0, "aw_dry": 0.05, "period": 24.0}


CYCLE_SETTINGS = {
    # saturated: growth far faster than hydrolysis, chains near L_MAX in both arms (no contrast expected)
    "saturated_50C": {},
    # hydrolysis-limited: hotter, slower condensation
    "hydrolysis_limited_80C": {"T": 80.0, "k_c": 0.1, "feed": 0.3},
    "hydrolysis_limited_90C": {"T": 90.0, "pH": 7.5, "k_c": 0.03, "feed": 0.3},
}


def test_c():
    out = {}
    for i, (name, extra) in enumerate(CYCLE_SETTINGS.items()):
        out[name] = _cycle_pair(dict(CYCLE_BASE, **extra), 31 + 10 * i)
        out[name]["overrides"] = extra
    lim = [v for k, v in out.items() if k.startswith("hydrolysis_limited")]
    out["pass"] = all(v["cycling_longer"] for v in lim)
    out["note"] = ("pass requires cycling > constant wet in the hydrolysis-limited settings; "
                   "the saturated setting is reported for completeness")
    return out


def _cycle_pair(base, seed):
    reps = 4
    aw_mean = (1 + base["aw_dry"]) / 2
    cyc = sim(reps, dict(base, depth=1.0, t_end=2000.0, no_cycle_average=True), seed)
    wet = sim(reps, dict(base, depth=0.0, t_end=2000.0, hydro_mult=aw_mean), seed + 1)
    def summ(o):
        h = o["length_hist"][:, 2:].double()
        Ls = torch.arange(2, chem.L_MAX + 1, dtype=torch.float64)
        ml = (h * Ls).sum(1) / h.sum(1).clamp(min=1)
        return {"mean_len": ml.tolist(), "mean_len_avg": float(ml.mean()),
                "chains_per_run": h.sum(1).tolist(), "frac_sites_L8": o["frac_sites_L8"].tolist(),
                "E2": o["E2"].tolist(), "conc_mstar": o["conc_mstar"].tolist()}
    c, w = summ(cyc), summ(wet)
    return {"mean_water_activity": aw_mean, "cycling": c, "constant_wet": w,
            "cycling_longer": c["mean_len_avg"] > w["mean_len_avg"]}


def test_d():
    reps = 4
    base = dict(CYCLE_BASE, depth=1.0)
    o_base = sim(reps, base, 41)
    o_wash = sim(reps, dict(base, tau=1.0), 42)
    o_nofeed = sim(reps, dict(base, feed=0.0), 43)
    f = lambda o: {"E1": o["E1"].tolist(), "E2": o["E2"].tolist(), "conc_mstar": o["conc_mstar"].tolist(),
                   "frac_sites_L8": o["frac_sites_L8"].tolist(), "t_end_h": o["t_end"].tolist()}
    return {"base": f(o_base), "tau_1h": f(o_wash), "no_feed": f(o_nofeed),
            "pass": (not any(o_wash["E2"].tolist())) and (not any(o_nofeed["E1"].tolist()))}


def test_f():
    u = torch.rand(256, chem.DIM, generator=torch.Generator().manual_seed(7)).to(dev)
    o1 = chem.simulate(u, gen(51), "screen", record=True)
    u2 = u.clone()
    u2[:, chem.NAMES.index("amph")] = torch.rand(256, generator=torch.Generator().manual_seed(8)).to(dev)
    o2 = chem.simulate(u2, gen(51), "screen", record=True)
    same = all(torch.equal(o1[k], o2[k]) for k in ("E1", "E2", "length_hist", "conc_mstar", "frac_sites_L8"))
    guard.check()
    return {"identical_with_amph_changed": same, "pass": same}


def test_g():
    from scipy.stats import qmc
    rows = []
    for fid, B in (("screen", 512), ("full", 64)):
        u = torch.tensor(qmc.Sobol(chem.DIM, scramble=True, seed=77).random(B), dtype=torch.float32, device=dev)
        torch.xpu.synchronize()
        t = time.time()
        o = chem.simulate(u, gen(61), fid, record=True)
        torch.xpu.synchronize()
        s = time.time() - t
        rows.append({"fidelity": fid, "sites_slots_steps": chem.FIDELITY[fid], "batch": B, "seconds": s,
                     "ms_per_run": 1000 * s / B, "steps": int(o["steps"][0]),
                     "E1": float(o["E1"].float().mean()), "E2": float(o["E2"].float().mean()),
                     "reached_end": float(o["reached_end"].float().mean())})
        guard.check()
    return {"rows": rows}


def main():
    """python validate4c1.py [out.json] [test letters, e.g. abcdf]"""
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "validate4c1.json")
    which = sys.argv[2] if len(sys.argv) > 2 else "abcdfg"
    res = {}
    for name, fn in (("a_steady_state", test_a), ("b_arrhenius", test_b), ("c_cycling", test_c),
                     ("d_washout_nofeed", test_d), ("f_inactive_variable", test_f), ("g_cost", test_g)):
        if name[0] not in which:
            continue
        t = time.time()
        res[name] = fn()
        res[name]["test_seconds"] = time.time() - t
        print(name, "pass" if res[name].get("pass", True) else "FAIL", f"({time.time()-t:.0f}s)", flush=True)
        ckpt.atomic_json(res, out)
    res["thermal"] = guard.summary()
    ckpt.atomic_json(res, out)
    print(json.dumps({k: v.get("pass") for k, v in res.items() if isinstance(v, dict)}))


if __name__ == "__main__":
    main()
