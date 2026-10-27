"""Stage 4c-3 validation tests 1-5 (spec docs/stage4c3_spec.md, section 3); test 6 (kill and
resume) is kill_resume_test4c3.sh with run4c3_map.py.

  (1) copy fidelity (unprimed and primed copies; primers are the dimers made by initiation):
      a fixed template, no hydrolysis or washout, M* clamped, T = 25 C (so the
      Arrhenius factor is 1 and r = k_c k_t thermo M*/V is known in closed form). The measured
      error rate of the logged copies (against the complement of the template segment they were
      copied from) equals eps within 3 binomial standard errors; the copy-length distribution
      matches the processivity formula q(l) = r / (r + K0 exp((l - 1) dH/R (1/T_m8 - 1/T)))
      (computed here independently of chem4c3): mean length within 3 standard errors and a
      chi-square test on the length histogram (p > 0.001);
  (2) no copying: k_t = 0 at REF_L2 -> L2_strict never;
  (3) random-copy control: copy_random as the main run at REF_L2 -> L2_strict never;
  (4) error threshold: at REF_L2, eps * L_H = 0.1, 0.3, 1, 3 -> the information share falls and
      L2_strict disappears at eps * L_H = 3;
  (5) strand separation, at an explicit condition with long copies (see test_5; at REF_L2
      copies are too short for duplexes to last a step), dry phase with temperature effects
      only, run length fixed: (a) T = T_m8 - 30 without cycling: templates end up paired and
      copying stalls (paired share > 0.5; share of copy starts in the last quarter < 0.1,
      counting released copies before the M* and slot limits, since slots fill up); (b) cycling
      that crosses T_m8 (T = T_m8 - 5, depth 1: wet T_m8 - 15, dry T_m8 + 5, resolved phases):
      fewer templates paired and copying continues (last-quarter share >= 0.15); (c) T = T_m8
      + 30 without cycling: copies fall off before they are long enough (mean copy length
      < L_H, less than half that in (a), and shorter than in (b)).
  (7) parasite takeover (plan 4c-3 milestone): see test_7.
  (8) null false-positive rate of L2_strict with random letters, 64 reps at two box points: test_8.
  Test 4 was redone at an explicit condition (see test_4): no natural L2 reference exists.
REF_L2 (reference condition where L2 occurs) is the u vector chosen by scan4c3.py (stage 1:
4096 Sobol points of the full space at screening fidelity; stage 2: the L2_strict points rerun
at validation fidelity, 16 reps; the highest L2_strict rate wins), read from scan4c3.json.
Tests 2-5 run at REF_L2 with the overrides stated in each test.

  python validate4c3.py [out.json] [tests, e.g. 12345]
"""
import json, math, os, sys, time
import numpy as np
import torch
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem4c3 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402

dev = torch.device("xpu")
guard = thermal.Guard(job="4c3 validation")
FID = (16, 256, 600)
REPS = 16
ARR80 = math.exp(60000 / 8.314 * (1 / 298.15 - 1 / 353.15))
BASE = dict(T=80.0, pH=7.5, Mg=0.005, salinity=0.25, s_exp=0.0, P=100.0, feed=10.0, tau=1e4, k_c=0.03 / ARR80,
            p0=1e-6, kappa=1e6, L_cat=3.0, lam=0.5, Ea_unc=60.0, r_Ea=0.5, depth=0.0, thermo=1.0, mineral_type=0.0,
            x_nuc=1.0, Y_UV=0.1, UV=1.0, Y_H2=0.0, H2=1e-6, k_UV=0.0, avail_time=1e8,
            L_t=6.0, k_t=1e5, eps=1e-3, dH_bp=35.0, T_m8=50.0)
_scan = ckpt.load_json(os.path.join(HERE, os.environ.get("LOT_REF_SCAN", "scan4c3.json"))) or {}   # REF_L2 source
SCAN_PICK = _scan.get("ref")
REF_U = SCAN_PICK["u"] if SCAN_PICK else [0.5] * chem.DIM
REF_L2 = dict(SCAN_PICK["physical"]) if SCAN_PICK else {}
REF_L2.pop("amph", None)


def run(override, seed, reps=REPS, control=True, fid=FID):
    """REF_L2 (as its u vector) with physical overrides."""
    u = torch.tensor([REF_U] * reps, dtype=torch.float32, device=dev)
    g = torch.Generator(device=dev)
    g.manual_seed(seed)
    o = chem.simulate(u, g, fid, record=False, override=dict(override))
    guard.check()
    strict = torch.zeros(reps, dtype=torch.bool)
    s_ctl = torch.full((reps,), float("nan"), dtype=torch.float64)
    s_m = torch.full((reps,), float("nan"), dtype=torch.float64)
    if control:
        strict, s_m, s_ctl = chem.paired_control_l2(u, seed, fid, o, override=override)
        guard.check()
    f = lambda k: [round(float(x), 5) for x in o[k]]
    return {"L2_strict": int(strict.sum()), "L2_loose": int(o["L2_loose"].sum()), "E1": int(o["E1"].sum()),
            "span": f("L2_span_turnovers"), "share": f("L2_share"), "share_last5": f("L2_share_last5"),
            "share_ctl": [None if math.isnan(x) else round(float(x), 5) for x in s_ctl],
            "S_main_rerun": [None if math.isnan(x) else round(float(x), 5) for x in s_m],
            "copies": f("copies"), "copy_len": f("copy_mean_len"), "copy_err": f("copy_error_rate"),
            "copies_q4": f("copies_q4"), "paired": f("paired_template_share"), "t_end_h": float(o["t_end"][0]),
            "top_len": f("L2_top_lineage_len")}


def mean(x):
    return float(np.mean(x))


# ---------- test 1 ----------

def test_1():
    T, Tm, dH = 25.0, 40.0, 30.0
    k_c, k_t, thermo, ms, eps = 0.01, 1e5, 1e3, 1000.0, 0.05
    tseq = 0b101100111000101101001110          # 24 letters
    ov = dict(BASE, T=T, T_m8=Tm, dH_bp=dH, k_c=k_c, k_t=k_t, thermo=thermo, eps=eps, L_t=8.0, p0=0.0,
              k_h_fixed=0.0, k_UV=0.0, tau=1e12, clamp_mstar=ms, t_end=0.006, copy_log=True,
              init_templates={"n": 32, "L": 24, "seq": tseq})
    u = torch.full((4, chem.DIM), 0.5, device=dev)            # test 1: explicit condition, not REF_L2
    g = torch.Generator(device=dev)
    g.manual_seed(101)
    o = chem.simulate(u, g, FID, override=ov)
    guard.check()
    log = {k: torch.cat([d[k] for d in o["copy_log"]]) for k in o["copy_log"][0]}
    tq, p, l, cq, l0 = log["tseq"], log["p"], log["l"], log["cseq"], log["l0"]
    n = l.numel()
    # errors of the copied letters against the complement of the template segment; the primer
    # letters (primed starts) must be exact complements
    errs, perr = 0, 0
    for j in range(64):
        tb = (torch.bitwise_right_shift(tq, p + j) & 1)
        cb = (torch.bitwise_right_shift(cq, j) & 1)
        errs += int(((l > j) & (j >= l0) & (cb == tb)).sum())   # equal letter = not the complement = error
        perr += int(((l > j) & (j < l0) & (l0 >= 2) & (cb == tb)).sum())
    pos = int((l - l0).sum())
    rate = errs / pos
    se = math.sqrt(eps * (1 - eps) / pos)
    ok_err = abs(rate - eps) <= 3 * se
    # processivity formula, independent of chem4c3: r and k_off from the parameters
    r = k_c * 1.0 * k_t * 1.0 * 1.0 * thermo * ms / 1000.0
    e = dH * 1000 / 8.314 * (1 / (Tm + 273.15) - 1 / (T + 273.15))
    q = lambda j: r / (r + 1e6 * math.exp((j - 1) * e))
    maxl = log["maxl"].tolist()
    exp_mean, var_sum, pmf_tot = 0.0, 0.0, np.zeros(65)
    for m, a in zip(maxl, l0.tolist()):
        surv = [1.0]                                         # P(length >= a + i), i = 0..m - a
        for j in range(a, m):
            surv.append(surv[-1] * q(j))
        pmf = np.array([surv[i] - (surv[i + 1] if i + 1 < len(surv) else 0.0) for i in range(len(surv))])
        if a == 1:
            pmf = pmf[1:] / pmf[1:].sum()                    # unprimed: condition on a released copy (l >= 2)
        ls = np.arange(max(a, 2), m + 1)
        mu = float((pmf * ls).sum())
        exp_mean += mu
        var_sum += float((pmf * ls ** 2).sum()) - mu ** 2
        pmf_tot[max(a, 2):m + 1] += pmf
    obs = l.double().numpy()
    z = (obs.sum() - exp_mean) / math.sqrt(var_sum)
    hist = np.bincount(obs.astype(int), minlength=65)[:65].astype(float)
    # chi-square over bins with expected >= 5 (the tail pooled)
    bins, eb, ob, acc_e, acc_o = [], [], [], 0.0, 0.0
    for j in range(2, 65):
        acc_e += pmf_tot[j]
        acc_o += hist[j]
        if acc_e >= 5:
            eb.append(acc_e); ob.append(acc_o); acc_e = acc_o = 0.0
    if acc_e > 0 and eb:
        eb[-1] += acc_e; ob[-1] += acc_o
    stat = float(sum((o_ - e_) ** 2 / e_ for o_, e_ in zip(ob, eb)))
    pval = float(chi2.sf(stat, max(len(eb) - 1, 1)))
    r_logged = float(log["r"].double().mean())
    ok_len = abs(z) <= 3 and pval > 1e-3
    return {"copies_logged": n, "primed": int((l0 >= 2).sum()), "positions": pos, "errors": errs,
            "primer_letter_mismatches": perr, "error_rate": rate, "eps": eps, "se": se,
            "r_formula": r, "r_logged_mean": r_logged, "dup_e": e,
            "mean_len_obs": float(obs.mean()), "mean_len_exp": exp_mean / n, "z_mean_len": z,
            "chi2": stat, "chi2_bins": len(eb), "chi2_p": pval,
            "len_hist_obs": hist[2:25].tolist(), "len_hist_exp": [round(x, 2) for x in pmf_tot[2:25]],
            "pass_error_rate": bool(ok_err and perr == 0), "pass_length": bool(ok_len),
            "pass": bool(ok_err and ok_len and perr == 0)}


# ---------- tests 2-5 ----------

def test_2():
    r = run({"k_t": 0.0}, 201)
    # the same without adsorption protection in the turnover time (4c-2 formula), for comparison
    chem.TURNOVER_PROTECT = False
    try:
        r_old = run({"k_t": 0.0}, 201)
    finally:
        chem.TURNOVER_PROTECT = True
    return {"run": r, "run_turnover_unprotected": r_old,
            "loose_strict_protected": [r["L2_loose"], r["L2_strict"]],
            "loose_strict_unprotected": [r_old["L2_loose"], r_old["L2_strict"]], "pass": r["L2_strict"] == 0}


def test_3():
    ref = run({}, 301)
    rnd = run({"copy_random": True}, 301)
    return {"ref": ref, "copy_random": rnd, "pass": rnd["L2_strict"] == 0}


def test_4():
    """Error threshold, positive control at an explicit condition (no natural L2 reference
    exists: scans of the full space and of the box found L2_strict only at noise level). A
    master template of 12 letters (32 per site) under resolved cycling across T_m8 (wet T_m8 -
    15, dry T_m8 + 5), M* clamped, no initiation, hydrolysis k_h 0.01 /h per bond and washout
    tau 50 h: the master lives ~7 h and is re-copied about every 0.1 h cycle (R >> 1); run 30 h,
    1200 steps. eps * L_H = 0.1, 0.3, 1, 3. Measured: the information share S (L2 detector) and
    the master-lineage share (chain mass, chemistry 0, L >= L_H, equal to a segment of the
    master or of its complement, at the end); L2_strict with the paired copy_random controls.
    Pass: L2_strict in most reps (>= 9/16) at eps L_H = 0.1, none at eps L_H = 3, and the master
    share at eps L_H = 3 below 0.2 x that at 0.1."""
    Tm = 60.0
    master = 0b101100111000
    # growth slowed 100-fold at the same copying r (k_c 1e-4, k_t 1e7): uncatalysed growth
    # (~10 letters/h at k_c 0.01) otherwise appends random letters to every copy
    # Positive control (redesign approved 26 Oct 2026): many copy cycles per master lifetime
    # (R >> 1): period 0.1 h, k_h 0.01 /h per bond, washout tau 50 h (a 12-mer lives ~7 h, ~70
    # cycles). Growth fix by default (k_c 1e-4, k_t 1e7, outside the layer-C range of k_t);
    # LOT_T4ORIG=1: the original k_c 0.01, k_t 1e5 (quantifies the growth limitation; report only).
    kc4, kt4 = (0.01, 1e5) if os.environ.get("LOT_T4ORIG") == "1" else (1e-4, 1e7)
    # LOT_T4V3=1 (fix after the first positive-control run): a 30 h run is shorter than the
    # 5 turnover times (~14 h each) that the L2 window needs, and growth at 0.1 /h still
    # appended ~1 letter per copy lifetime: run 100 h, 2400 steps, k_c 1e-6, k_t 1e9 (same r)
    v3 = os.environ.get("LOT_T4V3") == "1"
    if v3 and os.environ.get("LOT_T4ORIG") != "1":
        kc4, kt4 = 1e-6, 1e9
    t_end4, steps4 = (100.0, 2400) if v3 else (30.0, 1200)
    base = dict(BASE, T_m8=Tm, dH_bp=35.0, k_c=kc4, k_t=kt4, thermo=1e3, L_t=8.0, p0=0.0,
                k_h_fixed=0.01, k_UV=0.0, tau=50.0, clamp_mstar=1000.0, t_end=t_end4, init_factor=0.0,
                evap=1.0, f_dry=1.0, aw_dry=1.0, Ea_unc=60.0, r_Ea=0.5, T=Tm - 5, depth=1.0,
                period=0.1, no_cycle_average=True, return_state=True,
                init_templates={"n": 32, "L": 12, "seq": master})
    segs = set()
    for v in (master, (~master) & 0xFFF):
        for a in range(12):
            for b in range(a + chem.L_H, 13):
                segs.add((b - a, (v >> a) & ((1 << (b - a)) - 1)))
    out = {}
    fid = (16, 1024, steps4)
    for x in (0.1, 0.3, 1.0, 3.0):
        ov = dict(base, eps=x / chem.L_H)
        u = torch.full((REPS, chem.DIM), 0.5, device=dev)
        g = torch.Generator(device=dev)
        g.manual_seed(401)
        o = chem.simulate(u, g, fid, override=dict(ov))
        guard.check()
        st = o["state"]
        Ls, sq, ch = st["L"], st["seq"], st["chem"]
        m = (ch == 0) & (Ls >= chem.L_H)
        hit = torch.zeros_like(m)
        for (l, v) in segs:
            hit |= m & (Ls == l) & (sq == v)
        msh = float((Ls * hit).sum() / (Ls * m).sum().clamp(min=1))
        # diagnostics: sites where a master segment reaches the abundance threshold (>= 2 copies
        # of one sequence), and S per checkpoint (largest share over sites, mean over reps)
        ab_sites = 0
        for r_ in range(Ls.shape[0]):
            kk = (Ls[r_][hit[r_]] * (1 << 40) + sq[r_][hit[r_]]).tolist()
            ab_sites += int(any(kk.count(x) >= 2 for x in set(kk)))
        s_chk = o["L2_share_matrix"].max(1).values.mean(0).tolist()
        ov.pop("return_state")
        strict, s_m, s_c = chem.paired_control_l2(u, 401, fid, o, override=ov)
        guard.check()
        out[f"epsL{x:g}"] = {"eps": x / chem.L_H, "master_share": msh, "S_mean": float(o["L2_share_last5"].mean()),
                             "L2_loose": int(o["L2_loose"].sum()), "L2_strict": int(strict.sum()),
                             "margin_mean": float(torch.nan_to_num(s_m - s_c, nan=0.0)[o["L2_loose"]].mean())
                             if bool(o["L2_loose"].any()) else None,
                             "master_abundant_sites": ab_sites, "sites": int(Ls.shape[0]),
                             "S_per_checkpoint": [round(x, 4) for x in s_chk],
                             "tau_h_mean": float(o["L2_tau_h"].mean()), "t_end": float(o["t_end"][0]),
                             "copies": float(o["copies"].sum()), "copy_len": float((o["copy_mean_len"] * o["copies"]).sum())
                             / max(float(o["copies"].sum()), 1.0), "copy_err": float(o["copy_error_rate"].mean())}
    keys = ["epsL0.1", "epsL0.3", "epsL1", "epsL3"]
    ms = [out[k]["master_share"] for k in keys]
    sv = [out[k]["S_mean"] for k in keys]
    st_ = [out[k]["L2_strict"] for k in keys]
    ok = st_[0] >= REPS // 2 + 1 and st_[-1] == 0 and ms[-1] < 0.2 * ms[0]
    out["summary"] = {"eps_times_LH": [0.1, 0.3, 1, 3], "master_share": ms, "S_mean": sv, "L2_strict": st_,
                      "k_c": kc4, "k_t": kt4, "L2_KEY_LEN": chem.L2_KEY_LEN, "t_end": t_end4, "steps": steps4}
    out["pass"] = bool(ok)
    return out


def test_5():
    """Strand separation at an explicit condition (as test 1; REF_L2's copies are ~2.5 letters,
    too short for any duplex to outlive a step there): templates of 24 letters, M* clamped,
    no hydrolysis, washout or initiation, T_m8 = 40 C, dH_bp = 30 kJ/mol, r = 1e6 /h at 25 C,
    L_t = 12, 4096 slots per site, dry phase with temperature effects only."""
    Tm = 40.0
    base = dict(BASE, T_m8=Tm, dH_bp=30.0, k_c=0.01, k_t=1e5, thermo=1e3, eps=1e-3, L_t=12.0, p0=0.0,
                k_h_fixed=0.0, k_UV=0.0, tau=1e12, clamp_mstar=1000.0, t_end=0.006, init_factor=0.0,
                evap=1.0, f_dry=1.0, aw_dry=1.0, Ea_unc=60.0, r_Ea=0.5,
                init_templates={"n": 32, "L": 24, "seq": 0b101100111000101101001110})
    rows = {"a_below": dict(base, T=Tm - 30, depth=0.0),
            "b_cycling": dict(base, T=Tm - 5, depth=1.0, period=0.006 / 8, no_cycle_average=True),
            "c_above": dict(base, T=Tm + 30, depth=0.0)}
    out, s = {}, {}
    for name, ov in rows.items():
        u = torch.full((4, chem.DIM), 0.5, device=dev)
        g = torch.Generator(device=dev)
        g.manual_seed(502)
        o = chem.simulate(u, g, (16, 4096, 600), record=True, override=dict(ov))
        guard.check()
        hist = o["length_hist"][:, 12:].sum().item()
        cp = float(o["copies"].sum())
        tries = float(o["copy_tries"].sum())
        s[name] = {"paired_share": float(o["paired_template_share"].mean()), "copies": cp, "tries": tries,
                   "q4_share": float(o["copy_tries_q4"].sum()) / max(tries, 1.0),
                   "added_per_copy": float((o["copy_mean_added"] * o["copies"]).sum()) / max(cp, 1.0),
                   "copy_len": float((o["copy_mean_len"] * o["copies"]).sum()) / max(cp, 1.0),
                   "dropped": float(o["copies_dropped"].sum()), "chains_ge_12": hist}
    a, b, c = s["a_below"], s["b_cycling"], s["c_above"]
    ok_a = a["paired_share"] > 0.5 and a["q4_share"] < 0.1
    ok_b = b["paired_share"] < a["paired_share"] and b["q4_share"] >= 0.15
    # "long enough" = the lineage length L_H of the L2 detector. (A first criterion, fewer than
    # 1 letter added per copy, failed with 1.30: at T_m8 + 30 the Arrhenius factor raises r about
    # 5-fold, so q(1) ~ 0.83 while q(l) still falls by e^-1 per letter.)
    ok_c = c["copy_len"] < chem.L_H and c["copy_len"] < 0.5 * a["copy_len"] and c["copy_len"] < b["copy_len"]
    out["summary"] = s
    out["pass_a"], out["pass_b"], out["pass_c"] = bool(ok_a), bool(ok_b), bool(ok_c)
    out["pass"] = bool(ok_a and ok_b and ok_c)
    return out


def _species_share(st, specs):
    """Per site: count of chains equal to each species or its complement -> [sites, n_species]."""
    L, sq, ch = st["L"], st["seq"], st["chem"]
    cnt = []
    for (l, v) in specs:
        comp = (~v) & ((1 << l) - 1)
        m = (ch == 0) & (L == l) & ((sq == v) | (sq == comp))
        cnt.append(m.sum(1))
    return torch.stack(cnt, 1).double()


def test_7():
    """Parasite takeover: two seeded template species, 16 copies each per site, compete for
    activated monomers and slots under resolved temperature cycling across T_m8, with
    hydrolysis on (k_h 0.05 /h per bond; no washout, no initiation, M* clamped). Measured: the
    share of the short species among chains that are full copies of either species (or their
    complements), per site at the end. (a) short (L 8) vs long (L 24): the short share rises
    from 0.5 to >= 0.9; (b) control, two species of equal length (16 and 16): the share stays
    within 0.5 +/- 3 standard errors (over sites); (c) as (a) with eps = 0.02 per letter."""
    Tm = 60.0
    sA, sB = 0b10110011, 0b101101001110100110010110             # 8 and 24 letters
    sC, sD = 0b1011001110001011, 0b0110100111010010             # 16 and 16 letters
    base = dict(BASE, T_m8=Tm, dH_bp=35.0, k_c=0.01, k_t=1e5, thermo=1e3, eps=1e-3, L_t=8.0, p0=0.0,
                k_h_fixed=0.05, k_UV=0.0, tau=1e12, clamp_mstar=1000.0, t_end=2.0, init_factor=0.0,
                evap=1.0, f_dry=1.0, aw_dry=1.0, Ea_unc=60.0, r_Ea=0.5, T=Tm - 5, depth=1.0,
                period=2.0 / 16, no_cycle_average=True, return_state=True)
    rows = {"a_short_vs_long": (dict(base, init_templates={"n": [16, 16], "L": [8, 24], "seq": [sA, sB]}),
                                [(8, sA), (24, sB)]),
            "b_equal_length": (dict(base, init_templates={"n": [16, 16], "L": [16, 16], "seq": [sC, sD]}),
                               [(16, sC), (16, sD)]),
            "c_short_vs_long_eps0.02": (dict(base, eps=0.02, init_templates={"n": [16, 16], "L": [8, 24],
                                                                           "seq": [sA, sB]}), [(8, sA), (24, sB)])}
    out, s = {}, {}
    for name, (ov, specs) in rows.items():
        u = torch.full((4, chem.DIM), 0.5, device=dev)
        g = torch.Generator(device=dev)
        g.manual_seed(701)
        o = chem.simulate(u, g, (16, 4096, 600), override=dict(ov))
        guard.check()
        c = _species_share(o["state"], specs)
        tot = c.sum(1)
        sh = c[:, 0] / tot.clamp(min=1)
        ok_sites = tot > 0
        shv = sh[ok_sites]
        s[name] = {"first_species_share_mean": float(shv.mean()) if shv.numel() else None,
                   "se": float(shv.std() / math.sqrt(max(shv.numel(), 1))) if shv.numel() > 1 else None,
                   "sites_with_lineages": int(ok_sites.sum()), "full_copies_per_site_mean": float(tot.mean()),
                   "copies": float(o["copies"].sum()), "copy_len": float((o["copy_mean_len"] * o["copies"]).sum())
                   / max(float(o["copies"].sum()), 1.0)}
    a, b, c = s["a_short_vs_long"], s["b_equal_length"], s["c_short_vs_long_eps0.02"]
    ok_a = a["first_species_share_mean"] is not None and a["first_species_share_mean"] >= 0.9
    ok_b = b["first_species_share_mean"] is not None and abs(b["first_species_share_mean"] - 0.5) <= 3 * max(b["se"] or 0, 0.02)
    ok_c = c["first_species_share_mean"] is not None and c["first_species_share_mean"] >= 0.9
    out["summary"] = s
    out["pass_a"], out["pass_b"], out["pass_c"] = bool(ok_a), bool(ok_b), bool(ok_c)
    out["pass"] = bool(ok_a and ok_b and ok_c)
    return out


def _box_point(i):
    """u vector of point i of the box scan (scan4c3.py --box, Sobol seed 2026)."""
    from scipy.stats import qmc
    X = qmc.Sobol(d=chem.DIM, scramble=True, seed=2026).random(4096).astype(np.float32)
    ix = {n: chem.NAMES.index(n) for n in ("T", "T_m8", "depth", "k_t", "thermo", "x_nuc")}
    Tm = 10 + 80 * X[:, ix["T_m8"]]
    X[:, ix["T"]] = np.clip((Tm - 15 * X[:, ix["T"]]) / 150.0, 0.0, 1.0)
    X[:, ix["depth"]] = 0.5 + 0.5 * X[:, ix["depth"]]
    X[:, ix["k_t"]] = 0.75 + 0.25 * X[:, ix["k_t"]]
    X[:, ix["thermo"]] = 0.375 + 0.625 * X[:, ix["thermo"]]
    X[:, ix["x_nuc"]] = 0.5 + 0.5 * X[:, ix["x_nuc"]]
    return [float(v) for v in X[i]]


def test_8():
    """Null false-positive rate of L2_strict with random letters (eps = 0.5, eps L_H = 3), 64 reps,
    at the box reference #3493 and at box point #289 (another L2_loose point), with both
    copy_random controls for every candidate. Pass: <= 2/64 (1/32) at the margin in use at both
    points; reported: the smallest margin that meets the target at both."""
    out = {}
    need = []
    for name, i in (("p3493", 3493), ("p289", 289)):
        u = torch.tensor([_box_point(i)] * 64, dtype=torch.float32, device=dev)
        ov = {"eps": 3.0 / chem.L_H}
        g = torch.Generator(device=dev)
        g.manual_seed(801)
        o = chem.simulate(u, g, FID, override=dict(ov))
        guard.check()
        strict, s_m, s_c = chem.paired_control_l2(u, 801, FID, o, override=ov, both=True)
        guard.check()
        mk = o["L2_loose"] & ~torch.isnan(s_m)
        d = sorted((s_m - s_c)[mk].tolist(), reverse=True)
        m_need = (d[2] + 1e-6) if len(d) > 2 else 0.0          # >= 3 diffs above -> margin must exceed the 3rd
        need.append(m_need)
        out[name] = {"loose": int(o["L2_loose"].sum()), "strict": int(strict.sum()), "diffs_desc": d[:10],
                     "smallest_margin_le_2_of_64": m_need}
    out["margin_used"] = chem.L2_MARGIN
    out["L2_KEY_LEN"] = chem.L2_KEY_LEN
    out["smallest_margin_both_points"] = max(need)
    out["pass"] = bool(all(out[k]["strict"] <= 2 for k in ("p3493", "p289")))
    return out


def test_9():
    """Error threshold at a natural reference (REF_L2 from LOT_REF_SCAN, e.g. the box2 scan):
    eps * L_H = 0.1, 0.3, 1, 3 at REF_L2 -> S (L2_share_last5) falls and L2_strict is gone at 3."""
    out = {"REF": run({}, 901)}
    for x in (0.1, 0.3, 1.0, 3.0):
        out[f"epsL{x:g}"] = run({"eps": x / chem.L_H}, 901)
    keys = ["epsL0.1", "epsL0.3", "epsL1", "epsL3"]
    share = [mean(out[k]["share_last5"]) for k in keys]
    rate = [out[k]["L2_strict"] for k in keys]
    ok = share[0] > share[-1] and rate[-1] <= 1 and rate[0] > rate[-1]
    out["summary"] = {"eps_times_LH": [0.1, 0.3, 1, 3], "S_mean": share, "L2_strict": rate,
                      "REF_eps": REF_L2.get("eps"), "REF_L2_strict": out["REF"]["L2_strict"]}
    out["pass"] = bool(ok)
    return out


def main():
    if os.environ.get("LOT_KEYLEN"):                     # lineage keyed on a prefix (see chem4c3.L2_KEY_LEN)
        chem.L2_KEY_LEN = int(os.environ["LOT_KEYLEN"])
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "validate4c3.json")
    which = sys.argv[2] if len(sys.argv) > 2 else "12345"
    res = ckpt.load_json(path) or {}
    res["REF_L2"] = {k: v for k, v in REF_L2.items()}
    res["scan_pick"] = SCAN_PICK
    res["fidelity"] = FID
    res["reps"] = REPS
    for k, fn in (("1", test_1), ("2", test_2), ("3", test_3), ("4", test_4), ("5", test_5), ("7", test_7),
                  ("8", test_8), ("9", test_9)):
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
