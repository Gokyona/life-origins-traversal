"""Stage 4c-1: bottom-up polymer chemistry (levels E1 and E2), spec docs/stage4c1_spec.md.

Nothing living is put in. Each run is one environment (a point u in [0,1]^DIM); it has N
reaction sites (pores / droplets / pools) that share the environment and differ only by
chance. A site holds
  * monomer pools per chemistry c (0 = nucleotide-like, 1 = amino-acid-like): feedstock F_c,
    inactive monomers M_c and activated monomers M*_c (counts);
  * up to S polymers stored explicitly: length L (0 = empty slot) <= L_MAX, chemistry flag,
    and a bit-packed sequence over an alphabet of ALPHA letters (BITS bits per letter,
    monomer i at bits [i*BITS, (i+1)*BITS)). Letters are drawn uniformly at incorporation.
Chains only grow from monomers of their own chemistry.

Units: time in hours, amounts in molecules per site, V = 1000 molecules per site = 1 mM, so
"concentration" c = count / V in mM.

Processes per step (operator splitting, adaptive tau-leap, dt per run):
  1. linear monomer network per chemistry, sampled exactly over dt (4-state Markov chain
     F -> M <-> M* plus washout of every free species at 1/tau; feed arrivals are Poisson
     and spread over the step with the integrated transition matrix):
       feed      s_F * share_c                 share_0 = x_nuc, share_1 = 1 - x_nuc
       synthesis F -> M   k_syn = E,           E = Y_UV * UV + Y_H2 * H2
       activation M -> M* k_act = E * f_dry^d in the dry phase, E in the wet phase
       deactivation M* -> M  k_deact = 0.1/h * exp(-Ea/R (1/T - 1/298.15))
     then mixing: pools relax towards the mean over the run's sites at rate `mix`
     (monomers only; chains do not move between sites in this layer).
  2. chain growth: per chain end g = k_c * c_eff * (M*_c / V) * f_dry^d(dry phase), number of
     monomers added per step ~ Poisson(g dt) capped at G_CAP (and at L_MAX - L); the
     expected demand is thinned to the available M*_c.
  3. initiation M* + M* -> dimer: Poisson(k_c * c_eff * f_dry^d * (M*_c)^2 / (2V) * dt),
     needs a free slot.
  4. per chain: adsorbed to the mineral with probability theta (redrawn each step = fast
     exchange); washout 1/tau for chains that are not adsorbed; at most one bond break per
     step, rate per bond k_h * a_w * (1 - p_min if adsorbed) + k_UV * UV, at a uniform bond.
     Fragments of length 1 become inactive monomers M_c; a right fragment of length >= 2
     takes a free slot (if none is free its monomers return to M_c, counted as overflow).
c_eff = evaporation^d (dry phase only) * thermal-gradient factor * (1 + 9 theta).
Hydrolysis per bond (per hour):
  nucleotide-like: 1.4e-7/min (Li & Breaker 1999 reference: 37 C, pH 7.4, 5 mM Mg, 0.25 M K)
     * 10^(0.983 (pH - 7.4)) * 10^(0.07 (T - 37)) * (Mg / 0.005)^0.8 * (max(sal, 0.01) / 0.25)^s_exp
     (base catalysis only: acid catalysis below pH ~5 is not in the law; Mg outside 5-50 mM
     and salinity outside 0.03-3.16 M are extrapolations)
  amino-acid-like: Arrhenius line through t1/2 = 400 y at 25 C and 5.5 weeks at 100 C
     (Radzicka & Wolfenden 1996): Ea = 101.7 kJ/mol, no pH / ion dependence.
Cycling (square wave, period from layer A, depth d in [0, 1], first half wet, second dry):
  dry phase: a_w = 1 - d (1 - a_w,dry), concentration x evaporation^d, growth and activation
  x f_dry^d, temperature T + 10 d (wet phase T - 10 d, floored at 0 C). When the half period
  is shorter than the target step, the phases are averaged (rates time-averaged over a cycle;
  flagged `cycle_averaged`).
Liquid water: pressure must be >= saturation vapour pressure at the hottest phase temperature
  (Antoine equation for water, NIST/Perry constants, two ranges: 1-100 C and 99-374 C);
  otherwise E1 and E2 fail (flag `liquid`).
Time: real available time (1 - 1e6 years) cannot be stepped through. The monomer pools
  start at the expected state of the linear network (no chains) after feeding for the
  available time, so slow pool build-up is capped by the available time; pools are capped at
  POOL_MAX = 1e6 molecules (1 M; precipitation). A run then lasts
  t_end = min(5 * t_slow, available time), t_slow = slowest of the activated-monomer,
  chain-turnover and chain-growth (to length 8) time scales, where activation and growth
  count at most 4 chain lifetimes (beyond that chains are quasi-steady; see _time_end).
  Target T_STEPS steps per run (dt_target = t_end / T_STEPS); dt is further limited by bond
  breaking (P(break) <= 0.3 for a chain of L_MAX) and by phase switches; a run that needs more
  than MAX_STEP_FACTOR * T_STEPS steps stops
  early (flag `reached_end` False). `stationary`: chain mass in the last quarter of the run
  within 10% of the third quarter.
Detectors (E2: time averages over the last quarter of the run):
  E1: effective activated-monomer availability c_eff * A*/V reaches C_WORK (0.01 mM) at some
      time during the run ("reach a working concentration"); liquid water. The last-quarter
      average version is reported as E1_last_quarter.
      (No minimum count per site: the site volume is arbitrary, and thermophoretic or mineral
      concentration can make a sub-molecule mean count an effective concentration.)
      A* = free M* + (monomers incorporated into chains per hour) / (k_deact + 1/tau), the M*
      stock that would stand without polymerisation, so that fast polymerisation does not hide
      the supply of activated monomers.
  E2(L_min): fraction of sites holding a chain >= L_min is >= 0.5; liquid water. Primary
      L_min = 8, also 4 and 16. E2 is reported raw; nested maps use E1 & E2.
Layer A variable 9 (amphiphile-precursor fraction) is decoded but unused in this layer.
"""
import math
import numpy as np
import torch

A_NAMES = ["T", "P", "pH", "H2", "UV", "salinity", "Mg", "feed", "x_nuc", "amph", "tau", "mix",
           "evap", "thermo", "depth", "period", "mineral", "mineral_type", "avail_time"]
C_NAMES = ["k_c", "f_dry", "aw_dry", "Y_UV", "Y_H2", "Ea_deact", "theta", "p_min", "s_exp", "k_UV"]
NAMES = A_NAMES + C_NAMES
DIM = len(NAMES)          # 29
INACTIVE = ["amph"]

V = 1000.0                # molecules per site per mM
C_WORK = 0.01             # mM, working concentration of activated monomers (E1)
L_MAX = 32
ALPHA = 2
BITS = 1
G_CAP = 4                 # max monomers added to one chain in one step
L_MINS = (4, 8, 16)
K_DEACT_25 = 0.1          # per hour at 25 C
K_NUC_REF = 1.4e-7 * 60   # per hour, Li & Breaker reference
K_PEP_25 = math.log(2) / (400 * 8766.0)
K_PEP_100 = math.log(2) / (5.5 * 7 * 24.0)
EA_R_PEP = math.log(K_PEP_100 / K_PEP_25) / (1 / 298.15 - 1 / 373.15)   # K
FIDELITY = {"full": (256, 2048, 300), "screen": (64, 512, 300)}         # (sites, slots, T_STEPS)
MAX_STEP_FACTOR = 2
POOL_MAX = 1e6            # molecules per site per pool (1 M): beyond this, feedstock precipitates
BREAK_P = 0.3             # max probability of a bond break per step for a chain of L_MAX


def decode(u):
    """u [B, DIM] in [0,1] -> dict of float64 CPU tensors [B] (physical values)."""
    u = u.detach().double().cpu()
    lg = lambda i, lo, hi: 10 ** (lo + (hi - lo) * u[:, i])
    p = {
        "T": 150.0 * u[:, 0], "P": lg(1, 0, math.log10(1100)), "pH": 2 + 10 * u[:, 2],
        "H2": lg(3, -6, 0), "UV": lg(4, -4, 0), "salinity": 1.2 * u[:, 5], "Mg": lg(6, -4, -1),
        "feed": lg(7, -3, 1), "x_nuc": u[:, 8], "amph": 0.5 * u[:, 9], "tau": lg(10, 0, 5),
        "mix": lg(11, -4, 1), "evap": lg(12, 0, 3), "thermo": lg(13, 0, 8), "depth": u[:, 14],
        "period": lg(15, 0, 3), "mineral": lg(16, -3, 0),
        "mineral_type": torch.clamp(torch.floor(3 * u[:, 17]), max=2),
        "avail_time": lg(18, 0, 6) * 8766.0,
        "k_c": lg(19, -2, 2), "f_dry": lg(20, 0, 3), "aw_dry": lg(21, -3, math.log10(0.3)),
        "Y_UV": lg(22, -4, 0), "Y_H2": lg(23, -4, 0), "Ea_deact": 40 + 60 * u[:, 24],
        "theta": lg(25, -2, math.log10(0.9)), "p_min": 0.99 * u[:, 26], "s_exp": -1 + 2 * u[:, 27],
        "k_UV": lg(28, -6, -2),
    }
    return p


def p_sat_bar(T):
    """Saturation vapour pressure of water (bar), Antoine equation (mmHg constants)."""
    lo = 8.07131 - 1730.63 / (233.426 + T)
    hi = 8.14019 - 1810.94 / (244.485 + T)
    return 10 ** torch.where(T <= 100, lo, hi) * 133.322 / 1e5


def k_nuc(T, p):
    return (K_NUC_REF * 10 ** (0.983 * (p["pH"] - 7.4)) * 10 ** (0.07 * (T - 37))
            * (p["Mg"] / 0.005) ** 0.8 * (torch.clamp(p["salinity"], min=0.01) / 0.25) ** p["s_exp"])


def k_pep(T):
    return K_PEP_25 * torch.exp(-EA_R_PEP * (1 / (T + 273.15) - 1 / 298.15))


def phase_table(p):
    """Rates per phase state: index 0 wet, 1 dry, 2 cycle average. Each value [B, 3]."""
    d = p["depth"]
    T_ph = torch.stack([torch.clamp(p["T"] - 10 * d, min=0.0), p["T"] + 10 * d], 1)
    aw = torch.stack([torch.ones_like(d), 1 - d * (1 - p["aw_dry"])], 1)
    cev = torch.stack([torch.ones_like(d), p["evap"] ** d], 1)
    fd = torch.stack([torch.ones_like(d), p["f_dry"] ** d], 1)
    E = p["Y_UV"] * p["UV"] + p["Y_H2"] * p["H2"]
    theta = p["theta"] * p["mineral"] * torch.tensor([0.0, 1.0, 0.5], dtype=torch.float64)[p["mineral_type"].long()]
    cmin = 1 + 9 * theta
    t = {
        "kact": E[:, None] * fd,
        "kdeact": K_DEACT_25 * torch.exp(-p["Ea_deact"][:, None] * 1000 / 8.314
                                         * (1 / (T_ph + 273.15) - 1 / 298.15)),
        "knuc": torch.stack([k_nuc(T_ph[:, i], p) for i in range(2)], 1) * aw,
        "kpep": torch.stack([k_pep(T_ph[:, i]) for i in range(2)], 1) * aw,
        "ceff": cev * (p["thermo"] * cmin)[:, None],
        "gmult": fd,
    }
    for k in list(t):
        t[k] = torch.cat([t[k], t[k].mean(1, keepdim=True)], 1)
    return t, E, theta, T_ph


def _linear_q(E, kact, kdeact, w):
    """Generator of the linear monomer network, states F, M, M*, out. Inputs [B] float64."""
    B = E.shape[0]
    Q = torch.zeros(B, 4, 4, dtype=torch.float64)
    Q[:, 0, 1] = E
    Q[:, 1, 2] = kact
    Q[:, 2, 1] = kdeact
    Q[:, 0, 3] = w
    Q[:, 1, 3] = w
    Q[:, 2, 3] = w
    for i in range(3):
        Q[:, i, i] = -Q[:, i].sum(1)
    return Q


def _initial_pools(p, t, E, share):
    """Expected pools [B,2,3] (F, M, M*) after feeding for the available time with the
    cycle-averaged linear network and no chains, capped at POOL_MAX per pool."""
    Q = _linear_q(E, t["kact"][:, 2], t["kdeact"][:, 2], 1 / p["tau"])
    T = p["avail_time"]
    _, J = _expm_pair(Q, T)
    n = (p["feed"] * T)[:, None] * J[:, 0, :3]                       # [B,3] per unit share
    n = torch.clamp(n, max=POOL_MAX)
    return share[:, :, None] * n[:, None, :]


def _time_end(p, t, E, pools, n_steps):
    """Run length: 5 x the slowest of the activated-monomer, chain-turnover and chain-growth
    time scales, capped by the available time. Slow pool build-up is not stepped through:
    pools start at their state after the available time (_initial_pools)."""
    w = 1 / p["tau"]
    kact = t["kact"][:, 2]
    kdeact = t["kdeact"][:, 2]
    kbr = torch.maximum(t["knuc"][:, 2], t["kpep"][:, 2]) + p["k_UV"] * p["UV"]
    ms = pools[:, :, 2].sum(1)
    g = p["k_c"] * t["ceff"][:, 2] * t["gmult"][:, 2] * ms / V
    t_grow = 8 / torch.clamp(g, min=1e-300)
    t_ch = 1 / (4 * kbr + w)
    # activation and growth slower than 4 chain lifetimes cannot change whether chains persist
    # within the run (chains are quasi-steady); capping them keeps break-limited runs short
    t_slow = torch.maximum(t_ch, torch.minimum(torch.maximum(1 / (kact + kdeact + w), t_grow), 4 * t_ch))
    t_end = torch.minimum(5 * t_slow, p["avail_time"])
    return t_end, t_end / n_steps


def _expm_pair(Q, dt):
    """P = exp(Q dt) and J = int_0^dt exp(Q s) ds / dt for Q [B,4,4], dt [B] (float64 CPU)."""
    B = Q.shape[0]
    Z = torch.zeros(B, 8, 8, dtype=torch.float64)
    Z[:, :4, :4] = Q * dt[:, None, None]
    Z[:, :4, 4:] = torch.eye(4, dtype=torch.float64) * dt[:, None, None]
    E = torch.linalg.matrix_exp(Z)
    P = E[:, :4, :4]
    J = E[:, :4, 4:] / torch.clamp(dt, min=1e-300)[:, None, None]
    J = torch.where(dt[:, None, None] > 0, J, torch.eye(4, dtype=torch.float64).expand(B, 4, 4))
    return P, J


def _multinomial(n, probs, gen):
    """n [...] counts, probs [..., K] -> counts [..., K] by sequential binomials."""
    out = []
    rest = n
    left = torch.ones_like(n)
    for k in range(probs.shape[-1] - 1):
        q = torch.nan_to_num(torch.clamp(probs[..., k] / torch.clamp(left, min=1e-12), 0.0, 1.0), nan=0.0)
        x = torch.binomial(rest, q.contiguous(), generator=gen)
        out.append(x)
        rest = rest - x
        left = left - probs[..., k]
    out.append(rest)
    return torch.stack(out, -1)


def _rate(x):
    """Poisson rates must be finite and >= 0 (the XPU kernel aborts the process otherwise)."""
    return torch.nan_to_num(x, nan=0.0, posinf=1e30).clamp(min=0.0).contiguous()


def _mask(nbits):
    one = torch.ones_like(nbits)
    return torch.where(nbits >= 63, torch.full_like(nbits, -1), torch.bitwise_left_shift(one, nbits) - 1)


@torch.no_grad()
def simulate(u, generator, fidelity="screen", record=False, override=None, init=None):
    """Simulate B environments. Returns dict of per-run results (CPU tensors / numpy).

    override: dict of physical values replacing decoded ones (tensors [B] or scalars), plus
      special keys for validation only: "clamp_mstar" (hold M* per site at this count, split by
      x_nuc share), "t_end" (fixed run length in hours), "no_cycle_average" (always resolve
      phases), "empty_start" (pools start at zero), "k_h_fixed" (per-bond hydrolysis rate per
      hour, all phases), "hydro_mult" (multiply hydrolysis), "init_factor" (multiply initiation).
    init: {"dimers": n} start every site with n dimers of the dominant chemistry (validation).
    """
    dev = u.device
    gen = generator
    N, S, n_steps = FIDELITY[fidelity] if isinstance(fidelity, str) else fidelity
    B = u.shape[0]
    R = B * N
    p = decode(u)
    override = dict(override or {})
    special = {k: override.pop(k) for k in ("clamp_mstar", "t_end", "no_cycle_average", "empty_start",
                                            "k_h_fixed", "hydro_mult", "init_factor") if k in override}
    for k, v in override.items():
        p[k] = torch.as_tensor(v, dtype=torch.float64).expand(B).clone()
    tab, E, theta, T_ph = phase_table(p)
    if "k_h_fixed" in special:
        kf = torch.as_tensor(special["k_h_fixed"], dtype=torch.float64).expand(B)[:, None]
        tab["knuc"] = kf.expand(B, 3).clone()
        tab["kpep"] = kf.expand(B, 3).clone()
    if "hydro_mult" in special:
        hm = torch.as_tensor(special["hydro_mult"], dtype=torch.float64).expand(B)[:, None]
        tab["knuc"] = tab["knuc"] * hm
        tab["kpep"] = tab["kpep"] * hm
    init_f = torch.as_tensor(special.get("init_factor", 1.0), dtype=torch.float64).expand(B)
    share = torch.stack([p["x_nuc"], 1 - p["x_nuc"]], 1)            # [B,2]
    pools0 = _initial_pools(p, tab, E, share) if not special.get("empty_start") else torch.zeros(B, 2, 3, dtype=torch.float64)
    t_end, dt_target = _time_end(p, tab, E, pools0, n_steps)
    if "t_end" in special:
        t_end = torch.as_tensor(special["t_end"], dtype=torch.float64).expand(B).clone()
        dt_target = t_end / n_steps
    half = p["period"] / 2
    cycling = p["depth"] > 0
    resolved = cycling & ((half >= dt_target) | bool(special.get("no_cycle_average", False)))
    liquid = p["P"] >= p_sat_bar(T_ph.max(1).values)
    w = 1 / p["tau"]
    kuv = p["k_UV"] * p["UV"]

    f32 = lambda x: x.to(dev, torch.float32)
    run_of_row = torch.arange(R, device=dev) // N
    # state
    lam0 = f32(pools0)[run_of_row]                                       # [R,2,3]
    F = torch.poisson(_rate(lam0[:, :, 0]), generator=gen)
    M = torch.poisson(_rate(lam0[:, :, 1]), generator=gen)
    Ms = torch.poisson(_rate(lam0[:, :, 2]), generator=gen)
    L = torch.zeros(R, S, dtype=torch.int64, device=dev)
    seq = torch.zeros(R, S, dtype=torch.int64, device=dev)
    chem = torch.zeros(R, S, dtype=torch.int64, device=dev)
    if init and init.get("dimers"):
        nd = int(init["dimers"])
        L[:, :nd] = 2
        seq[:, :nd] = torch.randint(0, 1 << (2 * BITS), (R, nd), device=dev, generator=gen)
        chem[:, :nd] = (share[:, 1] > share[:, 0]).long().to(dev)[run_of_row][:, None]
    clamp = special.get("clamp_mstar")
    if clamp is not None:
        clamp = f32(share * torch.as_tensor(clamp, dtype=torch.float64).expand(B)[:, None])[run_of_row].round()
        Ms = clamp.clone()

    t = torch.zeros(B, dtype=torch.float64)
    win_lo = 0.75 * t_end
    acc = {k: torch.zeros(B, dtype=torch.float64) for k in ("time", "conc", "ms", "q3mass", "q4mass", "q3t")}
    acc_frac = torch.zeros(B, len(L_MINS), dtype=torch.float64)
    t_first_e1 = torch.full((B,), float("nan"), dtype=torch.float64)
    t_first_e2 = torch.full((B,), float("nan"), dtype=torch.float64)
    ms_series = torch.full((B, 33), float("nan"), dtype=torch.float64)
    next_rec = torch.zeros(B, dtype=torch.long)
    overflow = torch.zeros(B, dtype=torch.float64)
    overdraw = torch.zeros(B, dtype=torch.float64)
    steps = 0
    max_steps = MAX_STEP_FACTOR * n_steps

    while steps < max_steps:
        active = t < t_end * (1 - 1e-9)
        if not bool(active.any()):
            break
        # phase state per run
        k_half = torch.floor(t / half * (1 + 1e-12) + 1e-9)     # robust at phase boundaries
        ph = torch.where(resolved, (k_half % 2).long(), torch.where(cycling, 2, 0))
        to_switch = torch.where(resolved, (k_half + 1) * half - t, torch.full_like(t, float("inf")))
        g = lambda name: tab[name].gather(1, ph[:, None]).squeeze(1)
        kact, kdeact, knuc, kpep, ceff, gmult = (g(n) for n in ("kact", "kdeact", "knuc", "kpep", "ceff", "gmult"))
        # break limit from the chemistries that are fed (a chemistry with zero share has no chains)
        kbr_max = torch.maximum(torch.where(share[:, 0] > 0, knuc, 0.0), torch.where(share[:, 1] > 0, kpep, 0.0)) + kuv
        dt = torch.minimum(dt_target, BREAK_P / ((L_MAX - 1) * kbr_max))
        dt = torch.minimum(dt, to_switch)
        dt = torch.minimum(dt, t_end - t)
        dt = torch.where(active, torch.clamp(dt, min=0.0), torch.zeros_like(dt))
        dtr = f32(dt)[run_of_row]

        # 1. linear monomer network (exact over dt), states F, M, M*, out
        P, J = _expm_pair(_linear_q(E, kact, kdeact, w), dt)
        Pd, Jd = f32(P)[run_of_row], f32(J)[run_of_row]                  # [R,4,4]
        arrivals = torch.poisson(_rate(f32(p["feed"][:, None] * share * dt[:, None])[run_of_row]), generator=gen)
        newF = torch.zeros_like(F); newM = torch.zeros_like(M); newMs = torch.zeros_like(Ms)
        for c in range(2):
            tot = torch.zeros(R, 4, device=dev)
            for src, n in ((0, F[:, c]), (1, M[:, c]), (2, Ms[:, c])):
                tot += _multinomial(n, Pd[:, src, :], gen)
            tot += _multinomial(arrivals[:, c], Jd[:, 0, :], gen)
            newF[:, c], newM[:, c], newMs[:, c] = tot[:, 0], tot[:, 1], tot[:, 2]
        F, M, Ms = newF, newM, newMs
        mixf = f32(1 - torch.exp(-p["mix"] * dt))[run_of_row][:, None]
        for X in (F, M, Ms):
            mean = X.view(B, N, 2).mean(1, keepdim=True).expand(B, N, 2).reshape(R, 2)
            X += (mean - X) * mixf
            X.copy_(torch.clamp(torch.floor(X + torch.rand(X.shape, device=dev, generator=gen)), max=POOL_MAX))
        if clamp is not None:
            Ms = clamp.clone()

        # 2. growth
        ceffr = f32(ceff)[run_of_row][:, None]
        gm = f32(gmult)[run_of_row][:, None]
        kc = f32(p["k_c"])[run_of_row][:, None]
        occ = L > 0
        growable = occ & (L < L_MAX)
        n_grow = torch.stack([(growable & (chem == c)).sum(1) for c in range(2)], 1).float()
        lam = kc * ceffr * gm * (Ms / V) * dtr[:, None]
        demand = lam * n_grow
        lam = lam * torch.clamp(Ms / torch.clamp(demand, min=1e-9), max=1.0)
        lam_c = torch.gather(lam, 1, chem)                                    # [R,S]
        uu = torch.rand(R, S, device=dev, generator=gen)
        e = torch.exp(-lam_c)
        cdf = e.clone()
        term = e.clone()
        k_add = torch.zeros_like(L)
        for k in range(1, G_CAP + 1):
            k_add += (uu > cdf).long()
            term = term * lam_c / k
            cdf = cdf + term
        k_add = torch.where(growable, torch.minimum(k_add, L_MAX - L), torch.zeros_like(k_add))
        bits = torch.randint(0, 1 << (G_CAP * BITS), (R, S), device=dev, generator=gen)
        seq = seq | torch.bitwise_left_shift(bits & _mask(k_add * BITS), L * BITS)
        L = L + k_add
        used = torch.stack([(k_add * (chem == c)).sum(1) for c in range(2)], 1).float()
        inc_site = used.sum(1).view(B, N).mean(1).double().cpu()        # incorporated this step
        Ms = Ms - used
        overdraw += torch.clamp(-Ms, min=0).view(B, N, 2).sum((1, 2)).double().cpu()
        Ms = torch.clamp(Ms, min=0)

        # 4. per chain: adsorption, washout, one bond break
        th = f32(theta)[run_of_row][:, None]
        ads = torch.rand(R, S, device=dev, generator=gen) < th
        kb = torch.where(chem == 0, f32(knuc)[run_of_row][:, None], f32(kpep)[run_of_row][:, None])
        kb = torch.where(ads, kb * f32(1 - p["p_min"])[run_of_row][:, None], kb) + f32(kuv)[run_of_row][:, None]
        occ = L > 0
        pb = 1 - torch.exp(-(L - 1).clamp(min=0).float() * kb * dtr[:, None])
        pw = torch.where(ads, torch.zeros_like(kb), (1 - torch.exp(-f32(w)[run_of_row] * dtr))[:, None].expand(R, S))
        u2 = torch.rand(R, S, device=dev, generator=gen)
        wash = occ & (u2 < pw)
        v = (u2 - pw) / torch.clamp(1 - pw, min=1e-12)
        brk = occ & ~wash & (v < pb) & (L >= 2)
        j = 1 + torch.floor(torch.clamp(v / torch.clamp(pb, min=1e-30), max=0.999999) * (L - 1).float()).long()
        j = torch.where(brk, torch.minimum(j, L - 1), torch.zeros_like(j))
        L = torch.where(wash, torch.zeros_like(L), L)
        right_len = torch.where(brk, L - j, torch.zeros_like(L))
        right_seq = torch.bitwise_right_shift(seq, j * BITS) & _mask(right_len * BITS)
        seq = torch.where(brk, seq & _mask(j * BITS), seq)
        L = torch.where(brk, j, L)
        mono = (brk & (j == 1)).long() + (brk & (right_len == 1)).long()
        L = torch.where(brk & (j == 1), torch.zeros_like(L), L)
        M = M + torch.stack([(mono * (chem == c)).sum(1) for c in range(2)], 1).float()
        seq = torch.where(L == 0, torch.zeros_like(seq), seq)

        # 3 + allocation: right fragments (>= 2) and initiations take free slots
        need = brk & (right_len >= 2)
        free = L == 0
        fc = free.sum(1)
        n_right = need.sum(1)
        rate_init = kc * ceffr * gm * (Ms / V) * Ms / 2 * dtr[:, None] * f32(init_f)[run_of_row][:, None]
        n_init = torch.poisson(_rate(rate_init), generator=gen)
        n_init = torch.minimum(n_init, torch.floor(Ms / 2))
        room = torch.clamp(fc - n_right, min=0).float()
        tot_init = n_init.sum(1)
        scale = torch.clamp(room / torch.clamp(tot_init, min=1), max=1.0)
        n_init = torch.floor(n_init * scale[:, None])
        Ms = Ms - 2 * n_init
        inc_site = inc_site + 2 * n_init.sum(1).view(B, N).mean(1).double().cpu()
        free_pos = free.flatten().nonzero().squeeze(1)
        row_off = torch.cumsum(fc, 0) - fc
        rr, cc = need.nonzero(as_tuple=True)
        if rr.numel():
            rank = (torch.cumsum(need.long(), 1) - 1)[rr, cc]
            ok = rank < fc[rr]
            lost_len = right_len[rr[~ok], cc[~ok]]
            lost_chem = chem[rr[~ok], cc[~ok]]
            if lost_len.numel():
                M.index_put_((rr[~ok], lost_chem), lost_len.float(), accumulate=True)
                overflow += torch.zeros(B, dtype=torch.float64).index_add_(
                    0, run_of_row[rr[~ok]].cpu(), torch.ones(int((~ok).sum()), dtype=torch.float64))
            rr, cc, rank = rr[ok], cc[ok], rank[ok]
            tgt = free_pos[row_off[rr] + rank]
            Lf, sf, cf = L.view(-1), seq.view(-1), chem.view(-1)
            src_len, src_seq, src_chem = right_len[rr, cc], right_seq[rr, cc], chem[rr, cc]
            Lf[tgt] = src_len
            sf[tgt] = src_seq
            cf[tgt] = src_chem
        n_right_ok = torch.minimum(n_right, fc)
        ni_tot = n_init.sum(1).long()
        if int(ni_tot.sum()):
            rows = torch.repeat_interleave(torch.arange(R, device=dev), ni_tot)
            start = torch.cumsum(ni_tot, 0) - ni_tot
            within = torch.arange(rows.numel(), device=dev) - start[rows]
            rank = n_right_ok[rows] + within
            tgt = free_pos[row_off[rows] + rank]
            chem_new = (within >= n_init[rows, 0].long()).long()
            Lf, sf, cf = L.view(-1), seq.view(-1), chem.view(-1)
            Lf[tgt] = 2
            sf[tgt] = torch.randint(0, 1 << (2 * BITS), (rows.numel(),), device=dev, generator=gen)
            cf[tgt] = chem_new
        if clamp is not None:
            Ms = clamp.clone()

        # detectors and records
        t_new = t + dt
        # activated-monomer availability: free M* plus what polymerisation took up, i.e. the
        # stock it would have formed without consumption (incorporation flux / loss rate)
        flux = inc_site / torch.clamp(dt, min=1e-300)
        msite = Ms.sum(1).view(B, N).mean(1).double().cpu() + torch.where(dt > 0, flux / (kdeact + w), 0.0)
        conc = ceff * msite / V
        has = torch.stack([(L >= lm).any(1) for lm in L_MINS], 1).float().view(B, N, -1).mean(1).double().cpu()
        mass = L.sum(1).view(B, N).sum(1).double().cpu()
        wq4 = torch.clamp(t_new - torch.maximum(t, win_lo), min=0)
        wq3 = torch.clamp(torch.minimum(t_new, win_lo) - torch.maximum(t, 0.5 * t_end), min=0)
        acc["time"] += wq4
        acc["conc"] += wq4 * conc
        acc["ms"] += wq4 * msite
        acc["q4mass"] += wq4 * mass
        acc["q3mass"] += wq3 * mass
        acc["q3t"] += wq3
        acc_frac += wq4[:, None] * has
        e1_now = conc >= C_WORK
        e2_now = has[:, L_MINS.index(8)] >= 0.5
        t_first_e1 = torch.where(torch.isnan(t_first_e1) & e1_now & active, t_new, t_first_e1)
        t_first_e2 = torch.where(torch.isnan(t_first_e2) & e2_now & active, t_new, t_first_e2)
        while True:
            grid_t = t_end * next_rec.double() / 32
            due = (next_rec <= 32) & (t_new >= grid_t * (1 - 1e-12)) & (active | (next_rec == 0))
            if not bool(due.any()):
                break
            idx = due.nonzero().squeeze(1)
            ms_series[idx, next_rec[idx]] = msite[idx]
            next_rec[idx] += 1
        t = t_new
        steps += 1

    tw = torch.clamp(acc["time"], min=1e-300)
    conc_avg = acc["conc"] / tw
    ms_avg = acc["ms"] / tw
    frac_avg = acc_frac / tw[:, None]
    q4 = acc["q4mass"] / tw
    q3 = acc["q3mass"] / torch.clamp(acc["q3t"], min=1e-300)
    stationary = (q4 - q3).abs() <= 0.1 * torch.clamp(torch.maximum(q4, q3), min=1.0)
    e1_avg = liquid & (conc_avg >= C_WORK)
    e1 = liquid & ~torch.isnan(t_first_e1)          # reached at any time during the run
    out = {"E1": e1, "E1_last_quarter": e1_avg, "liquid": liquid, "reached_end": t >= t_end * (1 - 1e-9), "t_reached": t.clone(),
           "stationary": stationary, "cycle_averaged": cycling & ~resolved,
           "conc_mstar": conc_avg, "mstar_per_site": ms_avg, "t_end": t_end, "steps": torch.full((B,), steps),
           "overflow_events": overflow, "overdraw_molecules": overdraw}
    for i, lm in enumerate(L_MINS):
        out[f"frac_sites_L{lm}"] = frac_avg[:, i]
        out[f"E2_L{lm}"] = liquid & (frac_avg[:, i] >= 0.5)
    out["E2"] = out["E2_L8"]
    if record:
        out.update(_composition(L, seq, chem, F, M, Ms, B, N, R, run_of_row))
        out["t_first_E1"] = t_first_e1
        out["t_first_E2"] = t_first_e2
        out["mstar_series"] = ms_series
        out["theta_adsorbed"] = theta
    return out


def _composition(L, seq, chem, F, M, Ms, B, N, R, run_of_row):
    dev = L.device
    occ = L > 0
    runs = run_of_row[:, None].expand(R, L.shape[1])
    hist = torch.zeros(B, L_MAX + 1, device=dev)
    hist.index_put_((runs[occ], L[occ]), torch.ones(int(occ.sum()), device=dev), accumulate=True)
    mass_c = torch.zeros(B, 2, device=dev)
    mass_c.index_put_((runs[occ], chem[occ]), L[occ].float(), accumulate=True)
    pools = torch.stack([F, M, Ms], 1).view(B, N, 3, 2).sum(1)                # [B,3,2]
    classes = [(2, 3), (4, 7), (8, 15), (16, 32)]
    distinct = torch.zeros(B, len(classes))
    entropy = torch.zeros(B, len(classes))
    for ci, (lo, hi) in enumerate(classes):
        m = occ & (L >= lo) & (L <= hi)
        if not bool(m.any()):
            continue
        key = (runs[m] << 40) | (chem[m] << 39) | (L[m] << 33) | (seq[m] & ((1 << 33) - 1))
        uk, cnt = torch.unique(key, return_counts=True)
        r = (uk >> 40).cpu()                      # sums on CPU in float64: deterministic order
        cnt = cnt.cpu().double()
        distinct[:, ci] = torch.zeros(B, dtype=torch.float64).index_add_(0, r, torch.ones_like(cnt)).float()
        tot = torch.zeros(B, dtype=torch.float64).index_add_(0, r, cnt)
        pr = cnt / tot[r]
        entropy[:, ci] = torch.zeros(B, dtype=torch.float64).index_add_(0, r, -pr * torch.log2(pr)).float()
    return {"length_hist": hist.cpu(), "chain_mass_by_chem": mass_c.cpu(), "pools_F_M_Mstar_by_chem": pools.cpu(),
            "distinct_seq_by_class": distinct, "entropy_bits_by_class": entropy,
            "length_classes": classes}
