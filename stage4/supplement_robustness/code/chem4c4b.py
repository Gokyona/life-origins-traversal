"""P3-1b (Review 4, R4.3): copy of lot-s4/chem4c4.py (md5 b79360f6, not modified) with the P3-1b switches ported from
chem_p31b.py (md5 a31cae30 / 60c68e27, the 4c-3 branch): MULTI_START (renewal rule), SLOT_ROUND_FIX, POISSON_INV, all
default False (False = chem4c4 bit for bit). Compartment specifics: primers are counted per (site, compartment), a template
inside a vesicle extends at r * F_PERM, a copy inherits its template's compartment.
--- chem4c4 description follows ---
Stage 4c-4 (chem4c4 v2): amphiphiles and vesicles (level E3), spec docs/stage4c4_spec.md (md5 75191e2a).
Copy of chem4c3.py (md5 cab66933) with these changes:
  * twelve layer-C dimensions appended (54 in all): CVC0 (0.01-100 mM, log), gamma_T (0-0.05 /C),
    pKa_app (6-10), w_pH (0.3-2), Mg_crit (0.5-50 mM, log), sal_crit (0.03-3 M, log), T_mem (40-100 C),
    phi0 (1e-3-0.3, log), f_rehyd (1-10), beta (0-10), s_shear (0.01-10, log), s_cycle (0.01-10, log);
    k_ex = 10^-1.5 /h and f_perm = 0.1 fixed (geometric middles of 1e-3..1 and 0.01..1);
  * the amphiphile-precursor fraction amph (layer A variable 9) takes its share of the feed (the
    monomer chemistries share the rest); precursor -> free amphiphile A at the synthesis rate E;
    A, precursor and non-vesicle membrane wash out at 1/tau; A mixes between sites;
  * assembly: free A above CVC = CVC0 exp(gamma_T (T - 25)) (concentration x the evaporation factor
    in the dry phase) goes into membrane at k_as = 1e3 /h, membrane returns below it (first the
    non-vesicle membrane, then vesicle membrane in proportion; a vesicle below a0/2 dissolves);
  * stability window per run: |pH - pKa_app| <= w_pH, Mg <= Mg_crit, salinity <= sal_crit, hottest
    phase <= T_mem; outside it vesicles break at k_break = 1 /h (see SPEC_NOTES) and release contents;
  * up to NV = 8 vesicles per site (membrane a_j, alive, id, parent, generation, birth time) and a
    compartment label per chain slot (int8: 0 outside, j + 1 inside vesicle j); a0 = 100 amphiphile
    molecules (newborn membrane), enclosed volume fraction phi_j = phi0 (a_j / a0)^1.5 (<= 1);
  * formation (window): Poisson(k_form B_free / a0 dt), k_form = 10 /h, at most the free vesicle
    slots and floor(B_free / a0); each outside chain enters new vesicle k with probability phi0
    (disjoint intervals of the same uniform draw); rehydration at every dry -> wet switch: outside
    chains enter vesicle j with probability min(phi_j f_rehyd, 1) (sequentially over j); v2: in
    cycle-averaged steps with wet-dry cycling, as an averaged effect of n = dt / period cycles,
    probability 1 - (1 - min(phi_j f_rehyd, 1))^n;
  * growth: non-vesicle membrane is taken up at k_grow = 10 /h, split in proportion to
    a_j (1 + beta Pi_j); exchange between vesicles relaxes a_j towards the same shares at k_ex;
    Pi_j = (enclosed letters / phi_j) / (outside letters / outside volume + 1);
  * division: a_j >= 2 a0 divides at k_div = s_shear mix + s_cycle (wet-dry switches per hour);
    membrane halves, each enclosed chain goes to the new daughter with probability 1/2, paired
    strands go together (one draw per duplex); without a free slot the second daughter's
    membrane returns to B_free and its contents are released;
  * losses: vesicles wash out at 1/tau with their contents (enclosed chains are not washed out
    individually); break down outside the window (contents released, membrane to B_free);
  * chemistry inside: enclosed chains grow from M* x f_perm; template copying extends at
    r x f_perm inside a vesicle; primers must share the template's compartment; a copy inherits
    the template's label; catalysed ligation only among species of one compartment (cat4c4);
    initiation makes outside dimers; a right fragment of length >= 4 stays in its vesicle;
  * runs keep chem4c3's horizon, step count and dt (no run extension for vesicles: coordinator
    decision 28 Oct 2026, after the pre-flight overflowed with a 10-turnover extension); E3 is
    evaluated over that horizon and flagged E3_evaluable = False when the last 8 checkpoints span
    fewer than 5 vesicle turnover times (1 / (1/tau + k_break [outside window]));
  * E3 detector at the 16 lineage checkpoints (Vesicles.outputs; v2: window parallel to L2, the
    longest run of checkpoints ending at t_end with >= half the sites holding a vesicle at each);
    vesicle composition fields;
  * hard guard: if no run of the batch has amph > 0 the vesicle layer is skipped entirely and
    draws no random numbers (bit-identical to chem4c3); the vesicle layer draws only from its own
    generator (seed + 104729).
--- chem4c3 description follows ---
Stage 4c-3: template copying and heredity (level L2) on top of the 4c-2 chemistry (spec
docs/stage4c3_spec.md). Copy of chem4c2.py (md5 f5eefabc) with these changes:
  * five layer-C dimensions appended (42 in all): L_t (4-12), k_t (1e1..1e5, log), eps
    (1e-3..0.3, log), dH_bp (25-45 kJ/mol per pair), T_m8 (10-90 C);
  * the detector settings adopted by sens4c2 are the defaults: NS_RAF = 8, REFRESH = 1,
    CAT_AFTER_E1 = True;
  * template copying (step 2c, nucleotide-like chains only): a free chain of chemistry 0 with
    L >= L_t starts at most one copy per step, with probability 1 - exp(-lam dt). Two kinds of
    start (coordinator decision, 26 Oct 2026): unprimed, one activated monomer, rate
    r = k_c f_min k_t Arr_cat(T) g_dry c_eff (M*_0 / V) per hour (the same r is the extension
    rate per position), at the template's first letter (p = 0); and primed, a free chain of chemistry 0
    and length l_p in PRIMER_LENS (default 2 and 3; never a template itself) from the site that
    is the complement of a template segment, rate
    k_on n_match(l_p) with k_on = r / M*_0 and n_match the number of free primers that are the
    complement of the template's first l_p letters (k_on n_primer 2^-l_p for a uniform primer
    pool). Copies start at the template's first letter (spec change approved 26 Oct 2026; the
    draft spec's random start position made full-length copies a 1/L event, so no lineage could
    re-make itself: tests 4 and 7 failed for that reason). The primer becomes the start of
    the copy (in its own slot) and extension continues from l = l_p. lam is the sum of the
    rates; a start that finds no unused matching primer (two templates in one step) fails.
    Finite copy time (COPY_TIME, coordinator decision 26 Oct 2026): after starting a copy that
    adds n letters, a template cannot start another copy for n / r hours;
    Position l+1
    is added with probability q(l) = r / (r + k_off(l)), up to the template end. Each
    added letter is the complement of the template letter (0 <-> 1, same bit order as the
    template segment), wrong with probability eps. A copy of >= 2 letters (an unprimed one in a
    free slot, a primed one in its primer's slot) uses its added letters' activated monomers of
    chemistry 0 and stays paired with its template (pair[], dlen[]).
    Copies are made after bond breaks and washout (step 3); copies beyond the site's M*_0 are
    dropped, and copies and initiations share the free slots left after right fragments in
    proportion to their rates (copying: r dt per free template times the released share of the
    row's starts, not capped at one start per template and step; initiation: its Poisson
    count); the copy share is rounded stochastically. Without this, initiation filled every
    freed slot first and copying was starved in every condition tried. The whole copy is made
    within one step (processive copying faster than a step: stated simplification);
  * duplexes: k_off(l, T) = K0 exp(dG/RT), dG = (l - 1)(dH_bp - T dS_bp), dS_bp = dH_bp /
    T_m8 (so dG(T_m8) = 0 for every length), K0 = 1e6 /h; a duplex of length dlen dissociates
    in a step with probability 1 - exp(-k_off dt). In a cycle-averaged phase k_off is the mean
    of the wet and dry phase values. Paired strands cannot template, be ligated or catalyse
    (they are removed from the abundant-species lists of cat4c2); they grow by monomer addition
    as before; a bond break in either strand dissolves the duplex, and washout takes both;
  * special override "copy_random": copies get random letters (copying on, information
    transfer off; the paired control of L2_strict); "copy_log": return every copy made
    (validation test 1); "init_templates": {"n", "L", "seq"} start every site with n copies of
    one nucleotide-like template;
  * runs with copying on (k_t > 0) last at least 10 dimer lifetimes, as runs with catalysis on;
  * L2 detector: at 16 checkpoints (t_end * k/16, all runs, not gated by E1) the abundant
    nucleotide-like sequences of length >= L_H (>= 2 copies, top NX_L2) of every site are
    stored with their complements. At the end, for each checkpoint k the information
    share is the share of the site's chain mass (chemistry 0, L >= L_H) whose sequence or
    complement was stored at k. Window length in turnover times of an L_H-mer (1 / ((L_H - 1)
    (k_h (1 - theta p_min) + k_UV) + (1 - theta)/tau), theta the adsorbed fraction: adsorption
    protects from hydrolysis and washout; the L1 window uses the same turnover formula, unlike
    4c-2; TURNOVER_PROTECT = False restores the 4c-2 formula). L2_loose: E1, liquid water and a window with share >= 0.2 spanning
    >= 5 turnover times (longest window, best site). L2_strict: the L2_loose runs are rerun as a
    batch of their own, as they are and twice with copy_random (one with the same random
    numbers, one with the next seed); the rerun must have a lasting lineage and its statistic S
    (largest share over all sites at the latest checkpoint whose window spans >= 5 turnover
    times, the same in all three) must exceed the larger control S by >= L2_MARGIN
    (paired_control_l2). Lineages are recorded in every site (all N, not only the NS_RAF
    sampled sites of the L1 detector). L2 is linked to E1 only;
  * composition fields: copies made, mean copy length, measured error rate, share of templates
    paired at the end, copies in the last quarter, L2 span and share, length of the most
    abundant lineage, first L2 time and its order relative to first L1, share of primed copies.
  * cat4c3.py (cat4c2 with deterministic tie-breaking in abundant()) replaces cat4c2.py, and
    torch.use_deterministic_algorithms(True) is set, so a run gives the same result in every
    process (kill-and-resume test).
--- 4c-2 description follows ---
Stage 4c-2: catalysis and self-sustaining networks (level L1) on top of the 4c-1 polymer
chemistry (spec docs/stage4c2_spec.md). Copy of chem4c1.py (md5 19da9687) with these changes:
  * L_MAX = 64 (17.6% of 4c-1 chains sat at 32);
  * eight layer-C dimensions appended (37 in all): L_cat (3-8), p0 (1e-6..1, log), lam (0.3-30,
    log), kappa (1e1..1e6, log), Ea_unc (40-120 kJ/mol), r_Ea (0.3-1), m_clay and m_sulf
    (0.1-1000, log);
  * uncatalysed condensation (growth and initiation) gets an Arrhenius factor
    exp(-Ea_unc/R (1/T - 1/298.15)), k_c being its rate at 25 C (4c-1 had no T dependence here);
    catalysed ligation uses Ea_cat = r_Ea * Ea_unc, so it is never more T-sensitive;
  * mineral catalysis: every condensation (growth, initiation, catalysed ligation) is
    multiplied by f_min = 1 + m_type * mineral availability (m_type = 0 / m_clay / m_sulf);
  * catalysed ligation among the abundant species of each site (cat4c2.py): fixed-hash
    catalysis C: A + B -> AB with p_cat(L) = p0 (1 - exp(-(L - L_cat + 1)/lam)), rate kappa times
    the uncatalysed condensation per unit catalyst concentration (at 25 C), times
    exp(-Ea_cat/R (1/T - 1/298.15)). The species list and the catalysed triples are refreshed
    every REFRESH steps (counts and slot positions every step); catalysis starts in a run once
    E1 has been reached;
  * with catalysis on, runs last at least 10 dimer lifetimes (persistence needs time);
  * the bond-break step limit uses the longest chain present in the run (+ G_CAP), not L_MAX;
  * L1 detector: at 16 checkpoints (t_end * k/16) the abundant species (>= 2 copies, top
    NX_RAF) of NS_RAF sampled sites are stored; at the end the maxRAF of each checkpoint is
    computed with that checkpoint's chemistry, and the persistent RAF of a window of consecutive
    checkpoints ending at t_end is the maxRAF of the catalysed reactions (reactants, product and
    catalyst) present in the RAF at every checkpoint of the window (see _l1_from_history).
    L1_loose: persistent RAF >= 3 reactions over >= 5 product turnover times, E1 and liquid
    water. L1_strict additionally needs the RAF products to be more abundant at the end than in
    a paired run without catalysis (p0 = 0, same condition and seed): paired_control().
    L1 is NOT nested in E2 (owner/coordinator decision, 26 Oct 2026): it needs E1 only;
  * composition fields: RAF size and components, polymer-catalysed share of condensation
    events (run and RAF site), RAF product mass including chains that carry a RAF product as
    prefix, catalyst lengths, nucleotide share of RAF products, L1 time relative to E2.
Limitations: RAFs can only form among species with repeated copies, i.e. short oligomers with
alphabet 2 (single-copy sequences neither catalyse nor are catalytically ligated); kappa
saturates above about 1e4 (catalysis is limited by the instances available); uncatalysed
chain-chain ligation is not simulated (as in 4c-1; it only sets the unit of kappa).
Special overrides added: "salt_per_step" (re-hash catalysis every step: shuffled control).
--- 4c-1 description follows ---
Stage 4c-1: bottom-up polymer chemistry (levels E1 and E2), spec docs/stage4c1_spec.md.

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
import os
import numpy as np
import torch
import cat4c3 as cat4c2   # cat4c2 with deterministic tie-breaking (see cat4c3.py)
import raf4c2
import cat4c4   # compartment-aware catalysis (label in the species identity)

torch.set_num_threads(4)
# Bit-identical results in every process (kill-and-resume test): with copying, some GPU kernels
# (atomic accumulation) gave run-to-run differences at screening fidelity; the deterministic
# implementations remove them (checked: 3 processes x 2 runs identical, 26 Oct 2026).
torch.use_deterministic_algorithms(True, warn_only=True)
XPU_FRACTION = float(os.environ.get("LOT_XPU_FRACTION", "0.875"))
if torch.xpu.is_available():
    # VRAM cap: 0.875 of the 32 GB card (28 GB; raised from 0.5 on 27 Oct 2026 at the owner's request:
    # one GPU job at a time is enforced by the GPU lock in thermal.py). Without a cap the caching
    # allocator keeps ~26 GB reserved, and two such processes overflowed the card into system RAM.
    torch.xpu.set_per_process_memory_fraction(XPU_FRACTION)

NS_RAF = 8                # sites per run whose maxRAF / lineages are recorded at each checkpoint
REFRESH = 1               # steps between refreshes of the abundant-species list and triples
CAT_AFTER_E1 = True       # catalysis starts in a run once E1 has been reached (False: from t = 0)
# NS_RAF, REFRESH and CAT_AFTER_E1 are read at call time (4c-2 adopted 8 / 1 / True: sens4c2)
NX_RAF = 128              # abundant species per site in the RAF detector (dynamics use cat4c2.NX)
N_CHECK = 16              # RAF and lineage checkpoints per run
L_H = 6                   # minimum lineage length in the L2 detector
NX_L2 = 128               # abundant sequences per site and checkpoint in the L2 detector
L2_SHARE = 0.2            # information share that counts as persisting
L2_SPAN = 5.0             # turnover times a lineage must persist (L2_loose)
L2_MARGIN = 0.15          # share above the copy_random controls (L2_strict; null-calibrated, see validate4c3 test 8)
PRIMER_LENS = (2, 3)      # lengths of free chains that can prime a copy (read at call time)
L2_KEY_LEN = None         # lineage identity: None = the whole sequence; k = its first k letters
                          # (copies start at the template's first letter, and growth or ligation
                          # only append letters, so a prefix survives them); chains need L >= max(L_H, k)
COPY_TIME = True          # a template stays busy for l / r hours after starting a copy of l letters
# P3-1b switches (Review 4: all default False; False = chem4c4 bit for bit). Read at call time.
MULTI_START = False       # several copy starts per template and step (renewal rule)
SLOT_ROUND_FIX = False    # slot sharing cuts copies only where demand exceeds free slots (always on under MULTI_START)
POISSON_INV = False       # Poisson draws by exact CDF inversion (always on under MULTI_START)
SLOT_CUT_RANDOM = False    # Review 9 (as chem_p31b, Review 7 A2b): old-path slot cut in random order (aux generator)
MS_MELT_FLOOR = True       # Review 9 (A2c False): no half-period floor on t_melt under MULTI_START
MS_ORDER_SLOT = False      # Review 9 (A2d/A2e): MULTI_START copies ordered by template slot index instead of at random
EMPTY_CACHE_EVERY = int(os.environ.get("LOT_EMPTY_CACHE_EVERY", "0"))   # memory only (as chem_p31b): release the XPU cache every n steps
COPY_CONC_BINS = 4096      # Review 9 record: lineage bins for the copy concentration record (0 = off); no effect on dynamics
DUPLEX_HYDRO_FACTOR = 1.0  # A4: hydrolysis of paired strands x this factor (A4: 0.1; UV term unchanged)
COPY_CHUNK = 1 << 18      # templates per chunk in the letter computation of the copy step (memory only)
K0_OFF = 1e6              # per hour, duplex dissociation prefactor
# --- 4c-4 vesicle layer ---
NV = 8                    # vesicle slots per site
A0_VES = 100.0            # amphiphile molecules in a newborn vesicle membrane (spec leaves a0 open)
K_AS = 1e3               # per hour, assembly / disassembly around the CVC
K_FORM = 10.0             # per hour, vesicle formation from non-vesicle membrane
K_GROW = 10.0             # per hour, uptake of non-vesicle membrane by existing vesicles
K_EX = 10 ** -1.5         # per hour, membrane exchange between vesicles (fixed, geometric middle)
F_PERM = 0.1              # activated-monomer access inside vesicles (fixed, geometric middle)
K_BREAK = 1.0             # per hour, breakdown outside the stability window (spec: layer C; see notes)
L_RETAIN = 4              # chains of at least this length stay inside a vesicle

A_NAMES = ["T", "P", "pH", "H2", "UV", "salinity", "Mg", "feed", "x_nuc", "amph", "tau", "mix",
           "evap", "thermo", "depth", "period", "mineral", "mineral_type", "avail_time"]
C_NAMES = ["k_c", "f_dry", "aw_dry", "Y_UV", "Y_H2", "Ea_deact", "theta", "p_min", "s_exp", "k_UV"]
C2_NAMES = ["L_cat", "p0", "lam", "kappa", "Ea_unc", "r_Ea", "m_clay", "m_sulf"]
C3_NAMES = ["L_t", "k_t", "eps", "dH_bp", "T_m8"]
C4_NAMES = ["CVC0", "gamma_T", "pKa_app", "w_pH", "Mg_crit", "sal_crit", "T_mem", "phi0", "f_rehyd", "beta",
            "s_shear", "s_cycle"]
NAMES = A_NAMES + C_NAMES + C2_NAMES + C3_NAMES + C4_NAMES
DIM = len(NAMES)          # 54
INACTIVE = []

V = 1000.0                # molecules per site per mM
C_WORK = 0.01             # mM, working concentration of activated monomers (E1)
L_MAX = 64
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
        "L_cat": torch.clamp(3 + torch.floor(6 * u[:, 29]), max=8), "p0": lg(30, -6, 0),
        "lam": lg(31, math.log10(0.3), math.log10(30)), "kappa": lg(32, 1, 6), "Ea_unc": 40 + 80 * u[:, 33],
        "r_Ea": 0.3 + 0.7 * u[:, 34], "m_clay": lg(35, -1, 3), "m_sulf": lg(36, -1, 3),
        "L_t": torch.clamp(4 + torch.floor(9 * u[:, 37]), max=12), "k_t": lg(38, 1, 5),
        "eps": lg(39, -3, math.log10(0.3)), "dH_bp": 25 + 20 * u[:, 40], "T_m8": 10 + 80 * u[:, 41],
    }
    if u.shape[1] >= 54:                                   # 4c-4 layer C (absent in 42-dim u)
        p.update({"CVC0": lg(42, -2, 2), "gamma_T": 0.05 * u[:, 43], "pKa_app": 6 + 4 * u[:, 44],
                  "w_pH": 0.3 + 1.7 * u[:, 45], "Mg_crit": lg(46, math.log10(5e-4), math.log10(5e-2)),
                  "sal_crit": lg(47, math.log10(0.03), math.log10(3)), "T_mem": 40 + 60 * u[:, 48],
                  "phi0": lg(49, -3, math.log10(0.3)), "f_rehyd": 1 + 9 * u[:, 50], "beta": 10 * u[:, 51],
                  "s_shear": lg(52, -2, 1), "s_cycle": lg(53, -2, 1)})
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
        "arr_unc": torch.exp(-p["Ea_unc"][:, None] * 1000 / 8.314 * (1 / (T_ph + 273.15) - 1 / 298.15)),
        "arr_cat": torch.exp(-(p["r_Ea"] * p["Ea_unc"])[:, None] * 1000 / 8.314
                             * (1 / (T_ph + 273.15) - 1 / 298.15)),
    }
    for k in list(t):
        t[k] = torch.cat([t[k], t[k].mean(1, keepdim=True)], 1)
    # duplex: dG(l)/RT = (l - 1) * dup_e with dup_e = dH_bp / R (1/T_m8 - 1/T) per phase [B,2]
    # (negative below T_m8); a cycle-averaged phase uses the mean of the two k_off values
    t["dup_e"] = (p["dH_bp"] * 1000 / 8.314)[:, None] * (1 / (p["T_m8"] + 273.15)[:, None] - 1 / (T_ph + 273.15))
    return t, E, theta, T_ph


TURNOVER_PROTECT = True   # turnover times include adsorption protection (False: as in 4c-2)


def turnover_rate(L, k_h, k_uv, w, theta, p_min):
    """Loss rate (per hour) of a chain of length L: (L - 1)(k_h (1 - theta p_min) + k_UV) +
    (1 - theta) / tau, theta = adsorbed fraction (adsorbed chains are protected from hydrolysis by
    p_min and from washout). With TURNOVER_PROTECT False: (L - 1)(k_h + k_UV) + 1/tau (4c-2)."""
    if not TURNOVER_PROTECT:
        return (L - 1) * (k_h + k_uv) + w
    return (L - 1) * (k_h * (1 - theta * p_min) + k_uv) + (1 - theta) * w


def k_off(l, e_ph, ph):
    """Duplex dissociation rate (per hour) for duplex lengths l [n] (float), per-row exponent
    table e_ph [n,2] (wet, dry) and phase index ph [n] (0 wet, 1 dry, 2 cycle average)."""
    a = torch.exp(torch.clamp((l - 1)[:, None] * e_ph, max=60.0))         # [n,2]
    pick = torch.where(ph == 1, a[:, 1], a[:, 0])
    return K0_OFF * torch.where(ph == 2, a.mean(1), pick)


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
    g = p["k_c"] * t["ceff"][:, 2] * t["gmult"][:, 2] * t["arr_unc"][:, 2] * ms / V
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


POISSON_EXACT_MAX = 1e4   # above this rate a Poisson draw uses the normal approximation
POISSON_CAP = 1e7         # larger rates are capped (pools and slots are far smaller anyway)


INV_CAP = 200             # Review 6: hard cap on the inversion search for lam > 10 (error, never a silent fallback)
INV_STATS = {"max_iter": 0}


def _poisson_inv(lam, gen):
    """Exact Poisson draws by CDF inversion (float64), one torch.rand uniform per element from gen. lam <= 10: the walk
    from the mode with the pmf recurrence (as before). lam > 10 (Review 6): start at the Cornish-Fisher quantile
    x0 = round(lam + sqrt(lam) z + (z^2 - 1) / 6), z = Phi^-1(u), and step to the smallest x with u <= F(x),
    F(x) = Q(x + 1, lam) (gammaincc); the number of steps is logged (INV_STATS) and capped at INV_CAP (error)."""
    u = torch.rand(lam.shape, device=lam.device, dtype=torch.float64, generator=gen)
    l = lam.double()
    big = l > 10
    x = _poisson_walk(torch.where(big, torch.zeros_like(l), l), u)
    if bool(big.any()):
        bi = big.nonzero(as_tuple=True)
        lb, ub = l[bi], u[bi]
        z = math.sqrt(2.0) * torch.erfinv(torch.clamp(2 * ub - 1, min=-1 + 1e-16, max=1 - 1e-16))
        xb = torch.clamp(torch.round(lb + torch.sqrt(lb) * z + (z * z - 1) / 6), min=0.0)
        act = torch.ones_like(ub, dtype=torch.bool)
        it = 0
        while True:
            Fx = torch.special.gammaincc(xb + 1, lb)
            Fm = torch.where(xb > 0, torch.special.gammaincc(xb, lb), torch.zeros_like(lb))
            dn = act & (ub <= Fm)
            up = act & (ub > Fx)
            act = dn | up
            if not bool(act.any()):
                break
            it += 1
            if it > INV_CAP:
                raise RuntimeError(f"_poisson_inv: no convergence after {INV_CAP} steps (lam max {float(lb.max())})")
            xb = torch.where(dn, xb - 1, torch.where(up, xb + 1, xb))
        INV_STATS["max_iter"] = max(INV_STATS["max_iter"], it)
        x = x.clone()
        x[bi] = xb
    return x.to(lam.dtype)


def _poisson_walk(l, u):
    """The pmf-recurrence walk from the mode (the P3-1b Review 3 sampler) for given uniforms u."""
    pos = l > 0
    lp = torch.where(pos, l, torch.ones_like(l))
    m = torch.floor(lp)
    F0 = torch.where(pos, torch.special.gammaincc(m + 1, lp), torch.ones_like(l))
    p0 = torch.where(pos, torch.exp(-lp + m * torch.log(lp) - torch.lgamma(m + 1)), torch.ones_like(l))
    m = torch.where(pos, m, torch.zeros_like(m))
    x = m.clone()
    nmax = int(10 * math.sqrt(float(l.max())) + 100) if l.numel() else 0
    act = pos & (u > F0)                                   # walk up
    k, p, F = m.clone(), p0.clone(), F0.clone()
    for _ in range(nmax):
        if not bool(act.any()):
            break
        k = torch.where(act, k + 1, k)
        p = torch.where(act, p * lp / k, p)
        F = torch.where(act, F + p, F)
        x = torch.where(act, k, x)
        act = act & (u > F) & (p > 0)
    act = pos & (u <= F0) & (m > 0)                        # walk down
    p, F = p0.clone(), F0.clone()
    for _ in range(nmax):
        if not bool(act.any()):
            break
        F1 = F - p
        go = act & (u <= F1)
        F = torch.where(go, F1, F)
        p = torch.where(go, p * x / lp, p)
        x = torch.where(go, x - 1, x)
        act = go & (x > 0)
    return x


def _poisson(lam, gen, aux):
    """Poisson draws that stay fast for huge rates: exact up to POISSON_EXACT_MAX, normal
    approximation (from the auxiliary generator) above, rates capped at POISSON_CAP. A huge
    rate (initiation can reach 1e30 with extreme concentration factors) made the XPU sampler
    run for many minutes."""
    lam = torch.clamp(_rate(lam), max=POISSON_CAP)
    if POISSON_INV or MULTI_START:
        x = _poisson_inv(torch.clamp(lam, max=POISSON_EXACT_MAX), gen)
    else:
        x = torch.poisson(torch.clamp(lam, max=POISSON_EXACT_MAX), generator=gen)
    big = lam > POISSON_EXACT_MAX
    if bool(big.any()):
        z = torch.randn(lam.shape, device=lam.device, generator=aux)
        x = torch.where(big, torch.clamp(torch.round(lam + torch.sqrt(lam) * z), min=0.0), x)
    return x


def _mask(nbits):
    one = torch.ones_like(nbits)
    return torch.where(nbits >= 64, torch.full_like(nbits, -1), torch.bitwise_left_shift(one, nbits) - 1)


def _comp_key(L, seq):
    """Species keys of the complements (chemistry 0) of sequences seq of length L."""
    return cat4c2.species_key(torch.zeros_like(L), L, torch.bitwise_not(seq) & _mask(L * BITS))


def _l2_min_len():
    return max(L_H, L2_KEY_LEN or 0)


def _lineage_view(L, seq):
    """(length, sequence) that identifies a chain's lineage: the whole chain, or its prefix."""
    if not L2_KEY_LEN:
        return L, seq
    k = torch.full_like(L, int(L2_KEY_LEN))
    return k, seq & _mask(k * BITS)


def _lineage_keys(L, seq, chem):
    """Rows [n,S] -> [n, 2 NX_L2] keys of the abundant nucleotide-like sequences of length >= L_H
    (>= 2 copies, top NX_L2) and of their complements; empty entries are I64MAX."""
    Lk, sk = _lineage_view(L, seq)
    Lh = torch.where((chem == 0) & (L >= _l2_min_len()), Lk, torch.zeros_like(L))
    X = cat4c2.abundant(Lh, sk, chem, nx=NX_L2)
    v = X["valid"]
    big = torch.full_like(X["key"], cat4c2.I64MAX)
    k = torch.where(v, X["key"], big)
    kc = torch.where(v, _comp_key(X["L"], X["seq"]), big)
    out = torch.cat([k, kc], 1)
    if out.shape[1] < 2 * NX_L2:
        out = torch.cat([out, torch.full((out.shape[0], 2 * NX_L2 - out.shape[1]), cat4c2.I64MAX,
                                         dtype=out.dtype, device=out.device)], 1)
    return out


def _l2_shares(L, seq, chem, rows, keys):
    """Information share [n, N_CHECK]: for the sites `rows` [n], the share of chain mass
    (chemistry 0, L >= L_H) whose key is in that site's stored keys [n, N_CHECK, K]."""
    Ls, ss, cs = L[rows], seq[rows], chem[rows]
    m = (cs == 0) & (Ls >= _l2_min_len())
    Lk, sk = _lineage_view(Ls, ss)
    k = torch.where(m, cat4c2.species_key(cs, Lk, sk), torch.full_like(Ls, -1))
    mass = torch.where(m, Ls, torch.zeros_like(Ls)).float()
    tot = mass.sum(1)
    out = torch.zeros(rows.numel(), N_CHECK, device=L.device)
    for c in range(N_CHECK):
        srt = torch.sort(keys[:, c, :], dim=1).values
        pos = torch.searchsorted(srt, k).clamp(max=srt.shape[1] - 1)
        hit = torch.gather(srt, 1, pos) == k
        out[:, c] = (mass * hit.float()).sum(1) / tot.clamp(min=1.0)
    return out


@torch.no_grad()
def simulate(u, generator, fidelity="screen", record=False, override=None, init=None, probe=None):
    """chem4c3.simulate plus the 4c-4 vesicle layer (module docstring). Extra special overrides
    (validation): "amph_total" (initial amphiphile per site, no precursor), "amph_closed" (no
    amphiphile feed or washout, no vesicle washout), "window_close_t" (window closed from this
    time on), "window_open" (force the window open), "ves_log" (log formations and divisions),
    "init_vesicles" ({"a": [...] in units of a0, "n": [...] chains each, "L", "seq"}),
    "init_duplexes" ({"n", "L", "seq"} paired strands outside), "no_division", "no_formation",
    "ves_state" (return the final vesicle state); "init_vesicles" may carry "pairs" (duplexes per vesicle)."""
    return _simulate(u, generator, fidelity, record, override, init, probe)


@torch.no_grad()
def _simulate(u, generator, fidelity="screen", record=False, override=None, init=None, probe=None):
    """Simulate B environments. Returns dict of per-run results (CPU tensors / numpy).

    override: dict of physical values replacing decoded ones (tensors [B] or scalars), plus
      special keys for validation only: "clamp_mstar" (hold M* per site at this count, split by
      x_nuc share), "t_end" (fixed run length in hours), "no_cycle_average" (always resolve
      phases), "empty_start" (pools start at zero), "k_h_fixed" (per-bond hydrolysis rate per
      hour, all phases), "hydro_mult" (multiply hydrolysis), "init_factor" (multiply initiation).
    init: {"dimers": n} start every site with n dimers of the dominant chemistry (validation).
    probe: list over runs of None or (site index, iterable of species keys): the returned
      "probe_abundance" is the final total count of those species in that site (used for the
      paired no-catalysis control of the L1 self-sustaining test).
    """
    dev = u.device
    gen = generator
    cat4c2.HASH_DEVICE = dev
    aux = torch.Generator(device=dev)
    aux.manual_seed(generator.initial_seed() + 7919)
    N, S, n_steps = FIDELITY[fidelity] if isinstance(fidelity, str) else fidelity
    B = u.shape[0]
    R = B * N
    p = decode(u)
    override = dict(override or {})
    special = {k: override.pop(k) for k in ("clamp_mstar", "t_end", "no_cycle_average", "empty_start",
                                            "k_h_fixed", "hydro_mult", "init_factor", "salt_per_step",
                                            "copy_random", "copy_log", "init_templates", "return_state",
                                            "amph_total", "amph_closed", "window_close_t", "window_open", "ves_log",
                                            "init_vesicles", "init_duplexes", "no_division", "no_formation", "ves_state")
               if k in override}
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
    ves_on = "CVC0" in p and bool((p["amph"] > 0).any())          # hard guard: no amphiphiles -> chem4c3
    if ves_on:
        share = share * (1 - p["amph"])[:, None]                   # amphiphile precursors take their share
    pools0 = _initial_pools(p, tab, E, share) if not special.get("empty_start") else torch.zeros(B, 2, 3, dtype=torch.float64)
    t_end, dt_target = _time_end(p, tab, E, pools0, n_steps)
    # L1 persistence needs >= 5 turnover times of short RAF members: with catalysis on, runs last
    # at least 10 dimer lifetimes (1 / (k_h + k_UV + 1/tau)), still capped by the available time
    k_dimer = torch.maximum(tab["knuc"][:, 2], tab["kpep"][:, 2]) + p["k_UV"] * p["UV"] + 1 / p["tau"]
    t_end = torch.where((p["p0"] > 0) | (p["k_t"] > 0),
                        torch.minimum(torch.maximum(t_end, 10 / k_dimer), p["avail_time"]), t_end)
    if ves_on:
        in_win = ((p["pH"] - p["pKa_app"]).abs() <= p["w_pH"]) & (p["Mg"] <= p["Mg_crit"]) & \
                 (p["salinity"] <= p["sal_crit"]) & (T_ph.max(1).values <= p["T_mem"])
        if special.get("window_open"):
            in_win = torch.ones_like(in_win)
        # no run extension (coordinator decision 28 Oct 2026): horizon, steps and dt stay chem4c3's
    dt_target = t_end / n_steps
    if "t_end" in special:
        t_end = torch.as_tensor(special["t_end"], dtype=torch.float64).expand(B).clone()
        dt_target = t_end / n_steps
    half = p["period"] / 2
    cycling = p["depth"] > 0
    resolved = cycling & ((half >= dt_target) | bool(special.get("no_cycle_average", False)))
    liquid = p["P"] >= p_sat_bar(T_ph.max(1).values)
    w = 1 / p["tau"]
    kuv = p["k_UV"] * p["UV"]
    m_type = torch.stack([torch.zeros_like(p["m_clay"]), p["m_clay"], p["m_sulf"]], 1).gather(
        1, p["mineral_type"].long()[:, None]).squeeze(1)
    f_min = 1 + m_type * p["mineral"]
    cat_on = p["p0"] > 0
    salt_step = bool(special.get("salt_per_step", False))
    copy_random = bool(special.get("copy_random", False))
    copy_log = [] if special.get("copy_log") else None
    copy_on = p["k_t"] > 0

    f32 = lambda x: x.to(dev, torch.float32)
    run_of_row = torch.arange(R, device=dev) // N
    # state
    lam0 = f32(pools0)[run_of_row]                                       # [R,2,3]
    F = _poisson(lam0[:, :, 0], gen, aux)
    M = _poisson(lam0[:, :, 1], gen, aux)
    Ms = _poisson(lam0[:, :, 2], gen, aux)
    L = torch.zeros(R, S, dtype=torch.int64, device=dev)
    seq = torch.zeros(R, S, dtype=torch.int64, device=dev)
    chem = torch.zeros(R, S, dtype=torch.int64, device=dev)
    if init and init.get("dimers"):
        nd = int(init["dimers"])
        L[:, :nd] = 2
        seq[:, :nd] = torch.randint(0, 1 << (2 * BITS), (R, nd), device=dev, generator=gen)
        chem[:, :nd] = (share[:, 1] > share[:, 0]).long().to(dev)[run_of_row][:, None]
    it = special.get("init_templates")
    if it:
        ns_, Ls_, ss_ = (list(it[k]) if isinstance(it[k], (list, tuple)) else [it[k]] for k in ("n", "L", "seq"))
        j0 = 0
        for nt, lt, st in zip(ns_, Ls_, ss_):         # one or several template species
            L[:, j0:j0 + int(nt)] = int(lt)
            seq[:, j0:j0 + int(nt)] = int(st)
            chem[:, j0:j0 + int(nt)] = 0
            j0 += int(nt)
    pair = torch.full((R, S), -1, dtype=torch.int64, device=dev)    # partner slot of a paired strand
    busy = torch.zeros(R, S, device=dev)                             # hours a template is still busy copying
    dlen = torch.zeros(R, S, dtype=torch.int64, device=dev)         # duplex length (paired strands)
    rowS = torch.arange(R, device=dev)[:, None] * S
    lab = torch.zeros((R, S) if ves_on else (1, 1), dtype=torch.int8, device=dev)   # compartment label (4c-4)

    def unpair(m):
        """Dissolve the duplexes of the slots in bool mask m [R,S] (both strands)."""
        nonlocal pair, dlen
        m = m & (pair >= 0)
        if not bool(m.any()):
            return
        pf, df = pair.view(-1), dlen.view(-1)
        partner = (rowS + pair)[m]
        pf[partner] = -1
        df[partner] = 0
        pair = torch.where(m, torch.full_like(pair, -1), pair)
        dlen = torch.where(m, torch.zeros_like(dlen), dlen)
    clamp = special.get("clamp_mstar")
    if clamp is not None:
        clamp = f32(share * torch.as_tensor(clamp, dtype=torch.float64).expand(B)[:, None])[run_of_row].round()
        Ms = clamp.clone()

    if ves_on:
        VS = Vesicles(p, tab, T_ph, E, in_win, B, N, R, S, dev, run_of_row, special, generator, record, cycling)
        VS.init_chains(L, seq, chem, pair, dlen, lab, special)
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
    cat_prev = torch.zeros(B, dtype=torch.bool)
    cat_rows = torch.zeros(0, dtype=torch.long, device=dev)
    Xcat = None
    trip = None
    cat_row = torch.zeros(R, device=dev)
    grow_row = torch.zeros(R, device=dev)
    cond_cum = torch.zeros(B, dtype=torch.float64)
    cat_events = torch.zeros(B, dtype=torch.float64)
    cond_events = torch.zeros(B, dtype=torch.float64)
    sat_triples = 0
    raf_rows = torch.cat([torch.arange(b * N, b * N + min(NS_RAF, N)) for b in range(B)]).to(dev)
    raf_hist = [[[] for _ in range(min(NS_RAF, N))] for _ in range(B)]   # per run, site: list of (t, result)
    next_chk = torch.ones(B, dtype=torch.long)
    # copying counters (per run) and the L2 lineage record
    cp_n = torch.zeros(B, dtype=torch.float64)
    cp_len = torch.zeros(B, dtype=torch.float64)
    cp_ext = torch.zeros(B, dtype=torch.float64)
    cp_primed = torch.zeros(B, dtype=torch.float64)
    cp_try = torch.zeros(B, dtype=torch.float64)          # released copies before the M* / slot limits
    cp_long = torch.zeros(B, dtype=torch.float64)         # copies of length >= L_H
    tmpl_time = torch.zeros(B, dtype=torch.float64)       # time integral of template-capable chains per site
    cp_r = torch.zeros(B, dtype=torch.float64)            # sum of r (extension rate) over released copies
    cp_try_q4 = torch.zeros(B, dtype=torch.float64)
    cp_err = torch.zeros(B, dtype=torch.float64)
    cp_q4 = torch.zeros(B, dtype=torch.float64)
    cp_drop = torch.zeros(B, dtype=torch.float64)
    cc_bins = torch.zeros(B * COPY_CONC_BINS, dtype=torch.float64, device=dev)   # Review 9: copies by template lineage (bins)
    rc_sites = torch.zeros(B, dtype=torch.float64)   # P3-1b records: site-steps / copies lost to the rounding cut
    rc_lost = torch.zeros(B, dtype=torch.float64)
    ms_rec = {k: torch.zeros(B, dtype=torch.float64) for k in ("tsteps", "occ", "slots", "blocks", "primers", "any")}
    cp_sat = torch.zeros(B, dtype=torch.float64)
    n2 = N                                                # lineages are recorded in every site (low-variance S)
    l2_rows = torch.arange(R, device=dev).view(B, n2)
    l2_keys = torch.full((B, n2, N_CHECK, 2 * NX_L2), cat4c2.I64MAX, dtype=torch.int64, device=dev)
    l2_t = torch.zeros(B, N_CHECK, dtype=torch.float64)
    next_chk2 = torch.ones(B, dtype=torch.long)
    jpos = torch.arange(64, device=dev)

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
        # break limit for the longest chain present (plus what can grow in one step), not L_MAX
        l_top = torch.clamp(L.view(B, -1).max(1).values.cpu().double() + G_CAP, min=2.0, max=float(L_MAX))
        dt = torch.minimum(dt_target, BREAK_P / ((l_top - 1) * kbr_max))
        dt = torch.minimum(dt, to_switch)
        dt = torch.minimum(dt, t_end - t)
        dt = torch.where(active, torch.clamp(dt, min=0.0), torch.zeros_like(dt))
        dtr = f32(dt)[run_of_row]

        # 1. linear monomer network (exact over dt), states F, M, M*, out
        P, J = _expm_pair(_linear_q(E, kact, kdeact, w), dt)
        Pd, Jd = f32(P)[run_of_row], f32(J)[run_of_row]                  # [R,4,4]
        arrivals = _poisson(f32(p["feed"][:, None] * share * dt[:, None])[run_of_row], gen, aux)
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
        gm = f32(gmult * g("arr_unc"))[run_of_row][:, None]
        kc = f32(p["k_c"] * f_min)[run_of_row][:, None]
        occ = L > 0
        growable = occ & (L < L_MAX)
        n_grow = torch.stack([(growable & (chem == c)).sum(1) for c in range(2)], 1).float()
        if ves_on:                          # enclosed chains draw on M* x f_perm
            wgt = torch.where(lab > 0, F_PERM, 1.0)
            n_grow = torch.stack([((growable & (chem == c)).float() * wgt).sum(1) for c in range(2)], 1)
        lam = kc * ceffr * gm * (Ms / V) * dtr[:, None]
        demand = lam * n_grow
        lam = lam * torch.clamp(Ms / torch.clamp(demand, min=1e-9), max=1.0)
        lam_c = torch.gather(lam, 1, chem)                                    # [R,S]
        if ves_on:
            lam_c = torch.where(lab > 0, lam_c * F_PERM, lam_c)
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

        # 2b. catalysed ligation among abundant species (runs with p0 > 0 that have reached E1)
        cat_now = (cat_on & ~torch.isnan(t_first_e1)) if CAT_AFTER_E1 else cat_on.clone()
        if not torch.equal(cat_now, cat_prev):
            cat_rows = (cat_now.to(dev)[run_of_row]).nonzero().squeeze(1)
            cat_prev = cat_now.clone()
            Xcat = None
        if cat_rows.numel():
            # paired strands neither catalyse nor get ligated: hidden from the species lists
            Lr, sr, cr_ = torch.where(pair >= 0, 0, L)[cat_rows], seq[cat_rows], chem[cat_rows]
            salt = steps + 1 if salt_step else 0
            if Xcat is None or salt_step or steps % REFRESH == 0:
                Xcat = cat4c4.abundant_lab(Lr, sr, cr_, lab[cat_rows]) if ves_on else cat4c2.abundant(Lr, sr, cr_)
                trip = (cat4c4.catalytic_triples_lab if ves_on else cat4c2.catalytic_triples)(Xcat, f32(p["p0"]).double()[run_of_row[cat_rows]],
                                                p["L_cat"].long().to(dev)[run_of_row[cat_rows]],
                                                f32(p["lam"]).double()[run_of_row[cat_rows]], salt, ALPHA, BITS, L_MAX)
            else:
                Xcat = cat4c4.rematch_lab(Xcat, Lr, sr, cr_, lab[cat_rows]) if ves_on else cat4c2.rematch(Xcat, Lr, sr, cr_)
            rr_run = run_of_row[cat_rows].cpu()
            arr = tab["arr_cat"].gather(1, ph[:, None]).squeeze(1)
            r_cm = f32(p["k_c"] * f_min * p["kappa"] * arr * gmult * ceff ** 2 / V ** 2)[run_of_row[cat_rows]]
            phi = (Ms.sum(1) / torch.clamp(Ms.sum(1) + M.sum(1), min=1.0))[cat_rows]
            ev, sat = cat4c2.apply_catalysis(L, seq, chem, Ms, cat_rows, Xcat, trip, r_cm, r_cm * phi,
                                             dtr[cat_rows], gen, aux, _poisson, ALPHA, BITS, L_MAX)
            cat_events += torch.zeros(B, dtype=torch.float64).index_add_(0, rr_run, ev.double().cpu())
            cat_row.index_add_(0, cat_rows, ev.float())
            sat_triples += sat
        grow_row += k_add.sum(1).float()
        cond_events += grow_row.view(B, N).sum(1).double().cpu() - cond_cum
        cond_cum = grow_row.view(B, N).sum(1).double().cpu()

        # 2c'. duplex dissociation (the copies themselves are made in step 3, after breaks)
        phr = ph.to(dev)[run_of_row]                                      # [R]
        e_r = f32(tab["dup_e"])[run_of_row]                               # [R,2]
        top = (pair >= 0) & (pair > torch.arange(S, device=dev)[None, :])   # one entry per duplex
        rr, cc = top.nonzero(as_tuple=True)
        u_off = torch.rand(rr.numel(), device=dev, generator=gen)
        if rr.numel():
            ko = k_off(dlen[rr, cc].float(), e_r[rr], phr[rr])
            off = u_off < 1 - torch.exp(-ko * dtr[rr])
            m = torch.zeros(R, S, dtype=torch.bool, device=dev)
            m[rr[off], cc[off]] = True
            unpair(m)

        # 4. per chain: adsorption, washout, one bond break
        th = f32(theta)[run_of_row][:, None]
        ads = torch.rand(R, S, device=dev, generator=gen) < th
        kb = torch.where(chem == 0, f32(knuc)[run_of_row][:, None], f32(kpep)[run_of_row][:, None])
        if DUPLEX_HYDRO_FACTOR != 1.0:      # A4 (P3-1b)
            kb = torch.where(ads, kb * f32(1 - p["p_min"])[run_of_row][:, None], kb)
            kb = torch.where(pair >= 0, kb * DUPLEX_HYDRO_FACTOR, kb) + f32(kuv)[run_of_row][:, None]
        else:
            kb = torch.where(ads, kb * f32(1 - p["p_min"])[run_of_row][:, None], kb) + f32(kuv)[run_of_row][:, None]
        occ = L > 0
        pb = 1 - torch.exp(-(L - 1).clamp(min=0).float() * kb * dtr[:, None])
        pw = torch.where(ads, torch.zeros_like(kb), (1 - torch.exp(-f32(w)[run_of_row] * dtr))[:, None].expand(R, S))
        if ves_on:
            pw = torch.where(lab > 0, torch.zeros_like(pw), pw)
        u2 = torch.rand(R, S, device=dev, generator=gen)
        wash = occ & (u2 < pw)
        # a duplex is washed out as one: a strand goes if its partner goes
        wash = wash | ((pair >= 0) & wash.view(-1)[(rowS + pair.clamp(min=0)).view(-1)].view(R, S))
        v = (u2 - pw) / torch.clamp(1 - pw, min=1e-12)
        brk = occ & ~wash & (v < pb) & (L >= 2)
        unpair(wash | brk)                  # a bond break in either strand dissolves the duplex
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

        # 3 + allocation: right fragments (>= 2) take free slots first; template copies (2c) and
        # initiations share the remaining room in proportion to their numbers (both thinned)
        need = brk & (right_len >= 2)
        free = L == 0
        fc = free.sum(1)
        n_right = need.sum(1)
        n_right_ok = torch.minimum(n_right, fc)
        room = torch.clamp(fc - n_right, min=0).float()
        # 2c. template copying (chemistry 0 only) on the strands left after breaks and washout.
        # Starts per free template, all at the template's first letter: unprimed (one monomer,
        # rate r) or primed by a free chain of chemistry 0 and length in PRIMER_LENS that is the
        # complement of the template's first l_p letters: rate k_on * n_match(l_p), k_on = r / M*_0
        # (= k_on n_primer 2^-l_p for a uniform primer pool). At most one start per template and step.
        k_on = f32(p["k_c"] * f_min * p["k_t"] * g("arr_cat") * gmult * ceff / V)[run_of_row]          # [R]
        r_cp = k_on * Ms[:, 0]                                                                         # [R]
        tcap = (L >= f32(p["L_t"]).long()[run_of_row][:, None]) & (chem == 0)
        tmpl_time += dt * tcap.sum(1).view(B, N).double().mean(1).cpu()
        tmpl = tcap & (pair < 0)
        tmpl = tmpl & f32(copy_on.double()).bool()[run_of_row][:, None]
        ffree = None
        if COPY_TIME:
            if MULTI_START:     # share of the step in which a template is not still busy with an earlier copy
                ffree = torch.clamp(1 - busy / torch.clamp(dtr[:, None], min=1e-30), min=0.0, max=1.0)
            busy = torch.where(L == 0, torch.zeros_like(busy), torch.clamp(busy - dtr[:, None], min=0.0))
            tmpl = tmpl & (busy <= 0)
        plens = tuple(sorted(PRIMER_LENS))
        poff = {lp: sum(1 << k for k in plens if k < lp) for lp in plens}      # code offset per primer length
        NCODE = sum(1 << k for k in plens)
        isp = torch.zeros_like(L, dtype=torch.bool)
        pcode = torch.zeros_like(L)
        for lp in plens:
            isp = isp | (L == lp)
            pcode = torch.where(L == lp, poff[lp] + (seq & ((1 << lp) - 1)), pcode)
        prim = (chem == 0) & (pair < 0) & isp & f32(copy_on.double()).bool()[run_of_row][:, None]
        prim = prim & (L < f32(p["L_t"]).long()[run_of_row][:, None])        # a template is never also a primer
        prow = torch.arange(R, device=dev)[:, None].expand(R, S)
        NC1 = NV + 1 if ves_on else 1                     # compartments per site (primers pair within one)
        cmp = lab.long() if ves_on else torch.zeros_like(L)
        prc = prow * NC1 + cmp
        pcnt = torch.zeros(R * NC1 * NCODE, device=dev).index_add_(0, prc[prim] * NCODE + pcode[prim],
                                                            torch.ones(int(prim.sum()), device=dev)).view(R * NC1, NCODE)
        trr, tcc = tmpl.nonzero(as_tuple=True)
        nT = trr.numel()
        u_start = torch.rand(nT, device=dev, generator=gen)
        u_type = torch.rand(nT, device=dev, generator=gen)
        u_ppos = torch.rand(nT, device=dev, generator=gen)
        u_pos = torch.rand(nT, device=dev, generator=gen)
        u_ext = torch.rand(nT, 63, device=dev, generator=gen)
        u_err = torch.rand(nT, 64, device=dev, generator=gen)
        rnd_bits = torch.randint(0, 2, (nT, 64), device=dev, generator=gen)
        nc = torch.zeros(R, dtype=torch.long, device=dev)
        c_need = torch.zeros(R, device=dev)
        n_mono_start = torch.zeros(R, device=dev)
        n_mono_rel = torch.zeros(R, device=dev)
        n = 0
        ms_t = None
        if nT:
            Lt = L[trr, tcc]
            tsq = seq[trr, tcc]
            tcm = cmp[trr, tcc]                           # template compartment
            trc = trr * NC1 + tcm
            # matching primers for each primer length: copies start at the template's first
            # letter, so a primer must be the complement of the template's first l_p letters
            # (spec change approved 26 Oct 2026). Only position 0 is computed ([nT], not [nT, 64]).
            rho = []
            for lp in plens:
                segc = torch.bitwise_not(tsq) & ((1 << lp) - 1)
                cntp = torch.gather(pcnt[trc], 1, (poff[lp] + segc)[:, None]).squeeze(1) * (Lt >= lp).float()
                rho.append(k_on[trr] * cntp)
            r_m = r_cp[trr]
            if ves_on:                                     # extension inside draws on M* x f_perm
                r_m = torch.where(tcm > 0, r_m * F_PERM, r_m)
            lam_t = r_m + sum(rho)
            started = u_start < 1 - torch.exp(-lam_t * dtr[trr])
            cp_sat += torch.zeros(B, dtype=torch.float64).index_add_(
                0, run_of_row[trr].cpu(), ((lam_t * dtr[trr]) > 1).double().cpu())
            if MULTI_START:
                # P3-1b renewal rule (as chem_p31b, 4c-3 branch): cycle = 1 / lam + t_copy + t_melt; n ~ Poisson(dt_free /
                # cycle) capped at ceil(dt_free / (t_copy + t_melt)); per site the means are scaled where the expected
                # demand exceeds free slots (unprimed), matching primers in the template's compartment (primed) or
                # activated blocks; each copy is then drawn independently below.
                nT0 = nT
                dt_cp = dtr[trr]
                s_list = [1] + list(plens)
                a_s = [torch.zeros(nT, device=dev) for _ in s_list]
                e_s = [torch.zeros(nT, device=dev) for _ in s_list]
                p_cok_u = torch.ones(nT, device=dev)
                for jj in range(1, int(Lt.max())):
                    qj = r_m / (r_m + k_off(torch.full((nT,), float(jj), device=dev), e_r[trr], phr[trr]))
                    inj = jj < Lt
                    for i, s0 in enumerate(s_list):
                        if jj == s0:
                            a_s[i] = torch.ones(nT, device=dev)
                        if jj >= s0:
                            a_s[i] = torch.where(inj, a_s[i] * qj, torch.zeros_like(qj))
                            e_s[i] = e_s[i] + a_s[i]
                    if jj == 1:
                        p_cok_u = a_s[0].clone()          # an unprimed copy needs >= 2 letters to be kept
                eb = list(e_s)
                tcp = [e_ / torch.clamp(r_m, min=1e-30) for e_ in e_s]
                wt = [r_m] + list(rho)
                lam_s = torch.clamp(lam_t, min=1e-30)
                shr = [w_ / lam_s for w_ in wt]
                t_copy = sum(sh_ * tc_ for sh_, tc_ in zip(shr, tcp))
                kfl = k_off(Lt.float(), e_r[trr], phr[trr])
                t_melt = 1.0 / torch.clamp(kfl, min=1e-30)
                if MS_MELT_FLOOR:
                    t_melt = torch.where(phr[trr] == 2, torch.maximum(t_melt, f32(half)[run_of_row][trr]), t_melt)
                dt_free = dt_cp * (ffree[trr, tcc] if ffree is not None else 1.0)
                cyc = 1.0 / lam_s + t_copy + t_melt
                mu = torch.where(lam_t > 0, dt_free / cyc, torch.zeros_like(dt_free))
                d_slot = torch.zeros(R, device=dev).index_add_(0, trr, mu * shr[0] * p_cok_u)
                ms_dslot = d_slot
                f_slot = torch.clamp(room / torch.clamp(d_slot, min=1e-30), max=1.0)[trr]
                f_pr = []
                for k_, lp in enumerate(plens):
                    gk = trc * NCODE + poff[lp] + (torch.bitwise_not(tsq) & ((1 << lp) - 1))
                    dk = torch.zeros(R * NC1 * NCODE, device=dev).index_add_(0, gk, mu * shr[k_ + 1])
                    f_pr.append(torch.clamp(pcnt.view(-1)[gk] / torch.clamp(dk[gk], min=1e-30), max=1.0))
                mu_t = [mu * shr[0] * f_slot] + [mu * shr[k_ + 1] * f_pr[k_] for k_ in range(len(plens))]
                d_blk = torch.zeros(R, device=dev).index_add_(0, trr, sum(m_ * e_ for m_, e_ in zip(mu_t, eb)))
                f_blk = torch.clamp(Ms[:, 0] / torch.clamp(d_blk, min=1e-30), max=1.0)[trr]
                mu_t = [m_ * f_blk for m_ in mu_t]
                n_ms = _poisson(sum(mu_t), gen, aux)
                cap = torch.ceil(dt_free / torch.clamp(t_copy + t_melt, min=1e-30)).clamp(max=1e7)
                hit_cap = n_ms > cap
                n_ms = torch.minimum(n_ms, cap).long()
                lim_s = (f_slot < 1) & (mu * shr[0] * p_cok_u > 0)
                lim_p = torch.zeros(nT, dtype=torch.bool, device=dev)
                for k_ in range(len(plens)):
                    lim_p = lim_p | ((f_pr[k_] < 1) & (mu * shr[k_ + 1] > 0))
                lim_b = (f_blk < 1) & (mu > 0)
                tid = torch.repeat_interleave(torch.arange(nT, device=dev), n_ms)
                ordk = trr[tid] * (1 << 31) + torch.randint(0, 1 << 31, (tid.numel(),), device=dev, generator=gen)
                if MS_ORDER_SLOT:
                    ordk = trr[tid] * (1 << 31) + tcc[tid]
                tid = tid[torch.argsort(ordk, stable=True)]
                ms_t = {"run": run_of_row[trr], "occ": hit_cap, "slots": lim_s, "blocks": lim_b, "primers": lim_p}
                ms_cyc, ms_dtf, ms_tm = cyc[tid], dt_free[tid], t_melt[tid]
                wI = [m_[tid] for m_ in mu_t]
                trr, tcc, Lt, tsq, r_m = trr[tid], tcc[tid], Lt[tid], tsq[tid], r_m[tid]
                tcm, trc = tcm[tid], trc[tid]
                rho = [rh[tid] for rh in rho]
                nT = tid.numel()
                lam_t = sum(wI)
                started = torch.ones(nT, dtype=torch.bool, device=dev)
                u_type = torch.rand(nT, device=dev, generator=gen)
                u_carry = torch.rand(nT, device=dev, generator=gen)
                u_melt = torch.rand(nT, device=dev, generator=gen)
                u_ext = torch.rand(nT, 63, device=dev, generator=gen)
                u_err = torch.rand(nT, 64, device=dev, generator=gen)
                rnd_bits = torch.randint(0, 2, (nT, 64), device=dev, generator=gen) if copy_random else None
            x = u_type * lam_t
            l0 = torch.ones_like(Lt)
            wt0, wtk = (wI[0], wI[1:]) if MULTI_START else (r_m, rho)
            acc_r = wt0.clone()
            for lp, rh in zip(plens, wtk):
                l0 = torch.where(x >= acc_r, lp, l0)
                acc_r = acc_r + rh
            l0 = torch.where(started, l0, torch.zeros_like(l0))
            # start position: 0 for every copy (u_ppos is still drawn, so the random stream is
            # unchanged from the version that chose a primed position by weight)
            p0s = torch.zeros_like(Lt)
            # primer slot: the rank-th free primer of the matching sequence in the row
            segp = torch.bitwise_not(torch.bitwise_right_shift(tsq, p0s)) & _mask(torch.clamp(l0, min=1))
            offv = torch.zeros_like(l0)
            for lp in plens:
                offv = torch.where(l0 == lp, poff[lp], offv)
            gidx = trc * NCODE + offv + segp
            primed = l0 >= 2
            pslot = torch.full((nT,), -1, dtype=torch.long, device=dev)
            if bool(primed.any()):
                pr_r, pr_c = prim.nonzero(as_tuple=True)
                gs = prc[pr_r, pr_c] * NCODE + pcode[pr_r, pr_c]
                og = torch.argsort(gs, stable=True)
                gs_s, slot_s = gs[og], (pr_r * S + pr_c)[og]
                ug = torch.where(primed, gidx, torch.full_like(gidx, -1))
                ou = torch.argsort(ug, stable=True)
                us = ug[ou]
                first = torch.searchsorted(us, us, right=False)
                rk = torch.empty_like(ou)
                rk[ou] = torch.arange(nT, device=dev) - first
                st = torch.searchsorted(gs_s, gidx)
                have = pcnt.view(-1)[gidx.clamp(min=0)].long()
                okp = primed & (rk < have)
                pslot = torch.where(okp, slot_s[(st + rk).clamp(max=max(gs_s.numel() - 1, 0))], pslot)
                if MULTI_START:
                    ms_t["primers"][tid[primed & ~okp & started]] = True
                started = started & (~primed | okp)
            maxl = Lt - p0s
            lcp = torch.where(started, l0, torch.zeros_like(l0))
            alive = started & (lcp < maxl)
            for jj in range(1, 64):
                if not bool(alive.any()):
                    break
                act = alive & (lcp == jj)
                q = r_m / (r_m + k_off(torch.full((nT,), float(jj), device=dev), e_r[trr], phr[trr]))
                cont = act & (u_ext[:, jj - 1] < q) & (lcp < maxl)
                lcp = lcp + cont.long()
                alive = alive & ~(act & ~cont)
            mono = started & (l0 == 1)
            n_mono_start = torch.zeros(R, device=dev).index_add_(0, trr, mono.float())
            cok = started & (lcp >= 2)
            n_mono_rel = torch.zeros(R, device=dev).index_add_(0, trr, (cok & mono).float())
            cneed = torch.where(cok, lcp - l0, 0).float()
            cum = torch.zeros(R, device=dev).index_add_(0, trr, cneed)
            within_c = torch.cumsum(cneed, 0) - (torch.cumsum(cum, 0) - cum)[trr]
            if MULTI_START:
                ms_t["blocks"][tid[cok & (within_c > Ms[trr, 0])]] = True
            cok = cok & (within_c <= Ms[trr, 0])
            ntry = torch.zeros(B, dtype=torch.float64).index_add_(0, run_of_row[trr].cpu(), cok.double().cpu())
            cp_try += ntry
            cp_try_q4 += ntry * (t >= win_lo).double()
            nc = torch.zeros(R, dtype=torch.long, device=dev).index_add_(0, trr, (cok & mono).long())
            c_need = torch.zeros(R, device=dev).index_add_(0, trr, torch.where(cok, lcp - l0, 0).float())
            n = nT
        Ms_i = Ms.clone()
        Ms_i[:, 0] = torch.clamp(Ms_i[:, 0] - c_need, min=0.0)
        rate_init = kc * ceffr * gm * (Ms_i / V) * Ms_i / 2 * dtr[:, None] * f32(init_f)[run_of_row][:, None]
        n_init = _poisson(rate_init, gen, aux)
        n_init = torch.clamp(torch.minimum(n_init, torch.floor(Ms_i / 2)), min=0.0)
        tot_init = n_init.sum(1)
        # slot competition in proportion to rates (unprimed copies only; a primed copy takes its
        # primer's slot): copying demand = expected unprimed starts (r dt per free template, not
        # capped at one per step) times the share of unprimed starts released in this row
        if MULTI_START:
            # P3-1b: expected unprimed copies that need a slot x the released share; never below the copies made (M8b)
            w_c = (ms_dslot if ms_t is not None else torch.zeros(R, device=dev)) * n_mono_rel / torch.clamp(n_mono_start, min=1.0)
            w_c = torch.maximum(w_c, nc.float())
        else:
            w_c = r_cp * dtr * tmpl.sum(1).float() * n_mono_rel / torch.clamp(n_mono_start, min=1.0)
        w_fix = torch.maximum(w_c, nc.float())
        if SLOT_ROUND_FIX:
            w_c = w_fix
        frac_c = w_c / torch.clamp(w_c + tot_init, min=1e-30)
        u_keep = torch.rand(R, device=dev, generator=gen)
        keep_c = torch.minimum(nc, torch.floor(torch.minimum(room, w_c + tot_init) * frac_c + u_keep).long())
        # records only: copies the fixed rule (same u_keep) would keep that this rule cuts
        frac_f = w_fix / torch.clamp(w_fix + tot_init, min=1e-30)
        keep_f = torch.minimum(nc, torch.floor(torch.minimum(room, w_fix + tot_init) * frac_f + u_keep).long())
        lost_rc = torch.clamp(keep_f - keep_c, min=0).double()
        rc_sites += torch.zeros(B, dtype=torch.float64).index_add_(0, run_of_row.cpu(), (lost_rc > 0).double().cpu())
        rc_lost += torch.zeros(B, dtype=torch.float64).index_add_(0, run_of_row.cpu(), lost_rc.cpu())
        room_i = torch.clamp(room - keep_c.float(), min=0.0)
        scale = torch.clamp(room_i / torch.clamp(tot_init, min=1), max=1.0)
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
            if ves_on:
                lab.view(-1)[tgt] = torch.where(src_len >= L_RETAIN, lab[rr, cc], torch.zeros_like(lab[rr, cc]))
        if n:
            mono_ok = cok & mono
            crank = torch.cumsum(mono_ok.long(), 0) - 1
            crank = crank - (torch.cumsum(nc, 0) - nc)[trr]
            if SLOT_CUT_RANDOM and not MULTI_START:    # random rank within the row (aux generator)
                key = trr * (1 << 31) + torch.randint(0, 1 << 31, (nT,), device=dev, generator=aux)
                o_ = torch.argsort(key, stable=True)
                cs_ = torch.cumsum(mono_ok[o_].long(), 0) - 1 - (torch.cumsum(nc, 0) - nc)[trr[o_]]
                crank = torch.empty_like(crank)
                crank[o_] = cs_
            drop = mono_ok & (crank >= keep_c[trr])
            ok = cok & ~drop
            if MULTI_START:
                ms_t["slots"][tid[drop]] = True
            cp_drop += torch.zeros(B, dtype=torch.float64).index_add_(0, run_of_row[trr].cpu(), drop.double().cpu())
            seg = torch.bitwise_right_shift(tsq, p0s * BITS) & _mask(lcp.clamp(min=0) * BITS)
            comp = torch.bitwise_not(seg) & _mask(lcp.clamp(min=0) * BITS)
            # letters of the copies, in chunks of templates (memory: [chunk, 64] instead of [nT, 64])
            wbit = torch.bitwise_left_shift(torch.ones_like(jpos), jpos)[None, :]
            pack = lambda bits: (bits.long() * wbit).sum(1)
            eps_t = f32(p["eps"])[run_of_row[trr]]
            cseq = torch.empty_like(comp)
            nerr = torch.empty_like(lcp)
            for a0 in range(0, nT, COPY_CHUNK):
                sl = slice(a0, min(nT, a0 + COPY_CHUNK))
                inpos = jpos[None, :] < lcp[sl, None]
                ext = inpos & (jpos[None, :] >= l0[sl, None])     # letters added by copying (not the primer)
                errb = (u_err[sl] < eps_t[sl, None]) & ext
                nerr[sl] = errb.sum(1)
                if copy_random:
                    # random letters for the copied part; a primer keeps its own (complementary) letters
                    cseq[sl] = (comp[sl] & ~pack(ext)) | pack(rnd_bits[sl].bool() & ext)
                else:
                    cseq[sl] = comp[sl] ^ pack(errb)
            if copy_log is not None:
                copy_log.append({k: v[ok].cpu() for k, v in (("tseq", tsq), ("tlen", Lt), ("p", p0s), ("l", lcp),
                                                                ("l0", l0), ("maxl", maxl), ("cseq", cseq), ("r", r_m),
                                                                ("ph", phr[trr]), ("run", run_of_row[trr]))})
            if bool(ok.any()):
                okm, okp2 = ok & mono, ok & ~mono
                tgt = torch.zeros(nT, dtype=torch.long, device=dev)
                if bool(okm.any()):
                    tgt[okm] = free_pos[row_off[trr[okm]] + n_right_ok[trr[okm]] + crank[okm]]
                tgt[okp2] = pslot[okp2]                            # a primed copy grows in its primer's slot
                rr2, cc2, l2, s2, t2 = trr[ok], tcc[ok], lcp[ok], cseq[ok], tgt[ok]
                if COPY_CONC_BINS:                  # Review 9 record: copies by template lineage key (both strands one lineage)
                    kt_ = cat4c2.species_key(torch.zeros_like(Lt[ok]), Lt[ok], tsq[ok])
                    kc_ = _comp_key(Lt[ok], tsq[ok])
                    lk_ = torch.remainder(torch.minimum(kt_, kc_), COPY_CONC_BINS)
                    cc_bins.index_add_(0, run_of_row[rr2] * COPY_CONC_BINS + lk_, torch.ones(lk_.numel(), dtype=torch.float64, device=dev))
                Lf, sf, cf, pf, df = L.view(-1), seq.view(-1), chem.view(-1), pair.view(-1), dlen.view(-1)
                Lf[t2] = l2
                sf[t2] = s2
                cf[t2] = 0
                if ves_on:                                 # a copy inherits its template's compartment
                    lab.view(-1)[t2] = tcm[ok].to(torch.int8)
                if MULTI_START:
                    # as chem_p31b: copies before the last one of each template have melted (free chains); the last
                    # started g = U min(dt_free, cycle) before the end: still being made if g < its copy time (template
                    # busy for the rest, stays paired), else melted with probability 1 - exp(-(g - t_copy) / t_melt)
                    ii = torch.arange(nT, device=dev)
                    lastI = torch.full((nT0,), -1, dtype=torch.long, device=dev).scatter_reduce_(0, tid[ok], ii[ok], reduce="amax")
                    is_last = torch.zeros(nT, dtype=torch.bool, device=dev)
                    is_last[lastI[lastI >= 0]] = True
                    tcl_all = torch.where(ok, (lcp - l0).float(), 0.0) / torch.clamp(r_m, min=1e-30)
                    gap = u_carry * torch.minimum(ms_dtf, ms_cyc)
                    stay = (gap < tcl_all) | (u_melt >= 1 - torch.exp(-torch.clamp(gap - tcl_all, min=0.0)
                                                                         / torch.clamp(ms_tm, min=1e-30)))
                    carry = torch.clamp(tcl_all - gap, min=0.0)
                    lm = (is_last & stay)[ok]
                    pf[t2] = torch.where(lm, cc2, torch.full_like(cc2, -1))
                    df[t2] = torch.where(lm, l2, torch.zeros_like(l2))
                    pair[rr2[lm], cc2[lm]] = t2[lm] - rr2[lm] * S
                    dlen[rr2[lm], cc2[lm]] = l2[lm]
                    if COPY_TIME:
                        lb = is_last[ok]
                        busy[rr2[lb], cc2[lb]] = carry[ok][lb]
                else:
                    pf[t2] = cc2
                    df[t2] = l2
                    pair[rr2, cc2] = t2 - rr2 * S
                    dlen[rr2, cc2] = l2
                added = (l2 - l0[ok]).float()
                if COPY_TIME and not MULTI_START:   # the template is busy for the time the copy takes
                    busy[rr2, cc2] = added / torch.clamp(r_m[ok], min=1e-30)
                Ms[:, 0] -= torch.zeros(R, device=dev).index_add_(0, rr2, added)
                runs_c = run_of_row[rr2].cpu()
                one = torch.ones(rr2.numel(), dtype=torch.float64)
                cp_n += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, one)
                cp_primed += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, (~mono[ok]).double().cpu())
                cp_len += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, l2.double().cpu())
                cp_long += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, (l2 >= L_H).double().cpu())
                cp_r += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, r_m[ok].double().cpu())
                cp_ext += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, added.double().cpu())
                cp_err += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, nerr[ok].double().cpu())
                cp_q4 += torch.zeros(B, dtype=torch.float64).index_add_(0, runs_c, (t >= win_lo)[runs_c].double())
                grow_row += torch.zeros(R, device=dev).index_add_(0, rr2, added)
                inc_site = inc_site + torch.zeros(R, device=dev).index_add_(0, rr2, added).view(B, N).mean(1).double().cpu()
        if MULTI_START and ms_t is not None:
            rn = ms_t["run"].cpu()
            for k_ in ("occ", "slots", "blocks", "primers"):
                ms_rec[k_] += torch.zeros(B, dtype=torch.float64).index_add_(0, rn, ms_t[k_].double().cpu())
            ms_rec["tsteps"] += torch.zeros(B, dtype=torch.float64).index_add_(0, rn, torch.ones(rn.numel(), dtype=torch.float64))
            anyl = ms_t["slots"] | ms_t["blocks"] | ms_t["primers"]
            ms_rec["any"] += torch.zeros(B, dtype=torch.float64).index_add_(0, rn, anyl.double().cpu())
        ni_tot = n_init.sum(1).long()
        if int(ni_tot.sum()):
            rows = torch.repeat_interleave(torch.arange(R, device=dev), ni_tot)
            start = torch.cumsum(ni_tot, 0) - ni_tot
            within = torch.arange(rows.numel(), device=dev) - start[rows]
            rank = n_right_ok[rows] + keep_c[rows] + within
            tgt = free_pos[row_off[rows] + rank]
            chem_new = (within >= n_init[rows, 0].long()).long()
            Lf, sf, cf = L.view(-1), seq.view(-1), chem.view(-1)
            Lf[tgt] = 2
            sf[tgt] = torch.randint(0, 1 << (2 * BITS), (rows.numel(),), device=dev, generator=gen)
            cf[tgt] = chem_new
            if ves_on:
                lab.view(-1)[tgt] = 0
        if clamp is not None:
            Ms = clamp.clone()

        if ves_on:
            lab.masked_fill_(L == 0, 0)
            VS.step(t, dt, ph, L, seq, chem, pair, dlen, lab, unpair, Ms)
            L, seq = VS.L, VS.seq
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
        chk_due = cat_on & (next_chk <= N_CHECK) & (t_new >= t_end * next_chk.double() / N_CHECK * (1 - 1e-12))
        if CAT_AFTER_E1:
            next_chk[chk_due & torch.isnan(t_first_e1)] += 1       # no catalysis before E1: nothing to store
            chk_due = chk_due & ~torch.isnan(t_first_e1)
        if bool(chk_due.any()):
            runs_due = chk_due.nonzero().squeeze(1)
            rows_due = torch.cat([raf_rows[b * min(NS_RAF, N):(b + 1) * min(NS_RAF, N)] for b in runs_due.tolist()])
            Xc = cat4c2.abundant(L[rows_due], seq[rows_due], chem[rows_due], nx=NX_RAF)
            Xc = {k: v.cpu() for k, v in Xc.items() if k != "sidx"}
            pres = cat4c2.present_keys(L[rows_due], seq[rows_due], chem[rows_due])
            ns = min(NS_RAF, N)
            salt_now = (steps + 1) if salt_step else 0
            for qi, b in enumerate(runs_due.tolist()):
                for si in range(ns):
                    row = qi * ns + si
                    Xr = {k: v[row] for k, v in Xc.items()}
                    raf_hist[b][si].append((float(t_new[b]), Xr, pres[row], salt_now))
            next_chk[runs_due] += 1
        chk2 = (next_chk2 <= N_CHECK) & (t_new >= t_end * next_chk2.double() / N_CHECK * (1 - 1e-12))
        if bool(chk2.any()):
            runs2 = chk2.nonzero().squeeze(1)
            rows2 = l2_rows[runs2.to(dev)].reshape(-1)
            keys2 = _lineage_keys(L[rows2], seq[rows2], chem[rows2]).view(runs2.numel(), n2, -1)
            kidx = (next_chk2[runs2] - 1).to(dev)
            l2_keys[runs2.to(dev), :, kidx] = keys2
            l2_t[runs2, next_chk2[runs2] - 1] = t_new[runs2]
            if ves_on:
                VS.checkpoint(runs2, kidx, rows2, L, lab, chem)
            next_chk2[runs2] += 1
        t = t_new
        steps += 1
        if EMPTY_CACHE_EVERY and steps % EMPTY_CACHE_EVERY == 0 and dev.type == "xpu":
            torch.xpu.empty_cache()

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
    p["_theta_eff"] = theta                         # adsorbed fraction (turnover with protection)
    l1 = _l1_from_history(raf_hist, p, tab, w, kuv, B, e1 & liquid)
    # RAF site details: catalysed share of condensation events, product mass incl. grown chains
    raf_cat_share = torch.full((B,), float("nan"))
    raf_mass = torch.zeros(B)
    for b in range(B):
        if l1["RAF_members"][b] is None:
            continue
        si = l1["RAF_members"][b][0]
        row = b * N + si
        c_, g_ = float(cat_row[row]), float(grow_row[row])
        raf_cat_share[b] = c_ / max(c_ + g_, 1.0)
        specs = l1["RAF_product_specs"][b]
        if specs:
            occ_r = L[row] > 0
            hit = torch.zeros_like(occ_r)
            for (pc, pl, ps) in specs:
                hit |= occ_r & (chem[row] == pc) & (L[row] >= pl) & ((seq[row] & cat4c2._mask(torch.tensor(pl * BITS, device=dev))) == ps)
            raf_mass[b] = float(L[row][hit].sum())
    out_probe = None
    if probe is not None:
        out_probe = torch.zeros(B)
        for b, pr in enumerate(probe):
            if pr is None:
                continue
            row = b * N + pr[0]
            k = cat4c2.species_key(chem[row], L[row], seq[row])[L[row] > 0].cpu()
            want = torch.tensor(sorted(pr[1]), dtype=torch.int64)
            out_probe[b] = float(torch.isin(k, want).sum()) if want.numel() else 0.0
    out.update(l1)
    out["L1_loose"] = l1["L1_candidate"]
    out["RAF_site_catalysed_share"] = raf_cat_share
    out["RAF_product_mass"] = raf_mass
    out["probe_abundance"] = out_probe
    out["cat_events"] = cat_events
    out["cat_thinned_triples_batch"] = sat_triples
    out["polymer_catalysed_share"] = cat_events / torch.clamp(cat_events + cond_events, min=1.0)
    out["mineral_factor"] = f_min
    out["t_first_E2_raw"] = t_first_e2
    out["L1_minus_E2_time"] = l1["t_first_L1"] - t_first_e2
    out.update(_l2_from_record(L, seq, chem, pair, l2_rows, l2_keys, l2_t, next_chk2, t, tab, kuv, w, B, n2,
                               e1 & liquid, p))
    out["L2_loose"] = out["L2_candidate"]
    out["L2_order_vs_L1"] = torch.sign(l1["t_first_L1"] - out["t_first_L2"])   # +1: L2 first, -1: L1 first
    out["copies"] = cp_n
    out["copy_mean_len"] = cp_len / torch.clamp(cp_n, min=1.0)
    out["copy_error_rate"] = cp_err / torch.clamp(cp_ext, min=1.0)       # errors per copied (non-primer) letter
    out["copy_primed_share"] = cp_primed / torch.clamp(cp_n, min=1.0)
    out["copy_mean_added"] = cp_ext / torch.clamp(cp_n, min=1.0)       # letters added by copying per copy
    out["copy_tries"] = cp_try
    out["copies_long"] = cp_long
    # long copies (>= L_H) per template-capable chain per L_H-mer lifetime: a reproduction number
    out["templates_avg_per_site"] = tmpl_time / torch.clamp(t, min=1e-300)
    out["long_copies_per_template_per_lifetime"] = (cp_long / N) / torch.clamp(tmpl_time, min=1e-300) * out["L2_tau_h"]
    out["copy_r_mean"] = cp_r / torch.clamp(cp_n, min=1.0)                # mean extension rate r at copies (/h)
    out["dup_e"] = tab["dup_e"]                                           # [B,2] duplex exponent, wet / dry
    out["copy_tries_q4"] = cp_try_q4
    out["copies_q4"] = cp_q4
    out["copies_dropped"] = cp_drop
    out["round_cut_site_steps"] = rc_sites
    out["round_cut_copies_lost"] = rc_lost
    if COPY_CONC_BINS:                      # Review 9 record: Gini and top-10% share of copies over occupied lineage bins
        cb_ = cc_bins.view(B, COPY_CONC_BINS).cpu()
        g_, t_, n_ = torch.full((B,), float("nan"), dtype=torch.float64), torch.full((B,), float("nan"), dtype=torch.float64), torch.zeros(B, dtype=torch.float64)
        for b_ in range(B):
            x_ = torch.sort(cb_[b_][cb_[b_] > 0]).values
            n_[b_] = x_.numel()
            if x_.numel():
                i_ = torch.arange(1, x_.numel() + 1, dtype=torch.float64)
                g_[b_] = float((2 * (i_ * x_).sum() / (x_.numel() * x_.sum())) - (x_.numel() + 1) / x_.numel())
                k_ = max(1, int(math.ceil(0.1 * x_.numel())))
                t_[b_] = float(x_[-k_:].sum() / x_.sum())
        out["copy_conc_gini"], out["copy_conc_top10"], out["copy_conc_bins"] = g_, t_, n_
    if MULTI_START:     # P3-1b records (no effect on the dynamics)
        ts_ = torch.clamp(ms_rec["tsteps"], min=1.0)
        out["ms_template_steps"] = ms_rec["tsteps"]
        out["ms_occupancy_capped_share"] = ms_rec["occ"] / ts_
        for k_ in ("slots", "blocks", "primers", "any"):
            out[f"ms_limited_{k_}_share"] = ms_rec[k_] / ts_
    out["copy_saturated_template_steps"] = cp_sat
    lt = f32(p["L_t"]).long()[run_of_row][:, None]
    tm = (chem == 0) & (L >= lt)
    out["paired_template_share"] = ((tm & (pair >= 0)).sum(1).view(B, N).sum(1).double().cpu()
                                    / torch.clamp(tm.sum(1).view(B, N).sum(1).double().cpu(), min=1.0))
    out["copy_log"] = copy_log
    out.update(VS.outputs(t, t_end, l2_t, next_chk2, L, chem, lab, out) if ves_on else _e3_empty(B))
    if special.get("return_state"):
        out["state"] = {"L": L.cpu(), "seq": seq.cpu(), "chem": chem.cpu(), "N": N}
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
    classes = [(2, 3), (4, 7), (8, 15), (16, 31), (32, 64)]
    distinct = torch.zeros(B, len(classes))
    entropy = torch.zeros(B, len(classes))
    for ci, (lo, hi) in enumerate(classes):
        m = occ & (L >= lo) & (L <= hi)
        if not bool(m.any()):
            continue
        key = cat4c2.species_key(chem[m], L[m], seq[m])
        pair = torch.stack([runs[m], key], 1)
        uk, cnt = torch.unique(pair, dim=0, return_counts=True)
        r = uk[:, 0].cpu()                        # sums on CPU in float64: deterministic order
        cnt = cnt.cpu().double()
        distinct[:, ci] = torch.zeros(B, dtype=torch.float64).index_add_(0, r, torch.ones_like(cnt)).float()
        tot = torch.zeros(B, dtype=torch.float64).index_add_(0, r, cnt)
        pr = cnt / tot[r]
        entropy[:, ci] = torch.zeros(B, dtype=torch.float64).index_add_(0, r, -pr * torch.log2(pr)).float()
    return {"length_hist": hist.cpu(), "chain_mass_by_chem": mass_c.cpu(), "pools_F_M_Mstar_by_chem": pools.cpu(),
            "distinct_seq_by_class": distinct, "entropy_bits_by_class": entropy,
            "length_classes": classes}


def _l1_from_history(raf_hist, p, tab, w, kuv, B, gate):
    """L1 candidates from the checkpoint history. At every checkpoint of a sampled site the
    maxRAF of its abundant-species network is computed with that checkpoint's chemistry (hash
    salt). The persistent RAF of a window of consecutive checkpoints ending at the last one is
    the maxRAF of the catalysed reactions (reactants, product AND catalyst) present in the RAF at
    every checkpoint of the window, so a reshuffled chemistry cannot persist. Candidate: a
    persistent RAF >= 3 reactions over a window spanning >= 5 member turnover times (mean over
    the RAF's products of 1 / ((L - 1)(k_h + k_UV) + 1/tau)); the longest window is kept. The
    self-sustaining test needs the paired run without catalysis and is applied by the caller."""
    raf_any = torch.zeros(B, dtype=torch.bool)
    persist = torch.zeros(B, dtype=torch.bool)
    size = torch.zeros(B)
    comps = torch.zeros(B)
    span_turnovers = torch.zeros(B, dtype=torch.float64)
    t_first = torch.full((B,), float("nan"), dtype=torch.float64)
    members = [None] * B
    cat_L = [[] for _ in range(B)]
    nuc = torch.full((B,), float("nan"))
    specs_out = [None] * B
    for b in range(B):
        if not bool(gate[b]) or not raf_hist[b] or not raf_hist[b][0]:
            continue
        pp = (float(p["p0"][b]), int(p["L_cat"][b]), float(p["lam"][b]))
        best = None
        for si, hist in enumerate(raf_hist[b]):
            res = [cat4c2.site_raf(Xh, *pp, salt_h, ALPHA, BITS) for (_, Xh, _, salt_h) in hist]
            if res[-1].get("size", 0) >= 3:
                raf_any[b] = True
            if res[-1].get("size", 0) < 3:
                continue
            X_last = hist[-1][1]
            t_last = hist[-1][0]
            info = {}
            for x in range(len(hist)):
                for kk, cnt in res[x].get("member_counts", {}).items():
                    info[kk] = cnt
            Lk = {int(kk): (int(l), int(c)) for kk, l, c, v in zip(X_last["key"].tolist(), X_last["L"].tolist(),
                                                                  X_last["chem"].tolist(), X_last["valid"].tolist()) if v}
            tup = set(res[-1].get("cat_tuples", set()))
            for k in range(len(hist) - 2, -1, -1):
                tup &= res[k].get("cat_tuples", set())
                n_raf, kept = cat4c2.tuples_raf(tup, ALPHA)
                if n_raf < 3:
                    break
                prods = {t[2] for t in kept}
                k_member = []
                for kk in prods:
                    if kk in Lk:
                        Lm, cm = Lk[kk]
                        kh = float(tab["knuc"][b, 2] if cm == 0 else tab["kpep"][b, 2])
                        k_member.append(float(turnover_rate(Lm, kh, kuv[b], w[b], p["_theta_eff"][b], p["p_min"][b])))
                if not k_member:
                    continue
                turnover = 1.0 / max(sum(k_member) / len(k_member), 1e-300)
                span = (t_last - hist[k][0]) / turnover
                if best is None or span > best[0]:
                    counts = {kk: res[-1].get("member_counts", {}).get(kk, 0) for kk in prods}
                    best = (span, si, n_raf, kept, counts, hist[k][0], res[-1])
        if best is not None:
            span, si, n_raf, kept, counts, t0, rl = best
            size[b] = n_raf
            comps[b] = rl.get("components", 0)
            span_turnovers[b] = span
            persist[b] = span >= 5
            members[b] = (si, counts)               # products of the persistent RAF (self-sustaining test)
            cat_L[b] = rl.get("catalyst_L", [])
            nuc[b] = rl.get("nuc_share", float("nan"))
            t_first[b] = t0
            X_last = raf_hist[b][si][-1][1]
            kx = X_last["key"].tolist()
            specs_out[b] = [(int(X_last["chem"][i]), int(X_last["L"][i]), int(X_last["seq"][i]))
                            for i, kk in enumerate(kx) if kk in counts and bool(X_last["valid"][i])]
    return {"RAF_final": raf_any, "L1_candidate": persist, "RAF_size": size,
            "RAF_components": comps, "RAF_span_turnovers": span_turnovers, "t_first_L1": t_first,
            "RAF_members": members, "RAF_catalyst_L": cat_L, "RAF_nuc_share": nuc,
            "RAF_product_specs": specs_out}


def _l2_from_record(L, seq, chem, pair, l2_rows, l2_keys, l2_t, next_chk2, t_fin, tab, kuv, w, B, ns, gate, p):
    """L2 candidates from the lineage record (see module docstring). Returns per-run tensors."""
    dev = L.device
    rows = l2_rows.reshape(-1)
    sh = _l2_shares(L, seq, chem, rows, l2_keys.view(B * ns, N_CHECK, -1)).view(B, ns, N_CHECK).cpu().double()
    tau_h = 1.0 / turnover_rate(L_H, tab["knuc"][:, 2], kuv, w, p["_theta_eff"], p["p_min"])   # [B]
    span = (t_fin[:, None] - l2_t) / tau_h[:, None]                               # [B, N_CHECK]
    rec = torch.arange(N_CHECK)[None, :] < (next_chk2 - 1)[:, None]               # checkpoint recorded
    good = (sh >= L2_SHARE) & rec[:, None, :]
    sp = torch.where(good, span[:, None, :].expand(B, ns, N_CHECK), torch.full_like(sh, -1.0))
    best_k = sp.argmax(2)                                                         # earliest k = longest span
    best_sp = sp.max(2).values                                                    # [B, ns]
    site = best_sp.argmax(1)
    b_ix = torch.arange(B)
    k_star = best_k[b_ix, site]
    span_star = best_sp[b_ix, site].clamp(min=0.0)
    share_star = sh[b_ix, site, k_star]
    cand = gate & (best_sp[b_ix, site] >= L2_SPAN)
    # share over the last window of >= L2_SPAN turnover times (latest such checkpoint), best site
    ok5 = (span >= L2_SPAN) & rec
    k5 = torch.where(ok5, torch.arange(N_CHECK)[None, :].expand(B, N_CHECK), torch.full((B, N_CHECK), -1)).max(1).values
    share5 = torch.where(k5 >= 0, sh[b_ix, :, k5.clamp(min=0)].max(1).values, torch.zeros(B, dtype=torch.float64))
    t_first = torch.where(cand, l2_t[b_ix, k_star], torch.full((B,), float("nan"), dtype=torch.float64))
    # length of the most abundant lineage (chemistry 0, L >= L_H) in the best site at the end
    best_rows = l2_rows[b_ix.to(dev), site.to(dev)]
    Lb, sb, cb = L[best_rows], seq[best_rows], chem[best_rows]
    X = cat4c2.abundant(torch.where((cb == 0) & (Lb >= L_H), Lb, torch.zeros_like(Lb)), sb, cb, nx=1)
    top_len = torch.where(X["valid"][:, 0], X["L"][:, 0], torch.zeros_like(X["L"][:, 0])).cpu().double()
    return {"L2_candidate": cand, "L2_span_turnovers": span_star, "L2_share": share_star, "L2_site": site,
            "L2_check": k_star, "L2_share_last5": share5, "L2_share_matrix": sh, "t_first_L2": t_first,
            "L2_top_lineage_len": top_len, "L2_tau_h": tau_h}


def paired_control_l2(u, seed, fidelity, out, override=None, record=False, both=False):
    """L2_strict for runs with L2_loose. The candidate runs are rerun as a batch of their own:
    once as they are (main', seed `seed`) and with copy_random (c1, seed `seed`: identical random
    numbers, so main' and c1 differ in the copied letters only); for the runs with S(main') -
    S(c1) >= L2_MARGIN also with copy_random and seed `seed` + 1 (c2, an independent
    realisation; `both`: for every candidate). The statistic is the same in all reruns, so no
    side is selected: S = L2_share_last5, the largest information share over all sites at the
    latest checkpoint whose window spans >= L2_SPAN turnover times. L2_strict: L2_loose in the
    original run and in main', and S(main') - max(S(c1), S(c2)) >= L2_MARGIN (coordinator
    decision 26 Oct 2026, after a single control passed 3/16 runs with random letters at the
    box reference). Running c2 only where c1 alone would pass gives the same verdicts.
    Returns bool [B], S(main') and max(S(c1), S(c2)) (S(c1) where c2 was not run)."""
    B = u.shape[0]
    strict = torch.zeros(B, dtype=torch.bool)
    s_main = torch.full((B,), float("nan"), dtype=torch.float64)
    s_ctl = torch.full((B,), float("nan"), dtype=torch.float64)
    idx = out["L2_loose"].nonzero().squeeze(1).tolist()
    if not idx:
        return strict, s_main, s_ctl

    def rerun(ix, sd, rnd):
        g = torch.Generator(device=u.device)
        g.manual_seed(sd)
        ov = dict(override or {})
        if rnd:
            ov["copy_random"] = True
        return simulate(u[ix], g, fidelity, record=record, override=ov)
    m = rerun(idx, seed, False)
    c1 = rerun(idx, seed, True)
    sm, sc = m["L2_share_last5"].clone(), c1["L2_share_last5"].clone()
    need = [k for k in range(len(idx)) if both or (sm[k] - sc[k] >= L2_MARGIN)]
    if need:
        c2 = rerun([idx[k] for k in need], seed + 1, True)
        nk = torch.tensor(need)
        sc[nk] = torch.maximum(sc[nk], c2["L2_share_last5"])
    ii = torch.tensor(idx)
    s_main[ii] = sm
    s_ctl[ii] = sc
    strict[ii] = m["L2_candidate"] & (sm - sc >= L2_MARGIN)
    return strict, s_main, s_ctl


def paired_control(u, seed, fidelity, out, override=None, record=False):
    """L1_strict: for runs with L1_loose, rerun the same conditions and seed with p0 = 0 and
    require the RAF products (best site) to be more abundant with catalysis. Returns a bool
    tensor [B] and the two abundance vectors."""
    B = u.shape[0]
    strict = torch.zeros(B, dtype=torch.bool)
    ab_cat = torch.zeros(B)
    ab_ctl = torch.zeros(B)
    idx = out["L1_loose"].nonzero().squeeze(1).tolist()
    if not idx:
        return strict, ab_cat, ab_ctl
    probe = []
    for b in idx:
        si, counts = out["RAF_members"][b]
        probe.append((si, set(counts.keys())))
        ab_cat[b] = float(sum(counts.values()))
    g = torch.Generator(device=u.device)
    g.manual_seed(seed)
    ov = dict(override or {})
    ov["p0"] = 0.0
    c = simulate(u[idx], g, fidelity, record=record, override=ov, probe=probe)
    ab_ctl[idx] = c["probe_abundance"]
    strict = out["L1_loose"] & (ab_cat > ab_ctl)
    return strict, ab_cat, ab_ctl


# ---------------------------------------------------------------- 4c-4 vesicle layer ----------

E3_KEYS_F = ("E3_frac_sites_persistent", "E3_span_turnovers", "E3_frac_sites_enclosing", "t_first_E3", "E3_order_vs_L2",
             "vesicles_per_site", "membrane_share", "mean_vesicle_size", "enclosed_chain_share", "enclosed_L2_share",
             "divisions_per_lifetime", "max_generation", "vesicles_formed", "divisions", "free_amphiphile_per_site")


def _e3_empty(B):
    out = {"E3_loose": torch.zeros(B, dtype=torch.bool), "E3_strict": torch.zeros(B, dtype=torch.bool)}
    for k in E3_KEYS_F:
        out[k] = torch.zeros(B, dtype=torch.float64) if k not in ("t_first_E3", "E3_order_vs_L2") else \
            torch.full((B,), float("nan"), dtype=torch.float64)
    return out


class Vesicles:
    """State and per-step dynamics of the vesicle layer (module docstring). Amounts of amphiphile are
    continuous (molecules, deterministic expected values); vesicle and chain events are random,
    drawn from the layer's own generator so the chemistry's random numbers are untouched."""

    def __init__(self, p, tab, T_ph, E, in_win, B, N, R, S, dev, run_of_row, special, generator, record, cycling):
        self.B, self.N, self.R, self.S, self.dev, self.rr = B, N, R, S, dev, run_of_row
        f = lambda x: torch.as_tensor(x, dtype=torch.float64).to(dev, torch.float32)[run_of_row]
        self.special, self.record = special, record
        closed = bool(special.get("amph_closed", False))
        self.closed = closed
        amph_s = p["feed"] * p["amph"]
        if "amph_total" in special or closed:
            amph_s = torch.zeros_like(amph_s)
        self.s = f(amph_s)
        self.E = f(E)
        w = 1 / p["tau"]
        self.w = f(torch.zeros_like(w) if closed else w)
        self.mix = f(p["mix"])
        d = p["depth"]
        cev = torch.stack([torch.ones_like(d), p["evap"] ** d], 1)
        cev = torch.cat([cev, cev.mean(1, keepdim=True)], 1)
        Tt = torch.cat([T_ph, T_ph.mean(1, keepdim=True)], 1)
        self.cvc_tab = (p["CVC0"][:, None] * torch.exp(p["gamma_T"][:, None] * (Tt - 25.0)) * V / cev)   # [B,3] free A
        self.win_run = in_win.clone()
        self.close_t = float(special.get("window_close_t", float("inf")))
        sw = torch.where((p["depth"] > 0), 2 / p["period"], torch.zeros_like(p["period"]))
        kdiv = p["s_shear"] * p["mix"] + p["s_cycle"] * sw
        if special.get("no_division"):
            kdiv = torch.zeros_like(kdiv)
        self.kdiv, self.phi0, self.frehyd, self.beta = f(kdiv), f(p["phi0"]), f(p["f_rehyd"]), f(p["beta"])
        self.loss_run = (torch.zeros_like(w) if closed else w) + K_BREAK * (~in_win).double()
        # amphiphile pools per site: precursor FA, free A, non-vesicle membrane Bf
        if "amph_total" in special:
            tot = torch.as_tensor(special["amph_total"], dtype=torch.float64).expand(B).clone()
            FA0 = torch.zeros(B, dtype=torch.float64)
        else:
            Q = _linear_q(E, torch.zeros_like(E), torch.zeros_like(E), 1 / p["tau"])
            _, J = _expm_pair(Q, p["avail_time"])
            n = (p["feed"] * p["amph"] * p["avail_time"])[:, None] * J[:, 0, :2]
            FA0, tot = torch.clamp(n[:, 0], max=POOL_MAX), torch.clamp(n[:, 1], max=POOL_MAX)
        self.FA, self.A, self.Bf = f(FA0), f(tot), torch.zeros(R, device=dev)
        self.va = torch.zeros(R, NV, device=dev)
        self.alive = torch.zeros(R, NV, dtype=torch.bool, device=dev)
        self.vid = torch.full((R, NV), -1, dtype=torch.int64, device=dev)
        self.vpar = torch.full((R, NV), -1, dtype=torch.int64, device=dev)
        self.vgen = torch.zeros(R, NV, dtype=torch.int64, device=dev)
        self.vnext = torch.zeros(R, dtype=torch.int64, device=dev)
        self.n_form = torch.zeros(R, device=dev)
        self.n_div = torch.zeros(R, device=dev)
        self.ves_time = torch.zeros(R, device=dev)
        self.max_gen = torch.zeros(R, dtype=torch.int64, device=dev)
        self.cycling = cycling
        self.period_r = f(p["period"])
        self.wd_r = (p["depth"] > 0).to(dev)[run_of_row]        # wet-dry cycling present (depth > 0)
        self.prev_ph = None
        self.g = torch.Generator(device=dev)
        self.g.manual_seed(generator.initial_seed() + 104729)
        self.e3_n = torch.zeros(B, N, N_CHECK, dtype=torch.int16)
        self.e3_enc = torch.zeros(B, N, N_CHECK, dtype=torch.bool)
        self.log = {"formation": [], "division": []} if special.get("ves_log") else None
        self.rows = torch.arange(R, device=dev)

    # ---- helpers ----
    def _inside(self, lab, vmask):
        """bool [R,S]: chains inside the vesicles in vmask [R,NV]."""
        j = (lab.long() - 1).clamp(min=0)
        return (lab > 0) & vmask.gather(1, j)

    def _u_pair(self, u, pair):
        """Draws per slot shared by the two strands of a duplex (the lower slot's draw)."""
        S = self.S
        lo = torch.where(pair >= 0, torch.minimum(pair, torch.arange(S, device=self.dev)[None, :]),
                         torch.arange(S, device=self.dev)[None, :].expand_as(pair))
        return torch.gather(u, 1, lo)

    def _letters_in(self, L, lab):
        m = lab > 0
        idx = (self.rows[:, None] * NV + lab.long() - 1)[m]
        return torch.zeros(self.R * NV, device=self.dev).index_add_(0, idx, L[m].float()).view(self.R, NV)

    def init_chains(self, L, seq, chem, pair, dlen, lab, special):
        S = self.S
        it = special.get("init_duplexes")
        if it:
            j0 = int((L[0] > 0).sum())
            Ld, sd = int(it["L"]), int(it["seq"])
            for k in range(int(it["n"])):
                a, b = j0 + 2 * k, j0 + 2 * k + 1
                L[:, a] = Ld; L[:, b] = Ld
                seq[:, a] = sd; seq[:, b] = (~sd) & ((1 << Ld) - 1)
                chem[:, a] = 0; chem[:, b] = 0
                pair[:, a] = b; pair[:, b] = a
                dlen[:, a] = Ld; dlen[:, b] = Ld
        iv = special.get("init_vesicles")
        if iv:
            j0 = int((L[0] > 0).sum())
            Lv, sv = int(iv["L"]), int(iv["seq"])
            prs = iv.get("pairs", [0] * len(iv["a"]))
            for j, (arel, nch, npr) in enumerate(zip(iv["a"], iv["n"], prs)):
                self.va[:, j] = float(arel) * A0_VES
                self.alive[:, j] = True
                self.vid[:, j] = j
                for c in range(int(nch)):
                    L[:, j0] = Lv; seq[:, j0] = sv; chem[:, j0] = 0; lab[:, j0] = j + 1
                    j0 += 1
                for c in range(int(npr)):                      # duplexes inside vesicle j
                    a_, b_ = j0, j0 + 1
                    L[:, a_] = Lv; L[:, b_] = Lv; seq[:, a_] = sv; seq[:, b_] = (~sv) & ((1 << Lv) - 1)
                    chem[:, a_] = 0; chem[:, b_] = 0; lab[:, a_] = j + 1; lab[:, b_] = j + 1
                    pair[:, a_] = b_; pair[:, b_] = a_; dlen[:, a_] = Lv; dlen[:, b_] = Lv
                    j0 += 2
            self.vnext[:] = len(iv["a"])

    def _kill(self, vmask):
        self.alive &= ~vmask
        self.va = torch.where(vmask, torch.zeros_like(self.va), self.va)

    # ---- one step ----
    def step(self, t, dt, ph, L, seq, chem, pair, dlen, lab, unpair, Ms):
        dev, R, S = self.dev, self.R, self.S
        dtr = dt.to(dev, torch.float32)[self.rr]
        phr = ph.to(dev)[self.rr]
        g = self.g
        # 1. amphiphile pools (exact for the precursor, first order for A)
        k1 = self.E + self.w
        e1 = torch.exp(-k1 * dtr)
        k1c = k1.clamp(min=1e-30)
        FAss = self.s / k1c
        FAint = torch.where(k1 > 0, FAss * dtr + (self.FA - FAss) * (1 - e1) / k1c, self.FA * dtr)
        self.FA = torch.where(k1 > 0, FAss + (self.FA - FAss) * e1, self.FA + self.s * dtr)
        ew = torch.exp(-self.w * dtr)
        self.A = self.A * ew + self.E * FAint
        self.Bf = self.Bf * ew
        m = 1 - torch.exp(-self.mix * dtr)
        Am = self.A.view(self.B, self.N).mean(1, keepdim=True).expand(self.B, self.N).reshape(R)
        self.A = self.A + (Am - self.A) * m
        # 2. assembly around the CVC
        Aeq = self.cvc_tab.gather(1, ph[:, None]).squeeze(1).to(dev, torch.float32)[self.rr]
        fa = 1 - torch.exp(-K_AS * dtr)
        exc = self.A - Aeq
        mv = exc.clamp(min=0) * fa
        need = (-exc).clamp(min=0) * fa
        tb = torch.minimum(need, self.Bf)
        rem = need - tb
        totv = (self.va * self.alive).sum(1)
        take = torch.minimum(rem, totv)
        frac = take / totv.clamp(min=1e-30)
        self.va = self.va * (1 - frac[:, None])
        self.A = self.A - mv + tb + take
        self.Bf = self.Bf + mv - tb
        dis = self.alive & (self.va < 0.5 * A0_VES)
        if bool(dis.any()):
            lab.masked_fill_(self._inside(lab, dis), 0)
            self.A = self.A + (self.va * dis).sum(1)
            self._kill(dis)
        win = (self.win_run & (t < self.close_t)).to(dev)[self.rr]
        # 3. division (before washout so duplex bookkeeping uses the current pairs)
        elig = self.alive & (self.va >= 2 * A0_VES) & win[:, None]
        if bool(elig.any()) and bool((self.kdiv > 0).any()):
            ud = torch.rand(R, NV, device=dev, generator=g)
            div = elig & (ud < (1 - torch.exp(-self.kdiv * dtr))[:, None])
            if bool(div.any()):
                self._divide(t, div, L, lab, pair)
        # 4. rehydration: at the end of a dry phase (resolved cycling), or as an averaged effect of
        #    n = dt / period cycles in a cycle-averaged step with wet-dry cycling (chem4c4 v2):
        #    p_j -> 1 - (1 - p_j)^n, sequential capping as for one event
        reh_res = (self.cycling & (self.prev_ph == 1) & (ph == 0)).to(dev)[self.rr] if self.prev_ph is not None \
            else torch.zeros(R, dtype=torch.bool, device=dev)
        reh_avg = (phr == 2) & self.wd_r
        reh = (reh_res | reh_avg) & win & self.alive.any(1)
        if bool(reh.any()):
            if True:
                phi = (self.phi0[:, None] * (self.va / A0_VES) ** 1.5).clamp(max=1.0) * self.alive
                pj = (phi * self.frehyd[:, None]).clamp(max=1.0)
                ncyc = (dtr / self.period_r.clamp(min=1e-30))[:, None]
                pj = torch.where(reh_avg[:, None] & ~reh_res[:, None], 1 - (1 - pj) ** ncyc, pj)
                cum = torch.cumsum(pj, 1).clamp(max=1.0)
                u = self._u_pair(torch.rand(R, S, device=dev, generator=g), pair)
                outc = (L > 0) & (lab == 0) & reh[:, None]
                lo = torch.zeros_like(u)
                for j in range(NV):
                    hit = outc & (u >= lo) & (u < cum[:, j:j + 1]) & self.alive[:, j:j + 1]
                    lab.masked_fill_(hit, j + 1)
                    outc = outc & ~hit
                    lo = torch.maximum(lo, cum[:, j:j + 1].expand_as(lo))
        # 5. formation
        if not self.special.get("no_formation"):
            nf = NV - self.alive.sum(1)
            lamf = (K_FORM * self.Bf / A0_VES * dtr) * win * (nf > 0)
            if bool((lamf > 0).any()):
                nnew = torch.poisson(lamf.contiguous(), generator=g)
                nnew = torch.minimum(torch.minimum(nnew, nf.float()), torch.floor(self.Bf / A0_VES)).long()
                if bool((nnew > 0).any()):
                    self._form(t, nnew, L, lab, pair)
        # 6. growth from non-vesicle membrane and exchange between vesicles
        anyv = self.alive.any(1) & win
        if bool(anyv.any()):
            ench = self._letters_in(L, lab)
            outl = torch.where((lab == 0) & (L > 0), L, torch.zeros_like(L)).sum(1).float()
            phi = (self.phi0[:, None] * (self.va / A0_VES) ** 1.5).clamp(max=1.0) * self.alive
            outvol = (1 - phi.sum(1)).clamp(min=0.05)
            Pi = ench / phi.clamp(min=1e-30) / (outl / outvol + 1.0)[:, None]
            wj = self.va * (1 + self.beta[:, None] * Pi) * self.alive
            sw = wj.sum(1).clamp(min=1e-30)
            share = wj / sw[:, None]
            gro = self.Bf * (1 - torch.exp(-K_GROW * dtr)) * anyv
            self.va = self.va + gro[:, None] * share
            self.Bf = self.Bf - gro
            tot = (self.va * self.alive).sum(1)
            fx = (1 - torch.exp(-K_EX * dtr)) * anyv
            self.va = torch.where(self.alive, self.va + (share * tot[:, None] - self.va) * fx[:, None], self.va)
        # 7. breakdown outside the window (contents released, membrane to non-vesicle membrane)
        out_w = self.alive & ~win[:, None]
        if bool(out_w.any()):
            ub = torch.rand(R, NV, device=dev, generator=g)
            brk = out_w & (ub < (1 - torch.exp(-K_BREAK * dtr))[:, None])
            if bool(brk.any()):
                lab.masked_fill_(self._inside(lab, brk), 0)
                self.Bf = self.Bf + (self.va * brk).sum(1)
                self._kill(brk)
        # 8. washout with contents
        if not self.closed and bool(self.alive.any()):
            uw = torch.rand(R, NV, device=dev, generator=g)
            wsh = self.alive & (uw < (1 - torch.exp(-self.w * dtr))[:, None])
            if bool(wsh.any()):
                cm = self._inside(lab, wsh)
                unpair(cm)
                L.masked_fill_(cm, 0)
                seq.masked_fill_(cm, 0)
                lab.masked_fill_(cm, 0)
                self._kill(wsh)
        self.ves_time = self.ves_time + self.alive.sum(1) * dtr
        self.prev_ph = ph.clone()
        self.L, self.seq = L, seq

    def _form(self, t, nnew, L, lab, pair):
        dev, R, S = self.dev, self.R, self.S
        free = ~self.alive
        rank = torch.cumsum(free.long(), 1) - 1
        new = free & (rank < nnew[:, None])
        self.alive |= new
        self.va = torch.where(new, torch.full_like(self.va, A0_VES), self.va)
        ids = self.vnext[:, None] + rank
        self.vid = torch.where(new, ids, self.vid)
        self.vpar = torch.where(new, torch.full_like(self.vpar, -1), self.vpar)
        self.vgen = torch.where(new, torch.zeros_like(self.vgen), self.vgen)
        self.vnext = self.vnext + nnew
        self.Bf = self.Bf - nnew.float() * A0_VES
        self.n_form = self.n_form + nnew.float()
        # slot of the k-th new vesicle in each row
        slot_by_rank = torch.full((R, NV), -1, dtype=torch.long, device=dev)
        rr_, cc_ = new.nonzero(as_tuple=True)
        slot_by_rank[rr_, rank[rr_, cc_]] = cc_
        width = torch.minimum(self.phi0, 1.0 / nnew.clamp(min=1).float())
        u = self._u_pair(torch.rand(R, S, device=dev, generator=self.g), pair)
        k = torch.floor(u / width[:, None]).long()
        outc = (L > 0) & (lab == 0) & (nnew[:, None] > 0)
        enc = outc & (k < nnew[:, None])
        tgt = torch.gather(slot_by_rank, 1, k.clamp(max=NV - 1))
        lab.copy_(torch.where(enc, (tgt + 1).to(torch.int8), lab))
        if self.log is not None:
            n_out = outc.sum(1)
            for r in (nnew > 0).nonzero().squeeze(1).tolist():
                cnt = [int((enc[r] & (k[r] == q)).sum()) for q in range(int(nnew[r]))]
                self.log["formation"].append({"row": r, "n_out": int(n_out[r]), "enclosed": cnt,
                                              "phi": float(width[r])})

    def _divide(self, t, div, L, lab, pair):
        dev, R, S = self.dev, self.R, self.S
        free = ~self.alive
        nfree = free.sum(1)
        rf = torch.cumsum(free.long(), 1) - 1
        slot_by_rank = torch.full((R, NV), -1, dtype=torch.long, device=dev)
        rr_, cc_ = free.nonzero(as_tuple=True)
        slot_by_rank[rr_, rf[rr_, cc_]] = cc_
        rdiv = torch.cumsum(div.long(), 1) - 1
        has = div & (rdiv < nfree[:, None])
        d2 = torch.where(has, torch.gather(slot_by_rank, 1, rdiv.clamp(min=0, max=NV - 1)), torch.full_like(rdiv, -1))
        a_m = self.va.clone()
        half = a_m / 2
        # chains of dividing vesicles: half to the second daughter (released if it has no slot)
        j = (lab.long() - 1).clamp(min=0)
        inm = (lab > 0) & div.gather(1, j)
        u = self._u_pair(torch.rand(R, S, device=dev, generator=self.g), pair)
        mv = inm & (u < 0.5)
        newlab = torch.gather(d2, 1, j) + 1                                # 0 if released
        if self.log is not None:
            n_m = torch.zeros(R * NV, device=dev).index_add_(0, (self.rows[:, None] * NV + j)[inm], torch.ones(int(inm.sum()), device=dev)).view(R, NV)
            n_2 = torch.zeros(R * NV, device=dev).index_add_(0, (self.rows[:, None] * NV + j)[mv], torch.ones(int(mv.sum()), device=dev)).view(R, NV)
        lab.copy_(torch.where(mv, newlab.to(torch.int8), lab))
        # membrane: halves (second half to Bf if no slot)
        self.va = torch.where(div, half, self.va)
        d2m = torch.zeros_like(self.va).scatter_add_(1, d2.clamp(min=0), torch.where(has, half, torch.zeros_like(half)))
        born = (torch.zeros(R, NV, dtype=torch.long, device=dev).scatter_add_(1, d2.clamp(min=0), has.long()) > 0) & ~self.alive
        self.va = torch.where(born, d2m, self.va)
        self.Bf = self.Bf + torch.where(div & ~has, half, torch.zeros_like(half)).sum(1)
        # ids: first daughter keeps the slot with a new id, second daughter gets a free slot
        old_id, old_gen = self.vid.clone(), self.vgen.clone()
        nd = div.sum(1)
        id1 = self.vnext[:, None] + rdiv
        self.vid = torch.where(div, id1, self.vid)
        self.vpar = torch.where(div, old_id, self.vpar)
        self.vgen = torch.where(div, old_gen + 1, self.vgen)
        self.vnext = self.vnext + nd
        r2 = torch.cumsum(has.long(), 1) - 1
        par2 = torch.zeros_like(self.vid).scatter_add_(1, d2.clamp(min=0), torch.where(has, old_id, torch.zeros_like(old_id)))
        gen2 = torch.zeros_like(self.vgen).scatter_add_(1, d2.clamp(min=0), torch.where(has, old_gen + 1, torch.zeros_like(old_gen)))
        id2 = torch.zeros_like(self.vid).scatter_add_(1, d2.clamp(min=0), torch.where(has, self.vnext[:, None] + r2, torch.zeros_like(r2)))
        self.vid = torch.where(born, id2, self.vid)
        self.vpar = torch.where(born, par2, self.vpar)
        self.vgen = torch.where(born, gen2, self.vgen)
        self.vnext = self.vnext + has.sum(1)
        self.alive |= born
        self.n_div = self.n_div + nd.float()
        self.max_gen = torch.maximum(self.max_gen, (self.vgen * self.alive).max(1).values)
        if self.log is not None:
            # duplexes split across daughters (must be 0) and membrane before / after
            pp = pair >= 0
            partner_lab = torch.gather(lab, 1, pair.clamp(min=0))
            split = int((pp & (L > 0) & (lab != partner_lab)).sum())
            for r, jj in div.nonzero().tolist():
                k2 = int(d2[r, jj])
                a_after = float(self.va[r, jj]) + (float(self.va[r, k2]) if k2 >= 0 else float(half[r, jj]))
                self.log["division"].append({"row": r, "n": int(n_m[r, jj]), "n2": int(n_2[r, jj]), "slot2": k2,
                                             "a_before": float(a_m[r, jj]), "a_after": a_after})
            self.log.setdefault("duplex_split", []).append(split)

    # ---- detector records ----
    def checkpoint(self, runs2, kidx, rows2, L, lab, chem):
        n = runs2.numel()
        nv = self.alive[rows2].sum(1).view(n, self.N).to("cpu", torch.int16)
        enc = ((lab[rows2] > 0) & (L[rows2] >= L_H)).any(1).view(n, self.N).cpu()
        rc, kc = runs2.cpu(), kidx.cpu()
        self.e3_n[rc, :, kc] = nv
        self.e3_enc[rc, :, kc] = enc

    def outputs(self, t, t_end, l2_t, next_chk2, L, chem, lab, out):
        B, N = self.B, self.N
        nrec = next_chk2 - 1
        # E3 window (chem4c4 v2, parallel to L2): the longest run of consecutive checkpoints ending at
        # t_end in which >= half the sites hold >= 1 vesicle at every checkpoint; E3_loose if it spans
        # >= 5 vesicle turnovers; evaluable if checkpoint 1 -> t_end spans >= 5 turnovers
        pres = self.e3_n > 0                                              # [B, N, K]
        frac_k = torch.stack([pres[:, :, k:].all(2).float().mean(1) for k in range(N_CHECK)], 1)   # [B, K]
        okk = frac_k >= 0.5                                               # nested windows: monotone in k
        kmin = torch.where(okk.any(1), okk.float().argmax(1), torch.full((B,), N_CHECK - 1, dtype=torch.long))
        has_win = okk.any(1)
        persist = frac_k.gather(1, kmin[:, None]).squeeze(1).double()
        span = torch.where(has_win, (t - l2_t.gather(1, kmin[:, None]).squeeze(1)) * self.loss_run,
                           torch.zeros(B, dtype=torch.float64))
        evaluable = (nrec >= N_CHECK) & ((t - l2_t[:, 0]) * self.loss_run >= 5.0)
        loose = evaluable & has_win & (span >= 5.0)
        encf = self.e3_enc[:, :, N_CHECK - 1].float().mean(1).double()
        strict = loose & (encf >= 0.5)
        has = (self.e3_n > 0).float().mean(1) >= 0.5                      # [B, N_CHECK]
        kf = torch.where(has.any(1), has.float().argmax(1), torch.full((B,), -1))
        t_first = torch.where(kf >= 0, l2_t.gather(1, kf.clamp(min=0)[:, None]).squeeze(1),
                              torch.full((B,), float("nan"), dtype=torch.float64))
        rs = lambda x: x.view(B, N).sum(1).double().cpu()
        Lf = L.float()
        inside = lab > 0
        long0 = (chem == 0) & (L >= L_H)
        tot_amph = rs(self.A + self.Bf + (self.va * self.alive).sum(1))
        nves = self.alive.sum(1).float()
        res = {"E3_loose": loose, "E3_strict": strict, "E3_frac_sites_persistent": persist, "E3_span_turnovers": span,
               "E3_window_first_checkpoint": torch.where(has_win, kmin, torch.full_like(kmin, -1)).double(),
               "E3_frac_sites_enclosing": encf, "t_first_E3": t_first,
               "E3_order_vs_L2": torch.sign(out["t_first_L2"] - t_first),       # +1: E3 first
               "vesicles_per_site": rs(nves) / N,
               "membrane_share": rs((self.va * self.alive).sum(1)) / tot_amph.clamp(min=1e-30),
               "mean_vesicle_size": rs((self.va * self.alive).sum(1)) / rs(nves).clamp(min=1) / A0_VES,
               "enclosed_chain_share": rs(torch.where(inside, Lf, 0.0).sum(1)) / rs(Lf.sum(1)).clamp(min=1),
               "enclosed_L2_share": rs(torch.where(inside & long0, Lf, 0.0).sum(1)) /
               rs(torch.where(long0, Lf, 0.0).sum(1)).clamp(min=1),
               "divisions_per_lifetime": rs(self.n_div) / rs(self.ves_time).clamp(min=1e-30) / self.loss_run.clamp(min=1e-30),
               "max_generation": self.max_gen.view(B, N).max(1).values.double().cpu(),
               "vesicles_formed": rs(self.n_form), "divisions": rs(self.n_div),
               "free_amphiphile_per_site": rs(self.A) / N,
               "membrane_total_per_site": rs(self.Bf + (self.va * self.alive).sum(1)) / N,
               "vesicles_series": self.e3_n.sum(1).double(), "E3_evaluable": evaluable}
        if self.log is not None:
            res["ves_log"] = self.log
        if self.special.get("ves_state"):
            res["ves_state"] = {"va": self.va.cpu(), "alive": self.alive.cpu(), "vid": self.vid.cpu(),
                                "vpar": self.vpar.cpu(), "vgen": self.vgen.cpu(), "lab": lab.cpu(),
                                "A": self.A.cpu(), "Bf": self.Bf.cpu()}
        return res
