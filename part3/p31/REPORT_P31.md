# Part III, P3-1 — implementation and validation report (B360, lot-p3/p31/)

30 October 2026, Toronto time. Spec `docs/part3_p31_spec.md` (b4c441b), cold-gap change `docs/part3_cold_gap.md` §4 (9a9ee5a),
side-product inhibition (owner-approved 17:13). **No map and no long run were started.**

## Code versions
| File | md5 | Role |
|---|---|---|
| chem_p31.py (final) | cfc1e3db797789fadda313d038d0c981 | P3-1 + cold-gap (a)-(d) + inhibition (e) |
| chem_p31_pre_inhibition.py | 28d95aaa1fe49ddc95f7a05ce6fc67b4 | cold-gap, no inhibition |
| chem_p31_8b6da3c1.py | 8b6da3c1f1a005d2ff9fcea3d5ad8927 | same code as 28d95aaa minus header lines (tests 2, 4 run on it) |
| chem_p31_pre_coldgap.py | 52234744e10dc165e4d3f9b615989f90 | original P3-1 (first full validation) |
| chem4c3.py (reference copy) | cab66933cc8fcd93f5c5b3c78f780952 | Stage 4c-3, unchanged |
Full md5 list: publish/MD5SUMS.txt.

## Validation results, final code (cfc1e3db), validate_p31_final.json
| Test | Result | Numbers |
|---|---|---|
| 1 duplex thermodynamics | PASS | code vs NNDB ΔG25: 0.000 kcal/mol (A:U −13.64, G:C −16.21); vs Fang measured −13.6 / −15.2 (G:C NNDB 1.0 kcal/mol more stable than measured); GPU path vs reference ≤ 1.2e-5 on all 2-letter 8-mers (3 pairs, 4 temperatures); AS2U − AU = −1.10 per pair |
| 2 extension kinetics | PASS | 30 cases (3 pairs × 2 letters × c = 0.1-10 K_m), r̂ vs Michaelis-Menten, max |z| = 2.03, ~36 000 trials per case |
| 3 pulse integration | PASS | closed form vs explicit (20 sub-steps), max relative difference 0.17% (duplex & bond, τ_p 1 s / 1 min / 1 h) |
| 4 Stage 4 regression | PASS | PAIR=None vs stored chem4c3 results, 16 box2 runs (#1558, #3709 + 14), main + L1 + L2 paired controls: 85 fields bit-identical (also with per-step empty_cache on) |
| 5 activation calibration (new) | PASS | simulator reagent network vs calibration ODE: max 0.10 mM (22 C), 0.03 mM (−3 C, liquid); 2-h A*A at k_act low / fit / high = 0.43 / 2.71 / 4.26 mM brackets observed 2.0-3.3 |
| 6 cold acceptance (new) | PASS (B, C); A residual | see below |
| 7 error threshold | PASS | x = εL_copy 0.15/0.3/0.5/0.7/1/1.5/3: L2_strict 16/16/16/16/0/0/0 of 16, S 0.86/0.65/0.44/0.28/0.18/0.11/0.08 → loss near εL_copy ≈ 0.85 |
| 8 null false positives | PASS | random letters (ε = 0.5), 64 reps, both copy_random controls: L2_strict 0/64 |
| 9 kill and resume | PASS | 640 runs (AS2U, batch 128), kill -9 after 2 batches, resume: IDENTICAL |
First full validation on the original code (52234744): tests 1-5, 7 (threshold 0.5-0.7), 8 (0/64) passed; old test 6 (check) gave
cold factor ≤ 1.15 — the start of the cold-gap work.

### Test 6 (cold acceptance), cold_accept.json
One parameter set must fit B (≥+7 at 5 d, daily re-addition: 7% / 60%) and C (pre-activated ≥+7 at 24 h: 42% RT / 8% +3 C)
within ±5 points; A (≥+1 at 24 h: 65% / 87%) reported with its residual.
| Set | A (0.87) | B (0.60) | C (0.08) |
|---|---|---|---|
| best: K_I 5.04 mM, E_ext 20 | 0.597 (residual ×2.25 per letter, **unexplained**) | 0.593 | 0.114 |
| accepted region | K_I 4.5-6.4 mM with E_ext 19-20 kJ/mol (8 grid points) | | |
| without inhibition (E_ext 18.4) | 0.55 | 0.25 | 0.13 |
**Flag:** the accepted region sits at the edge of the approved ranges: E_ext at the top of 15-20, K_I above the fitted 4.0
(range 2-8). E_ext from the pre-activated control: 18.4 kJ/mol nominal, **10.5-31.4 kJ/mol** with yields ±5 points and room
temperature 20-25 C (eext_uncert.json), so 19-20 is well within its uncertainty.

## Implementation choices
(a) Copy starts are thinned by the first-extension probability (a binding that adds no letter leaves no product); without
this, futile primer bindings (k_on ~1e10 /h) use up 4c-3's one start per template per step.
(b) One activated pool of chemistry 0, letters 1:1: c = 0.5 c_eff M*/V per letter; one block per monomer (the bridged
dinucleotide's second nucleotide is not tracked).
(c) t_half_act is a layer-C dimension (2-4 h, now 2-8 h after the in-situ fit of 5.4 h).
(d) Pulses sit on top of the wet/dry square wave at T_hot = max(T_p, phase T); exact pulse-on time per step; resolved
(step boundaries) when the period ≥ dt_target or one pulse would give > 0.3 breaks on a 64-mer.
(e) Copy starts use the base-temperature time dt − h; duplex melting, bond breaks, activated-block and reagent loss use
the integrated hazard.
(f) act_a: u < 0.1 → 0, else 0.1-500 log; after the cold-gap change it is the reagent supply rate (mM/h).
(g) Validation fixtures: koff_fixed, copy_unlimited (test 2), clamp_primers and PRIMER_LENS = (5,) (tests 7, 8),
reagent_R0, reagent_additions, init_monomer_mM, solute_osm (cold-gap tests). Duplex stacks read 5'→3' along the copy;
mismatches count as pairs; symmetry correction not applied; NN values at 1 M NaCl (no salt/Mg correction).
(h) Inhibitor I = φ(T) × side flux: the calibrated side share y also holds monomeric *A; φ (inhibitory share) is
calibrated on the Fig. 5 impurities (0.120 at 22 C, 0.015 at −3 C; logistic, dE_φ 58 kJ/mol) as a fixed constant. I is
per run (mean field), washed out at 1/τ, not multiplied by c_eff (approved).
(i) Reagent pool per run (mean field over the run's sites), exact within a step for constant supply, monomer at the step
start (operator splitting); side reactions remove the monomer (to "out").
(j) The reagent calibration needs the activation step's own energy E_actv (44 kJ/mol; without it RMS 1.09 vs 0.29 mM):
layer C 0-80 (approved). Nu-pA (catalyst adduct) counted with AMP as unactivated; A*A = model block; *A, AppA, impurities =
side products. ONE-EXPERIMENT CALIBRATION.
(k) Temperature floor −10 C; wet phase floored at T_f = −3.72 C × salinity; liquid only at T ≥ T_f. Aitken's −3 C
liquid relied on ~1 osmolal organics: represented only by the fixture term solute_osm (−1.86 C per osmol/kg).
(l) E_ext layer C 15-20 kJ/mol linear (centre 17.5, not 18.4); K_I 2-8 mM log (centre 4.0).
(m) Memory: optional torch.xpu.empty_cache() every N steps (LOT_EMPTY_CACHE_EVERY; results unchanged, test 4).
(n) The 0 C floor of the 4c-3 decode is lifted only for P3-1 (PAIR set); PAIR=None is 4c-3 unchanged.

## Science flags
- With measured NN thermodynamics, 2-3-mer primers rarely extend (G/C 4-mer k_off ≈ 5e3 /h at 5 C vs r ≈ 30 /h; AS2U/AU
  far worse); copying needs very high primer concentration or longer primers. In 640 random AS2U points (kill/resume runs)
  22 made any copy and none a copy ≥ L_H (diagnostic only, not a map).
- The 24 h single-addition cold benefit (A) is not reproduced: residual ×2.2 per letter, recorded as unexplained.

## Cost and pre-flight
Pre-flight (preflight_p31_AS2U_final.json, final code, LOT_EMPTY_CACHE_EVERY=1, screen fidelity 64 sites x 512 slots x 300
steps, main run + L1 + L2 paired controls; caps: VRAM 27.9 GiB (0.875), host scope MemoryMax 24 GiB; rule <= 75%):
| Batch | worst VRAM reserved (dense set) | worst host peak | host growth over 5 batches | s/run random | s/run dense | pass |
|---|---|---|---|---|---|---|
| 256 | 6.2 GiB (22%) | 4.04 GiB (17%) | 3.0% | 0.35 | 1.16 | yes |
| 512 | 11.5 GiB (41%) | 4.39 GiB (18%) | 4.8% (flat over the last 3) | 0.31 | 0.81 | yes (chosen) |
Without empty_cache (first pre-flight, original code): reserved reached the cap (27.9 GiB at batch 128 dense, 6.8 GiB
allocated) -> fail; maps must run with LOT_EMPTY_CACHE_EVERY=1. Cost per map: 32 768 runs x 0.31 s ~ 2.8 h per pair on
random points (up to ~7.4 h if all points were as slow as the dense set); kill/resume runs (batch 128): 0.46 s/run.

## Temperatures
From gpu_temp.log (all P3-1 jobs, 15:47-18:48): GPU hottest sensor 57-70 C (max 70 C), pkg max 66 C, CPU max 55 C, power up to 136 W; no thermal pause (limits pkg 90 C, VRAM 95 C), never >= 85 C; no OOM; no driver spill (VRAM used per xpu-smi <= ~12 GiB).
