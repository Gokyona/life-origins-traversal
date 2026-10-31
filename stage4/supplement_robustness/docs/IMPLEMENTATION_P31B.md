# P3-1b implementation choices (written before any validation result was read)

30 October 2026, Toronto time. Spec `docs/part3_p31b_spec.md` (e02a9fe, fetched md5 91d0dbc9). Code `chem_p31b.py`, a copy of
`lot-p3/p31/chem_p31.py` (cfc1e3db, not modified).

## Section 1, MULTI_START

| # | Choice | Why / note |
|---|---|---|
| M1 | `MULTI_START` module flag, default True; False runs the chem_p31 code path unchanged (all new code sits behind `if MULTI_START`, the random stream is not touched when off) | spec: off = bit-identical |
| M2 | λ = the old start rate of the template (P3-1: unprimed r(t0)·q(1) + primed k_on[primer]·q(l_p), thinned as in P3-1 choice (a); 4c-3: r + k_on·n_match) | unchanged definition |
| M3 | t_copy = expected copy time from the processivity rule, weighted over the start types by their rates: a deterministic forward pass over positions gives, for each start length, the expected letters added and Σ P(letter j added)/r(t_j); P3-1: the first letter is certain (as in P3-1) | the spec's t_copy = l/r needs l before the copy is drawn; the mean is used |
| M4 | t_melt = 1/k_off of the full-length duplex (length = template length): 4c-3 `k_off(L)` of the step's phase; P3-1 NN k_off integrated over the step's pulse time (k_b (dt−h) + k_hot h)/dt; in cycle-averaged steps (phase index 2) t_melt = max(1/k̄_off, half period), k̄_off the cycle mean the code already uses | spec 1: "expected waiting time to the next phase that melts it, at least half a period" — read as max(half period, 1/mean k_off) |
| M5 | dt_free = copying time of the step × (1 − busy_before/dt), clipped to [0, 1]; a template is eligible under the same condition as before (free, unpaired, busy ≤ dt); P3-1 copying time is the base-temperature time dt − h | spec: "part of the step in which the template is not already busy" |
| M6 | μ = dt_free / (1/λ + t_copy + t_melt); per site the means of the start types are scaled down in proportion where the expected demand exceeds supply: free slots (`room`, after right fragments) for unprimed copies that keep ≥ 2 letters, matching free primers per (site, primer species), then activated blocks (M*₀ after growth) for all copies; n ~ Poisson(Σ scaled means) with the existing `_poisson`; n ≤ ceil(dt_free / (t_copy + t_melt)) | spec 1 limits; "scaled in proportion, as the growth step does" = scaling of the expected demand |
| M7 | each copy drawn independently: start type by the scaled type means, length by the existing processivity loop, errors as before; copies are processed in a random order within each site, so where the realised demand still exceeds supply (Poisson noise) the existing hard limits (primer ranks, cumulative M*₀, slot ranks) cut random copies, i.e. templates in proportion to their n | spec "every template's n is scaled down in proportion" |
| M8 | slot competition with initiation (4c-3 rule, copies vs initiations in proportion to their demands): the copy demand is the expected unprimed slot demand before slot scaling × the released share of unprimed copies in the row (the 4c-3 rule used r·dt per template, uncapped) | keeps 4c-3's slot sharing on the same footing |
| M9 | end of step: all copies but the last of each template (last = last kept copy in the random order) are released as free chains. The last copy started g = U·min(dt_free, cycle) before the end: if g < its copy time it is still being made (template busy for t_copy − g, duplex paired); else it has melted with probability 1 − exp(−(g − t_copy)/t_melt); an unmelted duplex stays paired and melts in the duplex step of later steps | spec "the copy made last ... stays paired ... if the duplex has not melted; the template is busy ... for the remaining copy time" |
| M10 | records (no effect on dynamics), MULTI_START only: `ms_template_steps`, `ms_occupancy_capped_share` (Poisson draw above the cap), `ms_limited_{slots,blocks,primers,any}_share` (template-steps whose mean was scaled down by that resource, or one of whose copies was cut by it); the old `copy_saturated_template_steps` is still reported | spec records |
| M11 | memory: copy entries are expanded with `repeat_interleave`; per-copy random numbers [n, 63] and [n, 64] as in the old per-template code; `rnd_bits` drawn only with copy_random under MULTI_START | |

## Section 1 validation fixtures and criteria (fixed before running)

| Test | Fixture / statistic | Pass |
|---|---|---|
| 1 off switch | (a) PAIR=None, MULTI_START=False, the 16 stored box2 runs of P3-1 test 4 (`regress_4c3_stored.pkl`, which chem_p31 cfc1e3db matched bit for bit), FID (16, 256, 600), seed 4401, main + L1 + L2 paired controls (both); (b) AS2U, 8 points that copy (P3-1 kill-test Sobol set, seed 3101, n 640, indices 287, 16, 324, 552, 31, 184, 328, 303), screen fidelity, seed 3110, main + L1 + L2 (both) controls, chem_p31 vs chem_p31b | 0 differing fields in (a) and (b) |
| 2 low-rate limit | one step per run from free templates (FID (64 sites, 64 slots, 1 step), 32 identical 12-mer templates per site, no growth / initiation / hydrolysis / washout, M* clamped, dt set so λ·dt = 0.02); statistic: copies per template-step, MULTI_START on vs off, same fixture; three regimes: 2a P3-1 GC, k_off fixed 100·r (start-limited: λ(t_copy + t_melt) ≈ 0.01); 2b P3-1 GC, k_off fixed = r (λ(t_copy + t_melt) ≈ 1.5); 2c PAIR=None (4c-3 rule: start rate = extension rate), r = 1e5/h, k_off = 1e6/h all lengths | each regime: |ratio − 1| ≤ 3% and |z| < 3 (spec: "where λ·dt ≤ 0.05 for all templates") |
| 3 rate sensitivity | P3-1 GC, one step, 16 templates all-G and 16 all-C per site (16 sites, 2048 slots), c set so r(C required) / r(G required) = 1/2 (c ≈ 0.165 mM), k_off fixed 100/h, dt 10 h (λ·dt ≥ 10 for both), copy_unlimited (blocks), counts by template from copy_log | ratio of copies within 5% of the renewal formula C_slow/C_fast computed in closed form in the test; with MULTI_START off 1.00 ± 0.02 |
| 4 resource limits | test-3 fixture with (a) blocks low (thermo ×10, M* clamp /10: same c, 1/10 the blocks, no copy_unlimited) and (b) slots low (64 slots per site); per-site totals from the end state | copies never exceed supply in any site (blocks used ≤ M*₀, copies ≤ free slots); ratio of copies between the two templates within 5% of the unlimited ratio |
| 5 kill and resume | `run_p31b.py` (copy of run_p31.py) with MULTI_START on, PAIR=None, the 38 coexistence points tiled, 5 batches of 32, FID (16, 256, 300), kill -9 after 2 batches, resume | IDENTICAL |
| 6 memory | `preflight_p31b.py`, PAIR=None, MULTI_START on, LOT_EMPTY_CACHE_EVERY=1, screen fidelity, batch 512, sets coex, box2, coex ×4 (coex = the 38 4c-5 point-study points, 4c-3 dimensions, tiled), main + L1 + L2 controls, scope MemoryMax 24 G | peak XPU reserved ≤ 80% of the 27.9 GiB cap; host peak ≤ 80% of 24 GiB; no growth (≤ 10%) over the later batches |

## Correction after the first validation run (chem_p31b md5 21085a25 -> new md5 below)
| # | Change | Found by |
|---|---|---|
| M8b | slot sharing with initiation under MULTI_START: copy demand = max(expected demand (M8), unprimed copies actually made). With M8 alone, the slot rule capped the copies at the expected demand rounded stochastically even when free slots were plentiful (no initiation), cutting realised copies above the mean; the spec cuts only where demand exceeds supply | test 2 (2a ratio 0.976, z −1.9; 2b, 2c far below the renewal expectation) and the test-3 record `limited_any_share` 0.19 with ample slots |

## Review 1 (reviewer: coordinator session, 30 October 2026 21:35 Toronto time; recorded before anything else was run)
Context: the run was resumed at 21:15, and all earlier outputs were moved to `void_2026-10-30/`. On the fresh rerun with chem_p31b 53711ca2, tests 1, 3 and 4 PASS and test 2 FAILS (2a, 2b, 2c). Tests 5 and 6 were not run. Test 6 had failed earlier on 21085a25 (peak reserved 27.70 GiB). M8b was applied after a failure without review.

| # | Decision |
|---|---|
| R1.1 | Test 2: 2b is a fixture error. λ(t_copy+t_melt) = 1.5 is not a low-rate case; occupancy is tested by test 3. 2b is dropped from test 2, and its numbers stay in the record. For 2a and 2c no code change yet: first a like-for-like CPU diagnosis that decomposes the on and off counts (time base, primer/okp filter, q_at thinning, busy/paired filter, slot/block limit). Any fix goes to the reviewer as a diff before it is applied. |
| R1.2 | M8b: the diff and its rationale go to the reviewer. Until it is approved, 53711ca2 is provisional and no results on it count. |
| R1.3 | Test 6: a smaller batch is allowed (operational setting). After test 2 is resolved, the pre-flight runs at batch 256, and at 128 if 256 fails. Criteria unchanged: peak ≤ 80% of cap, no growth over 5 batches, LOT_EMPTY_CACHE_EVERY=1. |

## Review 2 (reviewer: coordinator session, 30 October 2026 ~22:10 Toronto time; recorded before anything was run)
Diagnosis (review1_diag/, CPU, counters only): start counts, time base, q_at thinning, the primer, ≥2-letter and block filters are identical on and off. The whole test-2 gap is the slot-sharing step. The inherited rule keep_c = min(nc, floor(min(room, w_c + tot_init)·frac_c + U)) caps a site's realised unprimed copies at its stochastically rounded expected demand, even with ample room and no initiation (≈ half the copies at w_c ≈ 0.64 per site). M8b removed that cap under MULTI_START.

| # | Decision |
|---|---|
| R2.a | M8b accepted as the rule "copies are cut only where demand exceeds free slots" (spec §1). It is always on under MULTI_START. |
| R2.b | Test 2 statistic = copy tries after the start, primer, ≥2-letter and block filters, before the slot step (`copy_tries`). Post-slot copies are reported, not gated. Criterion unchanged (within 3%, |z| < 3). 2b dropped. |
| R2.c | 2c made low-rate: T raised above T_m8 so that λ(t_copy + t_melt) ≈ 0.02. Same criterion. |
| R2.d | New switch `SLOT_ROUND_FIX` (default False, so OFF is bit-identical to cfc1e3db) applies the M8b rule in the MULTI_START-off path. New test 2d: fixture 2a, SLOT_ROUND_FIX on, post-slot copies ON vs OFF within 3%, |z| < 3. §3 gets variant A2a (slot fix only) next to A2 (multi-start + slot fix). Records (no effect on dynamics, no random numbers used) in every run: `round_cut_site_steps` (site-steps where the rounding cut removed copies that the fixed rule would keep) and `round_cut_copies_lost` (the number of those copies). |
| R2.e | The rounding cut goes on the assumptions list once its size at the published points is known from the base runs. |
| R2.f | Then tests 1–4 (with 2d), test 5, and test 6 (batch 256, 128 if needed). If all pass, continue into §2 and §3 without waiting. Any failure: stop and report. |

### Implementation of Review 2 (before any result on it was read): chem_p31b 53711ca2 -> 192b9047, validate_p31b b438f5f2 -> 7055684d
| # | Choice |
|---|---|
| R2.d-impl | `SLOT_ROUND_FIX` sets w_c = max(w_c, nc) in the P3-1 and 4c-3 branches. Under MULTI_START it is always applied (M8b). Records: keep_f is computed with the fixed rule and the same u_keep, and `round_cut_copies_lost` = Σ max(keep_f − keep_c, 0), `round_cut_site_steps` = the number of site-steps with a loss. No random numbers are used, so OFF stays bit-identical. |
| R2.c-impl | Raising T alone cannot reach λ(t_copy+t_melt) ≈ 0.02 under the 4c-3 rule: k_off(1) = K0 = 1e6/h at any T, so q(1) = r/(r+K0) and λ·t_copy ≥ q(1) = 0.09 at r = 1e5/h. 2c therefore uses T = T_m8 + 10 °C **and** r = k_t = 1e4/h, which gives λ(t_copy+t_melt) ≈ 0.01, computed in the test from the model's k_off at dH_bp = 35 kJ/mol. Copy tries are ≈ 2e-4 per template-step, so 32 batches × 2048 runs × 64 sites are pooled. |
| R2.d-test | 2d: fixture 2a, seeds 7011/7012, SLOT_ROUND_FIX on, statistic = post-slot copies, criterion within 3%, |z| < 3. |

## Review 3 (reviewer: coordinator session, 30 October 2026 ~22:30 Toronto time; recorded before anything was run)
Finding (scratch diagnosis): torch.poisson on the XPU. 8 seeds × 4.2M draws: at λ = 0.02 the per-seed means are 0.920–1.050 × λ and the pooled mean is −1.1% (z −8.6); at λ = 0.5 it is −0.2% (z −9). On CPU the check is clean. This is the cause of the 2d failure (ratio 1.056, z 5.4). Test 1 failed only on the two new record keys; all 85 + 96 original fields are bit-identical.

| # | Decision |
|---|---|
| R3.1 | New switch `POISSON_INV` (default False, OFF bit-identical): exact CDF inversion from torch.rand on the same generator up to POISSON_EXACT_MAX, the existing normal approximation above. It is on whenever MULTI_START is on, and also available on its own. New test 2e (XPU), criteria: (i) torch.rand, 8 seeds × 4M draws: mean and variance pooled |z| < 3, between-seed scatter consistent with iid (chi-square p > 0.01); (ii) the inversion sampler at λ = 0.02, 0.5, 5, 50: the same two criteria plus pooled mean within 0.5%; (iii) record only: torch.poisson on XPU and CPU. 2d is rerun on 6 seeds pooled, same criterion. |
| R3.2 | Test 1 compares only the reference's keys; the new record fields are reported, not compared. |
| R3.3 | 2c fixture (r = k_t = 1e4/h, T = T_m8 + 10 °C) accepted. |
| R3.4 | §3 gets variant A0 = base + POISSON_INV only, at all points (sizes the torch.poisson effect on the published Stage 4 levels). |
| R3.5 | Then tests 1–4 (with 2d, 2e), test 5, test 6 (batch 256, 128 if needed). If all pass, continue into §2 and §3; report after test 6 with s/run and a firm ETA. Any failure: stop and report. |

### Implementation of Review 3 (before any result on it was read): chem_p31b 192b9047 -> a31cae30, validate_p31b 7055684d -> 6b0b844d
| # | Choice |
|---|---|
| R3.1-impl | `_poisson_inv`: float64 on the device. One torch.rand(float64) uniform per element from the same generator. The search starts at the mode m = floor(λ), with F = Q(m+1, λ) (torch.special.gammaincc) and pmf p = exp(−λ + m ln λ − lgamma(m+1)), and walks up (p ← p·λ/k) or down (p ← p·x/λ) to the smallest x with u ≤ F(x). At most 10√λ_max + 100 steps per direction. Used for every `_poisson` call when POISSON_INV or MULTI_START; the normal approximation above POISSON_EXACT_MAX is unchanged. |
| R3.1-test | 2e = `test_e` in validate_p31b.py. Seeds 7000–7007 × 2²² draws. Variance of the sample variance from μ4 (uniform 1/80; Poisson λ(1+3λ)). χ² with 8 degrees of freedom on the per-seed means and on the per-seed variances. 2d: seeds 7011+10k (on) / 7012+10k (off), k = 0..5, pooled. |

### Tests 5 and 6 harness (written before running)
| # | Choice |
|---|---|
| T5 | `run_p31b.py` (aec25b4b) is a new runner modelled on run_p31.py: PAIR=None, MULTI_START on, the 38 coex points tiled, 5 batches × 32, FID (16, 256, 300). Seed = f(tag, batch). Every per-run output of simulate ([B] and [B, ≤64]) plus both paired controls (L1; L2 with both=True) goes into one fsynced JSON line per batch. `kill_resume_p31b.sh` (3283ff83): uninterrupted run vs kill -9 after 2 batches + resume; all fields except seconds must be identical. |
| T6 | `preflight_p31b.py` (ccd93c30): the L2 paired control now runs with both=True, as in the runner (conservative). Sets coex, box2, coex ×4 (first + 5 later batches), screen fidelity, MULTI_START on (so POISSON_INV on), LOT_EMPTY_CACHE_EVERY=1, scope MemoryMax 24 G. Batch 256, then 128 if 256 fails. Queue: `queue_s1.sh`. |

## Review 4 (reviewer: coordinator session, 30 October 2026 23:25 Toronto time)
Status at review: tests 1–4 and 2e PASS (validate_p31b.json d22f26dc, chem_p31b a31cae30). Test 5 PASS (5 + 5 batches, 160 runs, IDENTICAL; killed after 2 batches). Test 6 running. §2 sources: S2_sources_p31b.md.

| # | Decision |
|---|---|
| R4.1 | Salt (§2.2), option 2. `SALT_CORR` = the Chen & Znosko (2013) RNA Na⁺ length form ΔG°37(2) = ΔG°37(1) − 0.1016·N·ln([Na⁺]/1.021 M), entering through ΔS, validity 71–1021 mM. Off for any run with Mg²⁺. The Mg/Na mixing gap is recorded (Owczarzy 2008 unreadable). Arteaga et al. (2023) is recorded but not applied. |
| R4.2 | §3 on the 7 L2 points (#1558, #3709, engine435/5704/10653 = A/B/C, rep1_1926, rep1_9823): variants base, A0, A1, A2, A2a, A4, A5 (shortest-period points first, capped as in the spec) and A6. 16 seeds, paired controls, the published confirm fidelity (16, 256, 600). Start when test 6 passes. Largest batch that passed test 6; points and variants packed into big batches. |
| R4.3 | Coexistence groups and A3 by a port: chem4c5 (lot-s4, never modified) copied to p31b/chem4c5b.py with MULTI_START, POISSON_INV and SLOT_ROUND_FIX as switches, all default off. Validation before any coex run: (a) all switches off, bit-identical to chem4c5 on 3 seeds at #3709, engine12538 and one point without amphiphiles (the 4c-5 test-1 set), all reference keys; (b) tests 2a, 2d and 2e on chem4c5b; (c) pre-flight at the 4c-5 worst case (record on), ≤ 80% of cap, no growth over 5 batches. Then the 10 coex groups with base, A0, A2, A2a, A3 and A4, 16 seeds, 4c-5 paired controls (β = 0 and random-copy), starting right after R4.2 with no GPU gap. Any failure: stop and report. |
| R4.4 | Detector perturbation and the round-cut / Poisson records as in the spec, on all §3 runs. |

## §3 implementation choices (written before any §3 run; chem_p31b a31cae30 -> 60c68e27, run_s3_l2.py c41e673b)
| # | Choice |
|---|---|
| S1 | Ablation switches in chem_p31b, all off by default (the off path is unchanged code; test 1 is rerun on 60c68e27 before §3): `DUPLEX_HYDRO_FACTOR` (A4: 0.1 × the hydrolysis rate of paired strands, the UV term unchanged); `LONG_WASH_FACTOR` / `LONG_WASH_MIN_LEN` (A6: chains with L ≥ 4 wash out at 0.5/τ); `DT_PERIOD_FRAC` (A5: at cycling points no cycle averaging and dt ≤ period/4; the step cap is raised to 1.2 × t_end/(period/4) + 2 n_steps). A1 = module constant G_CAP 4 → 16. A0 = POISSON_INV. A2 = MULTI_START (which brings the slot fix and POISSON_INV). A2a = SLOT_ROUND_FIX. |
| S2 | Points: #1558 and #3709 (validate_p31._box2_points, the test-1 coordinates); A, B, C = engine435, engine5704, engine10653 (confirm4c3.json u_model); rep1_1926, rep1_9823 (confirm4c3_rep1.json u_model). FID (16, 256, 600), 16 seeds per point, as in the published confirm runs. |
| S3 | One batch per variant: all 7 points × 16 seeds = 112 runs (variants are module switches, so they cannot share a batch). Seed 83100 for every variant (common random numbers across variants; the base is not seed-identical to the published runs, so it is compared statistically). L2 controls: paired_control_l2(both=True), i.e. two random-copy controls. All per-run outputs and the controls are stored (one JSON line per batch). |
| S4 | Base vs published: per point, published L2_strict k/16 vs base k/16, Fisher exact two-sided p, and the spec reading (robust if |Δ| < half the base value or < 3 of 16). |
| S5 | A5 cost estimate before it runs: the steps needed per run = t_end/(period/4) from the base runs' t_end, × the measured s per step. A5 runs at the 3 shortest-period L2 points only if the total is < 4 h (spec). Otherwise it is not run and that is reported. |
| S6 | Detector perturbation. Margin 0.10/0.15/0.20 is recomputed from the stored S_main and S_ctl (strict = candidate ∧ S_main − S_ctl ≥ margin). L_H 5/7 and persistence 4/6 turnovers need per-checkpoint shares the code does not store, so they are obtained by deterministic replay of the base batch with the same seed and the detector constant changed (the dynamics do not read L_H or L2_SPAN). The replay's dynamics outputs are checked to be identical to the base ("det_*" variants). The controls' false-positive rate at each threshold comes from the same runs. |

## R4.3 port: implementation choices (before any chem4c5b run; chem4c4b fa9e35b3, chem4c5b ae3d7d91, validate_4c5b b839dd27, run_s3_coex e957e926)
| # | Choice |
|---|---|
| C1 | chem4c4b = lot-s4/chem4c4.py + the chem_p31b 4c-3-branch code for MULTI_START, SLOT_ROUND_FIX, POISSON_INV (incl. records), and A4 DUPLEX_HYDRO_FACTOR. All default off. Compartments: primer supply and demand are per (site, compartment) (trc), the template's extension rate already includes F_PERM inside vesicles, a copy inherits its template's compartment (tcm reindexed per copy). chem4c5b = lot-s4/chem4c5.py importing chem4c4b; switches are set on chem4c4b. |
| C2 | Validation (a): lot-s4 chem4c5 vs chem4c5b, record on, switches off, seeds 101–103, 4 runs each at #3709, engine12538, #3709 amph 0, FID (16, 256, 600). All tensor outputs of the reference, the final state (L, seq, chem) and repr(L3_details) must be equal. The lot-s4 modules are loaded by path (read only). |
| C3 | Validation (b): the P3-1 GC fixture of 2a/2d does not exist in the 4c-4 chemistry, so 2a and 2d use the low-rate 4c-3-rule fixture of P3-1b test 2c (r = 1e4/h, T = T_m8 + 10 °C, amph 0), 32 batches × 2048 runs × 64 sites. 2a: copy tries ON vs OFF. 2d: post-slot copies ON vs OFF with SLOT_ROUND_FIX on. Both within 3%, |z| < 3. 2e = validate_p31b.test_e on chem4c4b._poisson_inv. |
| C4 | Validation (c): the validate4c5 test-8 worst case (9a sets + engine12538, window open, phi0 0.3, s_shear 10, f_rehyd 10, record on, L1 and L2 controls, paired_inherit on every evaluable run), MULTI_START on, screen fidelity, batch 256 (128 if 256 fails), 6 batches. Peak XPU reserved ≤ 80% of the cap, host ≤ 80% of 24 GiB, growth ≤ 10% after the first batch. |
| C5 | Coex runs: run_s3_coex.py = run4c5_5a.py with chem4c5b; variants base, A0, A2, A2a, A3 (override f_rehyd = 1, passed to the main run and both controls), A4. Same points, seeds (95000 + batch), screen fidelity, 16 seeds, batch 256, controls (β = 0, random-copy), gate and records as published, plus copies, round-cut and MULTI_START records per run. One directory per variant. |

## Review 5 (reviewer: coordinator session, 31 October 2026 01:30 Toronto time)
| # | Decision |
|---|---|
| R5.1 | Speed-up of _poisson_inv pre-approved on these conditions: identical arithmetic per element, only the active elements compacted per iteration. Before use, exact equality with the current function on ≥ 1e8 draws on the XPU, same seeds, λ ∈ {0.02, 0.5, 5, 50, mixed arrays}. Then test 2e(ii) and test 1 (reference keys) on the new code. If all are identical and pass, swap it in; otherwise stop and report. Finished A0 runs stay valid (identical draws). |
| R5.2 | Base vs published as binomial tests per L2 point; flag p < 0.05. |
| R5.3 | Light GPU load during the small L2 batches is acceptable; never two GPU jobs at once. |

### R5.1 implementation (before the equality check)
The cost is the number of loop iterations (each a host sync; up to ~5√λ_max ≈ 500 at λ_max = 1e4), not the element count. `_poisson_inv` (new) runs the same per-element update in blocks of 32 iterations between host checks. Finished elements are left untouched by the masked update, exactly as before. Before each block the active elements are compacted into index lists. The total iteration cap (nmax) is unchanged. The per-element operation sequence is unchanged, so the draws must be bit-identical; that is checked as R5.1 requires.

## Review 6 (reviewer: coordinator session, 31 October 2026 01:40 Toronto time), supersedes the R5.1 implementation
The Cornish–Fisher sampler is approved as proposed. For λ > 10: start at x0 = round(λ + √λ·z + (z² − 1)/6), z = Φ⁻¹(u), then exact inversion (the smallest x with u ≤ F(x), F(x) = Q(x + 1, λ) by gammaincc, float64). For λ ≤ 10 the walk is unchanged. Same uniform u per element (the first torch.rand draw, as before).

Gates, fixed before running:
- 2e(ii) at λ = 0.02, 0.5, 5, 50, 1e3 and 1e4, same criteria;
- ≥ 99.9999% agreement with the current sampler on 1e8 draws;
- the search always terminates: maximum iterations logged, hard cap (200) raising an error, never a silent fallback;
- test 1 (reference keys) passes.

If all pass:
- swap the sampler into chem_p31b and chem4c4b before chem4c5b validation (c);
- rerun L2 A0 and A2 (variants A0_r6, A2_r6; the old A0/A2 kept as records);
- use the new sampler for all coex runs.

If any gate fails: keep the current sampler, continue as queued, and report. The R5.1 blocked/compacted variant is not used.

## Review 7 (reviewer: coordinator session, 31 October 2026 02:25 Toronto time)
A2 (MULTI_START) loses L2 at all 7 points; this must be decomposed before anything is concluded. Approved diagnostics, each one 112-run batch at the 7 L2 points, after the L2 r6 reruns and before the chem4c5b validation:
- A2b: base with the slot cut in random order instead of slot-index order, no other change;
- A2c: MULTI_START with t_melt = 1/k̄_off (the same cycle-averaged k_off as the old rule), no half-period floor;
- A2d: MULTI_START with the half-period floor, but slots cut in slot-index order as in the old rule.

Also report:
- from the A2 records: per-step share of copies lost to slots, blocks and primers; melt-wait share; copies per template per step; mean copy length;
- A4: paired-strand break events in base;
- A5: the cost estimate.

Reading, fixed now: if A2c recovers L2, the collapse comes from the melt-wait rule. If A2b loses L2, the published L2 depended on slot-order priority (a numerical artefact). The spec rules are not changed by me. Coex may start afterwards as planned.

## Review 8 (reviewer: coordinator session, 31 October 2026 ~02:30 Toronto time)
A5 approved at the 5 feasible points. engine435 and engine5704 are recorded as "not testable at affordable cost (~105 h each)". A5 is the phase-resolved reference that both melt rules approximate. Two parts: A5_part1 (p3709, rep1_9823, engine10653) right after the Review 7 diagnostics; A5_part2 (rep1_1926, p1558) at the very end of the queue.

### Implementation of Reviews 7 and 8 (before running; chem_p31b 990c391a -> new md5 below; test 1 is rerun before the diagnostics)
| # | Choice |
|---|---|
| D1 | `SLOT_CUT_RANDOM` (A2b, old path only): the slot cut ranks a row's unprimed copies in a random order (random keys from the auxiliary generator `aux`, so the main random stream is unchanged) instead of slot-index (row-major) order. Everything else is as in base. |
| D2 | `MS_MELT_FLOOR` (default True; A2c False): under MULTI_START, drop the half-period floor on t_melt in cycle-averaged steps. t_melt = 1/k_off(L) of the step's phase index, which for phase 2 is the cycle mean of the two k_off values, as the old rule uses. |
| D3 | `MS_ORDER_SLOT` (A2d): under MULTI_START, the copies of a site are ordered by template slot index (stable) instead of at random, so primer, block and slot cuts fall in slot-index order as in the old rule. The random keys are still drawn (stream unchanged). |
| D4 | Records (no random numbers): `break_events` and `paired_break_events` per run (bond breaks, and those in strands paired at the time), counted in every run. |
| D5 | Variants in run_s3_l2: base_rec (base with the records, for A4), A2b, A2c, A2d, A5_part1 / A5_part2 (DT_PERIOD_FRAC 0.25, points as in Review 8). |
| D6 | Records added under MULTI_START for Review 7 (no random numbers; chem_p31b b93e363f): ms_melt_wait_share (mean t_melt/cycle over template-steps); ms_copies_drawn_per_template_step; ms_expected_per_template_step (before resource scaling); ms_lossS_{slots,primers,blocks} (expected copies removed by the proportional scaling); ms_lossH_{slots,primers,blocks} (copies cut at the hard limits). Variant A2_rec = A2 rerun with these records. Test 1 (queued) covers the off path of this version. |

## Review 9 (reviewer: coordinator session, 31 October 2026 ~03:10 Toronto time; recorded before any result on it)
| # | Decision |
|---|---|
| R9.1 | A2e = MULTI_START on, MS_MELT_FLOOR False, MS_ORDER_SLOT True; 7 L2 points × 16 seeds, the same records as A2/A2c. Reading, fixed now: if A2e recovers L2 to roughly the A2b/base level, L2 needs priority ordering plus copy counts, not rate differences; if it stays near A2c, slot-index priority alone does not explain it under multi-start. |
| R9.2 | Coex variants: base, A0, A2, A2a, A2b, A2c, A3, A4. Pre-flight on the batch size for A2c (it restores copy counts under multi-start). Store the L3 records needed for the detector perturbation for A2b (and base). |
| R9.3 | If cheap and without effect on the dynamics: a per-run record of the concentration of copies by template lineage (Gini, top-10% share) for base, A2b, A2c. It must pass a bit-identity check on the stored fields first; otherwise skip and say so. |

### Implementation of Review 9 (before running)
| # | Choice |
|---|---|
| R9-i1 | Copy concentration record (`COPY_CONC_BINS` = 4096, chem_p31b and chem4c4b, record only, no random numbers). Each released copy adds 1 to bin (lineage key mod 4096) of its run. Lineage key = min(species key of the template, species key of its complement), i.e. the template sequence with both strands as one lineage. At the end, per run: `copy_conc_top10` = share of copies in the top 10% of occupied bins, `copy_conc_gini` = Gini over occupied bins, `copy_conc_bins` = occupied bins. Bin collisions merge lineages (4096 bins against ≤ ~1e4 distinct template lineages per run), which biases both measures down slightly and equally across variants. Bit-identity: test 1, plus L2 reruns base_conc / A2b_conc / A2c_conc that must equal base / A2b / A2c in every previously stored field. |
| R9-i2 | chem4c4b gets SLOT_CUT_RANDOM, MS_MELT_FLOOR and MS_ORDER_SLOT exactly as chem_p31b (D1–D3), and the copy concentration record. chem4c5b validation (a) checks the off path. |
| R9-i3 | Pre-flight (validate_4c5b test c) runs under the A2c settings (MULTI_START, MS_MELT_FLOOR False), the heaviest multi-start variant (copy counts restored), at batch 256 (128 if 256 fails). The coex batch follows it. |
| R9-i4 | L3 detector perturbation (base and A2b coex runs): the chem4c5b.paired_inherit hook `PERTURB_SETS` re-evaluates local_inherit on the lineage records already in memory (main, main', random-copy; and for each control its own runs) with L3_IDENT 0.85/0.95, L3_FRAC 0.4/0.6, L3_UP 2/4 (one changed at a time). No new simulation, no random numbers. Perturbed strict = gate ∧ L3_loose ∧ main share ≥ 0.5 ∧ share(main') − share(c) ≥ 0.25. Limitation: only runs in the default paired mask (gate ∧ default L3_strict) get the paired shares, so a looser setting cannot add runs outside that mask. This is reported. L3 loose ρ thresholds 0.05/0.1/0.2 are recomputed from the stored ρ and p of main and controls. |
| R9-i5 | Code after Review 9: chem_p31b 44f30af0, chem4c4b 06aeb226, chem4c5b 6abf1b2b, run_s3_l2 f64b35eb, run_s3_coex a64447c3, validate_4c5b d7e1c71c, check_conc_identity.py. Queue queue_s3f.sh: A5 part 1 (running) → test 1 → L2 A2e, base_conc, A2b_conc, A2c_conc → concentration identity check → chem4c5b validation a/b/c (c under A2c) → coex base, A2b, A2c, A0, A2, A2a, A3, A4 → L2 detector replays → A5 part 2. |

## Review 10 (reviewer: coordinator session, 31 October 2026 ~05:25 Toronto time)
Pre-flight (c) under A2c failed at batch 256: peak reserved 25.26 / 26.94 GiB (allocated 11.06 / 12.28 GiB) against the 22.3 GiB limit. It was stopped after batch 1. Cause: chem4c4 has no EMPTY_CACHE_EVERY, so LOT_EMPTY_CACHE_EVERY=1 had no effect. Approved fix (option A): EMPTY_CACHE_EVERY ported into chem4c4b (chem4c4b_ec.py 61d9a4f4 → chem4c4b.py). It calls torch.xpu.empty_cache() every n steps (env LOT_EMPTY_CACHE_EVERY, default 0). Memory only; it is the fix reviewed for P3-1 in chem_p31b.

Gates:
- chem4c5b validation (a), bit-identity vs lot-s4 chem4c5, with LOT_EMPTY_CACHE_EVERY=1. Test 1 is a chem_p31b test, so (a) is the check of the chem4c4b path.
- Then (c) at batch 256, same criteria (reserved ≤ 22.3 GiB, no growth over 5 batches), 128 if 256 fails.
- The coex batch is whichever passes.
- The queue then resumes unchanged.

The failed (c) record is kept in logs/validate_4c5b_c_failed_06aeb226.json.
