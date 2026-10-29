# Stage 4c-4 validation: amphiphiles and vesicles (Intel Arc Pro B70)

**Status:** complete, 28 October 2026, 23:19 Toronto time (validation passed 02:08); report: [stage4c4_vesicles_b70.md](stage4c4_vesicles_b70.md) · **Specification:** [stage4c4_spec.md](../docs/stage4c4_spec.md), with the implementation choices in its §8 · **Code:** chem4c4.py (v2, md5 b79360f6), validate4c4.py (765d5e41), check3709_4c4.py (b2a60bdf); kept on the benchmark machine and to be added to [stage4/](../stage4/README.md) · **Data:** result JSONs with md5 sums on the benchmark machine, to be added to the repository

## Bottom line so far

- The vesicle layer passes its validation tests. Without amphiphiles it is bit-identical to 4c-3.
- At the heredity point #3709, heredity (L2) and persistent vesicles that hold chains (E3 strict) occur together in 16 of 16 runs, and all lineage mass is inside vesicles. This needs membrane constants that allow bilayers at pH 4.6 and 0.77 M salt (see below). Whether such constants are realistic, and how common this overlap is, is what the map has to show.

## Memory pre-flight (test 9)

The worst case is a batch made only of the most template-rich heredity points (#3709, A, C), with the membrane window open and all vesicle slots in use.

| Batch | Result |
|---|---|
| 1024 | failed (GPU out of memory after about 4 s) |
| 512 | failed |
| 256 | **pass**: allocated 14.0 GiB (50% of the 27.9 GiB cap), host peak 4.46 of 24 GiB |

- The failures do not come from the vesicle layer. At these points nearly every chain slot holds a template, and the per-template arrays of the copy step (from 4c-3) exceed the cap at batch 512 and above. chem4c3 has the same exposure; random batches were never this dense.
- The map will use batch 256. Cost at the densest points: 2.64 s per run, about 4 times chem4c3 on random focused-region points.

## Design changes found during the pre-flight (before any map)

The pre-flight diagnostics showed two rules that made E3 impossible at exactly the heredity points. Both were changed before any validation result (spec §8, rows 8, 10, 14):

1. **No enclosure.** Vesicles formed in the first steps, before chains existed, filled all 8 slots and persisted (slow washout). Rehydration, the other route for enclosure, only acted at resolved dry-to-wet switches, and the heredity points are cycle-averaged. Enclosed chain share was 0.000. Now rehydration in cycle-averaged runs acts as an average over the cycles in each step.
2. **E3 not evaluable.** E3 used only the last 8 checkpoints, which needs a horizon of 10 τ. Now the E3 window is built like the L2 window: the longest run of checkpoints ending at the last one, at least 5 vesicle turnovers.

A run extension to 10 vesicle turnovers was also dropped, so the vesicle layer changes neither the horizon nor the time step.

## Validation tests

| Test | v1 | v2 |
|---|---|---|
| 1. No amphiphiles: bit-identical to chem4c3 | pass | pass |
| 2. CVC: membrane = total − CVC above the CVC, 0 below | pass | not affected |
| 3. Stability window: vesicles only inside; breakdown at k_break when a bound is crossed | pass | pass |
| 4. Encapsulation at formation is binomial with p = φ_j (χ²) | pass | pass |
| 5. Division: binomial ½ split, duplexes together, membrane conserved | pass | not affected |
| 6. Osmotic competition: loaded vesicles gain membrane with β > 0, no difference with β = 0 | pass | not affected |
| 7. Growth–division steady state; no division without shear or cycling | pass | pass |
| 8. Kill and resume bit-identical | – | pass (killed after 3 engine batches; runs, selections, compositions, curves and final maps identical) |
| 9. Memory pre-flight | pass at batch 256 | |

Criteria were written into the output file before each test ran.

## Heredity point #3709 with vesicles (v2)

16 runs, same seeds and identical horizon (44,720 h) in both arms. #3709 lies outside the membrane window at the default (mid-range) membrane constants: its pH is 4.57 and its salinity 0.77 M. The check therefore uses membrane constants inside their layer C ranges that allow bilayers there: amphiphile fraction 0.3, CVC₀ 0.1 mM, pKa_app 6, w_pH 2, sal_crit 3 M, φ₀ 0.1, f_rehyd 3.

| | chem4c3 | chem4c4 v2 |
|---|---|---|
| L2 strict / 16 | 16 | 16 |
| L2_long / 16 | 16 | 16 |
| Margin | 0.613 | 0.341 |
| Copy length | 63.2 | 62.6 |
| E3 evaluable / loose / strict | – | 16 / 16 / 16 |
| Enclosed share of chain mass | – | 1.000 |
| Enclosed share of L2 lineage mass | – | 1.000 |
| Vesicles per site | – | 7.84 |

- Heredity survives enclosure. The margin over the random-copy controls drops from 0.61 to 0.34, probably because each site's lineages are split over about 8 compartments.
- This shows that L2 and E3 can coexist in the model. It does not show that they do so for realistic membranes: fatty-acid bilayers usually need a pH near their apparent pKa and low salt, and #3709 is acidic and salty. The map over all 54 dimensions will show how much of the heredity region overlaps the membrane window.

Result files (md5): validate4c4.json (v1) fe0bf451, validate4c4_v2.json a6a56fdf, check3709_4c4.json f39224ed, with the failed and stopped pre-flight records.

## Full map (54 dimensions)

28 October 2026, 02:09 to 12:33 Toronto time (20,724 s). Ground truth 1000 × 16 runs, COMPASS 32,000 runs, batch 256. Result file md5 30bdb8c5.

| Level | Ground-truth fraction (P > 0.5) | COMPASS runs positive | COMPASS accuracy | 4c-3 full map |
|---|---|---|---|---|
| E1 (activated feedstock) | 0.450 | 18,572 | 0.975 | 0.466 |
| E2 (long chains) | 0.107 | 4,429 | 0.957 | 0.085 |
| L1 (catalytic network) | 0 | 21 | – | 0.001 |
| L2 strict (heredity) | 0 | 55 | – | 0 |
| L2_long | 0 | 1 | – | 0 |
| E3 loose / strict (persistent vesicles holding chains) | 0.002 / 0.002 | 43 / 41 | – | – |
| E3 strict ∩ L2_long | 0 | 0 | – | – |

- **Vesicles are rare.** Ground truth has 2 conditions of 1000 with E3, both at 16 of 16 seeds; all their chains are enclosed. Both also have E1 and E2, and neither has L1 or L2.
- **E3 is never found together with heredity** in the full map: 0 in ground truth and 0 in 32,000 COMPASS runs. In the E3 runs, copying is active but copies are short (median about 2.7 letters).
- E3 is evaluable (≥ 5 vesicle turnovers) in 69% of ground-truth runs.
- Importance for E1: τ 0.26, T 0.13, Y_UV 0.10, feedstock 0.10, thermal cycling 0.09, UV 0.09, H₂ 0.07. For E2: τ, T, feedstock, Y_UV, UV, H₂, pH. E3, L1 and L2 are too rare for a classifier; their accuracy only reflects the all-negative baseline.
- E2 is higher than in 4c-3 (0.107 against 0.085). **This is not caused by vesicles.** A paired rerun of all 16,000 ground-truth runs with the same points and seeds but no amphiphiles (e2attr4c4.json, md5 30104396) gives E2 0.1127 without amphiphiles against 0.1056 with them. The 50 runs that reach E2 only with amphiphiles have no vesicles at all. Amphiphiles slightly lower E2, because they take a share of the feedstock. The difference from 4c-3 comes from the point set: the ground-truth points are drawn from the same random stream, but with 54 instead of 42 dimensions, so they are different points.
- The two ground-truth E3 conditions are both mild (about 21 °C) with pH inside the window of a matching pKa_app (6.5 and 8.6), Mg²⁺ well below Mg_crit, and an amphiphile fraction just above 0.1.

## Membrane scan at the heredity points (spec §9a)

28 October 2026, 13:54 to 14:55 Toronto time. Result file scan4c4mem.json, md5 46eed70d. For each heredity point, environment and copying constants are fixed. 64 sets of membrane constants are drawn within their layer C ranges, restricted to sets under which the point lies inside the membrane stability window, with an amphiphile fraction ≥ 0.1. Each set runs with 4 seeds and both random-copy controls. The reference is chem4c3 at the same point, 16 runs with the same seed.

| Point | Share of membrane space compatible | E3 strict | L2_long | Both | Sets with both (≥ 3 of 4 seeds) | Enclosed share of L2 lineage | Margin | chem4c3 reference: L2_long / 16, margin |
|---|---|---|---|---|---|---|---|---|
| #3709 (38.7 °C, pH 4.57) | 0.004 | 1.00 | 0.28 | **0.28** | 0.23 | 0.998 | 0.25 | 16, 0.51 |
| #1558 (7.6 °C, pH 7.96) | 0.18 | 1.00 | 0.64 | **0.64** | 0.59 | 1.000 | 0.28 | 9, 0.25 |
| C (26.7 °C, pH 4.73) | 0.014 | 1.00 | 0.86 | **0.86** | 0.77 | 1.000 | 0.34 | 14, 0.31 |
| rep1_1926 (7.8 °C, pH 6.73) | 0.11 | 0.98 | 0.99 | **0.98** | 0.98 | 0.984 | 0.51 | 16, 0.52 |

- **Heredity and persistent vesicles occur together at all four points** where a membrane can be stable. Wherever the membrane window allows it, E3 strict follows almost always, and heredity survives in 28–98% of runs, with essentially all lineage mass inside vesicles and copies of 52–60 letters.
- **Enclosure costs little at three points.** At C and rep1_1926 the joint rate is close to the rate without membranes. At #1558 the reference reached only 9 of 16 at this seed (16 of 16 at its original seed), so the 64% with membranes is no lower. At #3709, the warmest and most acidic point, enclosure costs a lot (28% against 100%), and a large amphiphile fraction is the main cause of failure there (logistic coefficient −3.3 on standardised constants).
- **The three acidic points (A, B and rep1_9823, pH < 4) cannot hold any membrane** in this model, because the window needs pH ≥ pKa_app − w_pH ≥ 4. This follows from the layer C ranges, which are already more permissive than measured fatty-acid and alcohol mixtures (vesicles from about pH 6.5; Jordan et al. 2019).
- Membrane constants that matter, from logistic fits per point: the amphiphile fraction (#3709, #1558), φ₀ and f_rehyd (more enclosure helps at #1558 and C), and pKa_app and γ_T at C. At rep1_1926 nearly every set succeeds, so the fit there is weak.
- The share of membrane-constant space compatible with each point is small (0.4–18%). Which membrane chemistries are realistic therefore decides much of the answer.

## Focused map box3 (spec §9b)

28 October 2026, 14:55 to 17:31 Toronto time (9,361 s). Result md5 246114b0; confirmations md5 3e414487. Region: the 4c-3 box2 heredity region, restricted to the membrane stability window and an amphiphile fraction ≥ 0.1. Points outside are masked and not simulated.

| | Runs | Masked | E1 | E2 | L2 strict | L2_long | E3 strict | **E3 strict ∩ L2_long** |
|---|---|---|---|---|---|---|---|---|
| Ground truth 500 × 16 | 8,000 | 97.4% | 190 | 112 | 0 | 0 | 64 | **0** |
| COMPASS | 16,000 | 85.5% | 2,113 | 1,020 | 9 | 2 | 340 | **2** |

**Confirmations (16 new runs each):**

| Point | Environment | E3 strict | L2_long | Both | Margin | Copy length | Enclosed share of L2 lineage |
|---|---|---|---|---|---|---|---|
| **engine12538** | T 4.5 °C, T_m8 13.0 °C, period 63 h, pH 7.58, Mg²⁺ 1.0 mM, salinity 0.74 M, τ 11,178 h; k_t 2.8×10⁴, ε 0.0013 | 16 | 16 | **16** | 0.57 | 61.9 | 1.0 |
| engine8869 | T 1.7 °C, pH 6.87, period 10 h, τ 88,210 h | 16 | 0 | 0 | 0.07 | 59.3 | – |

- **engine12538 is the first point found by the search itself, not chosen by hand, where heredity and persistent filled vesicles coexist in every run.** Membrane constants there: amphiphile fraction 0.40, CVC₀ 0.20 mM, pKa_app 8.06 ± 1.10, Mg_crit 29 mM, sal_crit 2.76 M, T_mem 41 °C, φ₀ 0.042, f_rehyd 8.4, β 9.7. It is a cold, near-neutral pool, like #1558 and rep1_1926.
- engine8869 was a single-run fluke for heredity and is not confirmed.
- Only 2.6% of the region is feasible under masking, so the ground truth is too thin for a map of the joint level. The conditional version (§9c) addresses this.

## Focused map box3 v2 (spec §9c)

28 October 2026, 17:40 to 22:55 Toronto time; confirmations until 23:19. Result md5 bd6d27f7, confirmations md5 dba35a53. Region as in box3, but the membrane constants are drawn inside the stability window, given the environment (§9c). Environments where no membrane can be stable stay masked. This is a map conditional on a stable membrane.

| | Runs | Infeasible | E1 | E2 | L1 | L2 strict | L2_long | E3 loose | E3 strict | **E3 strict ∩ L2_long** |
|---|---|---|---|---|---|---|---|---|---|---|
| Ground truth 500 × 16 | 8,000 | 73.4% | 1,756 | 626 | 0 | 4 | 0 | 240 | 208 | **0** |
| COMPASS | 16,000 | 59.6% | 5,647 | 2,258 | 8 | 18 | 5 | 640 | 604 | **5** |

**Confirmations (16 new runs each):**

| Point | E3 strict | L2 strict | **Both** | Margin | Copy length | Enclosed share of L2 lineage |
|---|---|---|---|---|---|---|
| **engine11230** | 16 | 16 | **16** | 0.36 | 62.3 | 1.0 |
| **engine15233** | 16 | 16 | **16** | 0.57 | 61.3 | 1.0 |
| **engine15539** | 16 | 13 | **13** | 0.43 | 58–62 | 0.86 |
| engine8083 | 16 | 7 | 7 | 0.28 | 59.4 | 0.998 |
| engine4527 | 16 | 1 | 1 | – | 58–62 | – |

- **The search found three more points where heredity and filled vesicles coexist in most or all runs** (engine11230 and engine15233 in 16 of 16, engine15539 in 13 of 16). Together with engine12538 from box3, that makes four search-found points. engine8083 is marginal, and engine4527 was a fluke.
- At every confirmed point, E3 strict holds in all 16 runs, so the limit is heredity, not the membrane. Essentially all lineage mass is inside vesicles.
- The ground truth again has no joint hit in 8,000 runs. The joint level is rare even when the membrane is guaranteed stable: COMPASS found 5 candidates among 6,460 feasible runs.

## Next

1. The 4c-4 report.
2. Stage 4c-5 ([spec](../docs/stage4c5_spec.md)): coupling, L3 and L4 detection at these points and in the box3 v2 region.
