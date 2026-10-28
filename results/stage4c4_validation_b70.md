# Stage 4c-4 validation: amphiphiles and vesicles (Intel Arc Pro B70)

**Status:** validation complete (all tests pass), 28 October 2026, 02:08 (Toronto time) · **Specification:** [stage4c4_spec.md](../docs/stage4c4_spec.md), with the implementation choices in its §8 · **Code:** chem4c4.py (v2, md5 b79360f6), validate4c4.py (765d5e41), check3709_4c4.py (b2a60bdf); kept on the benchmark machine and to be added to [stage4/](../stage4/README.md) · **Data:** result JSONs with md5 sums on the benchmark machine, to be added to the repository

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
- E2 is higher than in 4c-3 (0.107 against 0.085). A paired rerun of the same points without amphiphiles is running to find out whether vesicles cause this.

## Next

1. E2 with and without amphiphiles on the same points.
2. Membrane scan at the heredity points (spec §9a). The three acidic points (pH < 4) are outside every membrane window, so the scan uses the other four.
3. Focused map of heredity and membranes together (box3, spec §9b), with confirmations.
4. The 4c-4 report.
