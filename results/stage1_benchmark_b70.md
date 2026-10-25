# Stage 1 benchmark: ordered testing vs random search (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Wall time:** 38 min · **Code:** [`stage1/`](../stage1) at commit `cf9ed7c`

## Bottom line

**This run does not show that the ordered engine beats random search.** On this test problem, the best methods at every budget are the simplest ones: one seed per condition, placed either at random or along a Sobol sequence. The boundary-refinement engine was the worst method from 4,000 runs upward. We publish the result as it is, together with what it teaches us about the next design.

## Setup

- **Model:** the stage-0 stochastic protocell chemostat (`stage1/model.py`): 6 conditions (Da, μL, food inflow, spontaneous formation rate, temperature, and a distractor with no effect), 800 time steps per run. A run counts as "life" when the master replicator is amplified far above the abiotic background.
- **Ground truth:** 4,000 random held-out conditions × 128 seeds each, computed by brute force (3 s on the B70). Life fraction 0.186. The **hard set** is the 620 conditions whose true success probability is stochastic but decidable (0.1–0.4 or 0.6–0.9).
- **Metric:** every method's tested points are turned into a map with the same surrogate (distance-weighted kNN). Accuracy = share of held-out conditions classified correctly (life if p > 0.5), as a function of the number of stochastic simulation runs spent. Mean of 5 repetitions, budget 128,000 runs.
- **Command:** `python bench.py --device xpu --budget 128000 --reps 5 --n-test 4000 --gt-seeds 128`

## Methods

| method | where points go | seeds per condition |
|---|---|---|
| random_1seed | uniform random | 1 |
| random_8seed | uniform random | 8 |
| sobol_1seed | Sobol sequence | 1 |
| sobol_8seed | Sobol sequence | 8 |
| sobol_sprt | Sobol sequence | sequential test, 4–24 |
| engine_no_monotone | Sobol + boundary refinement (active learning) | sequential test, 4–24 |
| engine_full | as above + monotone inference | sequential test, 4–24 |
| engine2_no_monotone | Sobol + boundary refinement | 1 |
| engine2_full | as above + monotone inference | 1 |

## Results

Accuracy on all 4,000 test conditions, by simulation runs spent:

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k |
|---|---|---|---|---|---|---|---|
| random_1seed | 0.940 | 0.948 | 0.953 | 0.958 | 0.961 | 0.965 | 0.966 |
| random_8seed | 0.911 | 0.928 | 0.941 | 0.949 | 0.955 | 0.961 | 0.967 |
| sobol_1seed | **0.946** | **0.948** | **0.955** | **0.960** | **0.963** | **0.966** | 0.968 |
| sobol_8seed | 0.919 | 0.938 | 0.947 | 0.954 | 0.960 | 0.965 | 0.969 |
| sobol_sprt | 0.918 | 0.928 | 0.940 | 0.947 | 0.956 | 0.961 | 0.967 |
| engine_no_monotone | 0.918 | 0.892 | 0.929 | 0.941 | 0.948 | 0.956 | 0.960 |
| engine_full | 0.918 | 0.894 | 0.933 | 0.943 | 0.945 | 0.953 | 0.962 |
| engine2_no_monotone | 0.935 | 0.938 | 0.947 | 0.952 | 0.958 | 0.965 | **0.970** |
| engine2_full | 0.935 | 0.938 | 0.947 | 0.952 | 0.958 | 0.965 | **0.970** |

Accuracy on the 620 hard (stochastic) conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k |
|---|---|---|---|---|---|---|---|
| random_1seed | 0.797 | 0.815 | 0.835 | 0.856 | 0.857 | 0.878 | 0.885 |
| random_8seed | 0.752 | 0.770 | 0.818 | 0.833 | 0.854 | 0.876 | 0.896 |
| sobol_1seed | **0.816** | **0.823** | **0.841** | **0.864** | **0.870** | 0.880 | 0.892 |
| sobol_8seed | 0.750 | 0.798 | 0.824 | 0.855 | **0.870** | **0.885** | **0.907** |
| sobol_sprt | 0.759 | 0.772 | 0.800 | 0.824 | 0.856 | 0.875 | 0.905 |
| engine_no_monotone | 0.759 | 0.680 | 0.781 | 0.808 | 0.839 | 0.853 | 0.877 |
| engine_full | 0.759 | 0.683 | 0.789 | 0.822 | 0.828 | 0.859 | 0.883 |
| engine2_no_monotone | 0.782 | 0.784 | 0.805 | 0.825 | 0.847 | 0.870 | 0.896 |
| engine2_full | 0.782 | 0.784 | 0.805 | 0.825 | 0.847 | 0.870 | 0.896 |

## Findings

1. **Sobol ordering helps a little; that part of the idea holds.** Sobol beats uniform random at equal seeds almost everywhere (e.g. hard set at 128k: 0.907 vs 0.896 with 8 seeds; 0.892 vs 0.885 with 1 seed). Differences are around one point and partly within rep-to-rep noise.
2. **Spending few seeds per condition and covering more conditions wins at small and medium budgets.** At 2k–16k runs, one seed per condition beats every multi-seed method. Repeated seeds only pay off at the largest budgets on the hard set.
3. **Boundary refinement, as implemented, hurts.** Both engine variants that refine the boundary are below plain Sobol at almost every budget, with a sharp dip at 4k runs (hard set 0.68). Our working explanation: the kNN uncertainty score sends points into regions that are *intrinsically* random (true p near 0.5), where extra samples cannot sharpen the map, while the rest of the space loses coverage. This needs to be checked by looking at where the engine places its points.
4. **Monotone inference is safe but not useful here.** In `engine_full` it labelled about 2,200 conditions per repetition without simulating them, saving about 11k runs (145k vs 156k spent); spot checks found 0–3 contradictions (one repetition hit the limit of 3 and switched inference off). Accuracy did not improve. With single-seed results (`engine2_full`) it never fired, so `engine2_full` equals `engine2_no_monotone`.
5. **The test problem is probably too easy.** Every method reaches about 0.94–0.97 overall accuracy within a few thousand runs, so there is little room for any search strategy to show an advantage. The shared kNN surrogate may also be the limiting factor.

## Next steps

- Diagnose finding 3: plot where boundary refinement places points versus the true boundary.
- Replace the uncertainty score with one that targets *where the decision is uncertain because of missing data*, not where the outcome is intrinsically random (e.g. disagreement between bootstrap surrogates, or expected change in the map).
- Build a harder benchmark: more dimensions (8–12), a thinner and curved life region, and a smaller life fraction, closer to what the real question looks like.
- Try a stronger shared surrogate (Gaussian-process classifier or small neural network) so that all methods are limited by the data they collect, not by the map-maker.

## Raw data

Full learning curves (`results_b70.json`) and the run log (`bench_b70.log`) are on the benchmark machine and will be added to this folder.
