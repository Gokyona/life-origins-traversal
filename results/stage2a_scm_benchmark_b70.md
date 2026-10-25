# Stage 2a benchmark: COMPASS on a stochastic corrector protocell population (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Run:** 14:41–15:03 Toronto time (22 min) · **Code:** [`stage2/`](../stage2) with the stage 1b harness, run with `--fast`

## Bottom line

**COMPASS reaches with 8,000 simulations the accuracy that random or Sobol sampling reaches with 64,000: about 8× fewer simulations.** At equal budget (64k) it makes about half as many mistakes (0.7% vs 1.3% of test conditions wrong). The saving matches the stage 1 protocell model (about 8×) and is larger than on the RAF model (about 2×), consistent with the earlier finding that the engine gains most when the life boundary is thin.

## The model

A well-mixed population of 128 protocells after the stochastic corrector model of Szathmáry & Demeter (1987); see [`stage2/README.md`](../stage2/README.md). Eight conditions: division size, parasite advantage, mutation rate, number of essential gene types, degradation, food, and two distractors. The target is level L2 of the life detector: at least half of the cells still functional after at least 15 divisions per cell slot. Before the benchmark the model reproduced the known qualitative behaviour: parasites win in large cells, very small cells lose genes at division, and there is a critical mutation rate.

Ground truth: 4,000 random conditions × 64 simulations (147 s). Life fraction 0.208. Hard set (success probability 0.1–0.4 or 0.6–0.9): only 105 conditions (2.6%), because a population of 128 cells averages out much of the randomness. Hard-set numbers below are therefore noisier than in earlier stages.

Command: `python bench_scm.py --device xpu --surrogate es --fast --budget 64000 --reps 5 --n-test 4000 --gt-seeds 64 --methods random_1seed,sobol_1seed,engine4_soft --out results_scm.json`.

## Results

Accuracy on all 4,000 test conditions (mean ± standard deviation over 5 repetitions):

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| random_1seed | 0.961 ± 0.003 | 0.968 ± 0.005 | 0.976 ± 0.001 | 0.981 ± 0.001 | 0.985 ± 0.001 | 0.987 ± 0.001 |
| sobol_1seed | 0.960 ± 0.003 | 0.973 ± 0.001 | 0.976 ± 0.001 | 0.981 ± 0.001 | 0.984 ± 0.001 | 0.987 ± 0.001 |
| **COMPASS (engine4_soft)** | **0.964 ± 0.004** | **0.980 ± 0.002** | **0.987 ± 0.001** | **0.990 ± 0.001** | **0.992 ± 0.001** | **0.993 ± 0.001** |

Accuracy on the 105 hard conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| random_1seed | 0.630 ± 0.046 | 0.690 ± 0.045 | 0.716 ± 0.044 | 0.758 ± 0.029 | 0.789 ± 0.040 | 0.792 ± 0.027 |
| sobol_1seed | 0.632 ± 0.038 | 0.701 ± 0.060 | 0.724 ± 0.029 | 0.728 ± 0.021 | 0.770 ± 0.025 | 0.789 ± 0.045 |
| **COMPASS (engine4_soft)** | 0.632 ± 0.026 | **0.722 ± 0.036** | **0.806 ± 0.027** | **0.802 ± 0.029** | **0.836 ± 0.013** | **0.859 ± 0.011** |

At 64k per repetition (all / hard): COMPASS 0.994/0.848, 0.992/0.857, 0.993/0.876, 0.994/0.867, 0.993/0.848; Sobol 0.988/0.829, 0.986/0.743, 0.988/0.848, 0.988/0.790, 0.987/0.733.

Time per repetition: baselines 47–58 s; COMPASS 134–140 s (12 refinement rounds).

## Findings

1. **About 8× fewer simulations.** COMPASS at 8k (0.987) equals both baselines at 64k (0.987) on all conditions; on the hard set it is already ahead at 8k (0.806 vs 0.789–0.792 at 64k).
2. **Half the errors at equal budget.** At 64k, 0.7% of test conditions are misclassified by COMPASS versus 1.3% by the baselines.
3. **More stable.** The hard-set spread across repetitions at 64k is ±0.011 for COMPASS versus ±0.027 and ±0.045 for the baselines.
4. **The saving tracks the shape of the boundary.** Across the three models tested so far: protocell chemostat about 8× (hard set 8% of the space), RAF about 2× (23%), stochastic corrector about 8× (2.6%). The thinner the uncertain region, the more it pays to concentrate tests on it.

## Next

Stage 2b: make division follow from the cell's own chemistry (membrane made inside the cell) and add the L3 detector.

## Raw data

`results_scm.json` and `bench_scm.log` are on the benchmark machine and will be added to this folder.
