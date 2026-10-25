# Stage 2b benchmark: COMPASS on coupled protocells, level L3 (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Run:** 15:17–15:39 Toronto time (22 min) · **Code:** [`stage2/protocell.py`](../stage2/protocell.py) with the stage 1b harness, run with `--fast`

## Bottom line

**On the most complex model so far, COMPASS reaches with 16,000 simulations the accuracy that random or Sobol sampling reaches with 64,000: about 4× fewer simulations.** At equal budget (64k) it misclassifies 1.4% of test conditions against 2.4–2.5% for the baselines, about 40% fewer errors. On the hardest conditions the baselines barely improve with more simulations, while COMPASS keeps improving.

## The model

A well-mixed population of 128 stochastic-corrector protocells (as in stage 2a) in which **division follows from the cell's own chemistry**: one gene type makes membrane from food, the membrane sets how much a cell can hold, a cell bursts when overfilled, and a cell divides only when its membrane has doubled. Ten conditions: division size, parasite advantage, mutation rate, number of gene types, degradation, food, membrane synthesis rate, membrane turnover, and two distractors. The target is level **L3** of the life detector: L2 (at least half the cells functional after at least 15 divisions per cell slot) plus homeostasis (content in balance with the membrane, few bursts). Validation is in [`stage2/README.md`](../stage2/README.md).

Ground truth: 4,000 random conditions × 64 simulations (151 s). Life fraction (L3) 0.160. Hard set (success probability 0.1–0.4 or 0.6–0.9): 103 conditions (2.6%), so hard-set numbers are noisy.

Command: `python bench_protocell.py --device xpu --surrogate es --fast --budget 64000 --reps 5 --n-test 4000 --gt-seeds 64 --methods random_1seed,sobol_1seed,engine4_soft --out results_protocell.json`.

## Results

Accuracy on all 4,000 test conditions (mean ± standard deviation over 5 repetitions):

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| random_1seed | 0.947 ± 0.003 | 0.954 ± 0.002 | 0.958 ± 0.002 | 0.962 ± 0.003 | 0.969 ± 0.003 | 0.976 ± 0.002 |
| sobol_1seed | 0.945 ± 0.002 | 0.953 ± 0.004 | 0.959 ± 0.001 | 0.962 ± 0.001 | 0.969 ± 0.001 | 0.975 ± 0.002 |
| **COMPASS (engine4_soft)** | 0.946 ± 0.002 | **0.959 ± 0.002** | **0.964 ± 0.002** | **0.977 ± 0.003** | **0.983 ± 0.001** | **0.986 ± 0.001** |

Accuracy on the 103 hard conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| random_1seed | 0.619 ± 0.039 | 0.649 ± 0.037 | 0.664 ± 0.010 | 0.645 ± 0.029 | 0.678 ± 0.042 | 0.666 ± 0.030 |
| sobol_1seed | 0.621 ± 0.014 | 0.625 ± 0.030 | 0.670 ± 0.031 | 0.662 ± 0.033 | 0.668 ± 0.030 | 0.703 ± 0.039 |
| **COMPASS (engine4_soft)** | 0.606 ± 0.016 | **0.660 ± 0.025** | 0.666 ± 0.010 | **0.705 ± 0.055** | **0.750 ± 0.032** | **0.790 ± 0.035** |

At 64k per repetition (all / hard): COMPASS 0.986/0.805, 0.985/0.776, 0.986/0.825, 0.985/0.815, 0.984/0.728; Sobol 0.978/0.737, 0.974/0.650, 0.972/0.660, 0.973/0.737, 0.976/0.728; random 0.976/0.621, 0.975/0.699, 0.977/0.689, 0.977/0.679, 0.973/0.640.

Time per repetition: baselines 50–59 s; COMPASS 140–146 s (12 refinement rounds).

## Findings

1. **About 4× fewer simulations.** COMPASS at 16k (0.977) matches both baselines at 64k (0.975–0.976) on all conditions.
2. **About 40% fewer errors at equal budget.** At 64k: 1.4% of conditions wrong for COMPASS, 2.4–2.5% for the baselines.
3. **Only COMPASS keeps learning the ambiguous regions.** On the hard set the baselines go from 0.62 to 0.67–0.70 between 2k and 64k; COMPASS goes from 0.61 to 0.79.
4. **A harder boundary.** All methods are less accurate than on the stage 2a model (0.986 vs 0.993 for COMPASS at 64k): the ten-dimensional L3 boundary, which includes a window in membrane synthesis rate (too slow and too fast both fail), is harder to draw. Savings so far at equal accuracy: protocell chemostat about 8×, RAF about 2×, stochastic corrector (2a) about 8×, coupled protocell (2b, L3) about 4×; every one is above the stage 2 target of 2×.

## Next

Stage 2d/2e: multi-fidelity COMPASS (cheap short runs to screen, long runs to confirm near the boundary) and maps for L1, L2 and L3 from the same simulations, to compare which conditions each level of life needs.

## Raw data

`results_protocell.json` and `bench_protocell.log` are on the benchmark machine and will be added to this folder.
