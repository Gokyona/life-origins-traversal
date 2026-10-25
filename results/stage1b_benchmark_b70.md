# Stage 1b benchmark: disagreement-driven ordered testing vs random search (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU, ~90 W, 2800 MHz) · **Wall time:** 24 min · **Code:** [`stage1b/`](../stage1b) at commit `6bb7cd6`

> **Update (stage 1c):** the caveat below was confirmed. With a fair, early-stopping map-maker, all methods improve a lot and the engine's advantage becomes "same accuracy with about 8× fewer simulations" rather than a permanent gap. See [stage1c_benchmark_b70.md](stage1c_benchmark_b70.md).

## Bottom line

**First positive result, with one caveat to check.** On a harder 10-condition problem, the revised engine (`engine3`) beats random search and Sobol sampling at every budget from 4,000 runs upward: about **+1.5 points on all conditions** and **+8 to +12 points on the hard, stochastic conditions**, with small spread across 5 repetitions. The caveat: the baselines barely improve between 2,000 and 256,000 runs, which suggests the shared map-maker, not the data, is limiting them. The advantage is real under this map-maker; whether it survives a map-maker that keeps learning from more data is the next test.

## What changed from stage 1

Stage 1 ([report](stage1_benchmark_b70.md)) was negative. Three changes followed from it:

1. **Harder problem** (`stage1b/model10.py`): 10 conditions (adds parasite fitness, base hydrolysis rate, wet/dry cycling, two distractors). Life fraction 0.108; the hard set (true success probability 0.1–0.4 or 0.6–0.9) has 666 of 8,000 test conditions.
2. **Shared neural map-maker** (`stage1b/surrogate.py`): every method's data is turned into a map by the same small neural network (2×128, 600 training steps), replacing kNN.
3. **New refinement rule** (`engine3`): one run per condition, 30% of each batch along the Sobol sequence, 70% where five bootstrap networks **disagree most**. Disagreement marks missing data; stage 1's p(1−p) rule chased conditions whose outcome is random by nature. `engine3_explore50` uses a 50/50 split.

Ground truth: 8,000 held-out random conditions × 128 runs each (3 s on the B70). Mean ± standard deviation over 5 repetitions. Command: `python bench1b.py --budget 256000 --reps 5 --n-test 8000 --gt-seeds 128`.

## Results

Accuracy on all 8,000 test conditions, by simulation runs spent:

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k | 256k |
|---|---|---|---|---|---|---|---|---|
| random_1seed | 0.944 ± 0.006 | 0.945 ± 0.006 | 0.944 ± 0.006 | 0.944 ± 0.006 | 0.945 ± 0.005 | 0.945 ± 0.005 | 0.945 ± 0.004 | 0.944 ± 0.005 |
| sobol_1seed | 0.943 ± 0.003 | 0.946 ± 0.006 | 0.944 ± 0.004 | 0.944 ± 0.004 | 0.944 ± 0.004 | 0.944 ± 0.006 | 0.944 ± 0.005 | 0.945 ± 0.004 |
| sobol_8seed | 0.936 ± 0.002 | 0.938 ± 0.003 | 0.942 ± 0.003 | 0.947 ± 0.004 | 0.948 ± 0.006 | 0.948 ± 0.006 | 0.948 ± 0.005 | 0.948 ± 0.005 |
| **engine3** | 0.943 ± 0.003 | **0.958 ± 0.001** | **0.960 ± 0.002** | **0.963 ± 0.001** | **0.963 ± 0.002** | **0.964 ± 0.000** | **0.963 ± 0.001** | **0.960 ± 0.002** |
| engine3_explore50 | 0.943 ± 0.003 | **0.958 ± 0.002** | **0.960 ± 0.001** | 0.961 ± 0.002 | 0.961 ± 0.001 | 0.962 ± 0.002 | 0.962 ± 0.001 | **0.960 ± 0.001** |

Accuracy on the 666 hard (stochastic) conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k | 256k |
|---|---|---|---|---|---|---|---|---|
| random_1seed | 0.705 ± 0.038 | 0.711 ± 0.045 | 0.707 ± 0.039 | 0.705 ± 0.039 | 0.705 ± 0.032 | 0.704 ± 0.029 | 0.706 ± 0.027 | 0.702 ± 0.027 |
| sobol_1seed | 0.694 ± 0.026 | 0.712 ± 0.037 | 0.698 ± 0.024 | 0.698 ± 0.022 | 0.699 ± 0.028 | 0.702 ± 0.035 | 0.702 ± 0.030 | 0.703 ± 0.026 |
| sobol_8seed | 0.675 ± 0.012 | 0.675 ± 0.007 | 0.687 ± 0.015 | 0.715 ± 0.021 | 0.724 ± 0.036 | 0.723 ± 0.040 | 0.723 ± 0.036 | 0.721 ± 0.033 |
| **engine3** | 0.694 ± 0.026 | **0.800 ± 0.005** | **0.808 ± 0.008** | **0.819 ± 0.008** | **0.818 ± 0.018** | **0.826 ± 0.005** | **0.820 ± 0.006** | 0.804 ± 0.015 |
| engine3_explore50 | 0.694 ± 0.026 | 0.794 ± 0.015 | 0.805 ± 0.007 | 0.805 ± 0.013 | 0.809 ± 0.010 | 0.810 ± 0.010 | 0.816 ± 0.010 | **0.811 ± 0.013** |

(At 2k all Sobol-based methods share the same first points, so they tie by construction.)

Time per repetition: baselines 12–16 s; engines 120–129 s (17 refinement rounds each, five network fits per round).

## Findings

1. **The engine wins, and consistently.** From 4k runs on, both engine variants beat all three baselines on both metrics. On the hard set the gap (about 0.80–0.83 vs 0.70–0.72) is several times larger than the spread across repetitions; every single engine repetition at 256k (0.779–0.827) is above the best baseline mean.
2. **The engine reaches its level with ~4k runs; the baselines never reach it.** In this setting ordered testing is worth more than a 60-fold increase in random sampling.
3. **Caveat: the baselines are flat from 2k to 256k.** More random data does not improve their maps, so the network fit (fixed 600 steps, fixed size) is the bottleneck. The engine's points, concentrated near the boundary, make the most of that limited map-maker. A fair final verdict needs a map-maker whose training grows with the data.
4. **Engine accuracy peaks around 64k and slips slightly at 256k** (hard set 0.826 → 0.804). Likely the same cause: the last rounds pile very many points into a thin band, and a fixed-budget fit then under-weights the rest of the space. The 50/50 explore variant degrades less.
5. **Explore fraction matters little** (30% vs 50%): 30% is slightly better at medium budgets, 50% slightly more stable at the largest.

## Next steps

- Re-run with a map-maker whose training steps scale with the number of points (and a larger network), to confirm finding 1 is not an artefact of finding 3.
- Cap or anneal refinement in late rounds to avoid the 256k dip.
- Stage 2: replace the well-mixed chemostat with a spatial protocell model (membranes, growth, division), where each run is far more expensive and fewer, better-placed tests matter even more.

## Raw data

`results_1b.json` and `bench_1b.log` are on the benchmark machine and will be added to this folder.
