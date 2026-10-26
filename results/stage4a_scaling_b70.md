# Stage 4a: does COMPASS still work in 20–30 dimensions? (Intel Arc Pro B70)

**Date:** 2026-10-25/26 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Run:** 23:59–00:46 Toronto time (47 min) · **Plan:** [stage4_plan.md](../docs/stage4_plan.md), step 4a · **Code:** [stage4/s4a_scaling.py](../stage4/s4a_scaling.py) · **Data:** [stage4/s4a_scaling_results.json](../stage4/s4a_scaling_results.json)

## Bottom line

The stage 4 environment space has about 20–30 dimensions: the 13 environmental variables, several of which need more than one number, plus uncertain constants. COMPASS had been tested up to 12. On a test model padded to 20 and 30 dimensions, COMPASS still beats Sobol sampling at every budget from 4k on.

- **Savings:** about 3.6–4.0× fewer simulations when the extra dimensions are irrelevant, and about 2.8–2.9× when they all have small real effects. The 10-dimensional original gives 4.4×.
- **Overhead:** the engine's cost does not grow with dimension, at about 75 s per 64k-simulation run.
- **No change needed:** no screening step is needed before the stage 4 maps.

## Design

The model is the coupled protocell model of [stage 2b](stage2b_protocell_benchmark_b70.md) at level L3, with 10 inputs. Extra inputs are appended to reach D = 20 or 30, in two kinds:

- **Distractor:** the model ignores the extra inputs.
- **Weak:** each extra input j multiplies the food supply by 1 − 0.15·(u_j − 0.5)·c_j, with fixed random signs c_j. One input moves food by at most ±7.5%. Together they spread importance over many dimensions, as in a real environment map where many variables matter a little.

The first 10 coordinates of the test conditions are the stage 2b test set. Ground truth is 4,000 conditions × 64 simulations per configuration. Methods are Sobol sampling with one simulation per point, and COMPASS (`engine4_soft`, fast mode), with 3 repetitions each and a budget of 64,000. The D = 10 figures are repetitions 0–2 of the [stage 3e](stage3e_engine_b70.md) run, which has the same seeds and ground truth.

## Results

**Accuracy at 64k simulations** (mean ± s.d. over 3 repetitions):

| configuration | hard set | Sobol, all | COMPASS, all | Sobol, hard | COMPASS, hard |
|---|---|---|---|---|---|
| D = 10 (original) | 107 | 0.975 ± 0.001 | **0.987 ± 0.001** | 0.667 ± 0.016 | **0.779 ± 0.023** |
| D = 20, distractors | 107 | 0.967 ± 0.001 | **0.980 ± 0.001** | 0.617 ± 0.013 | **0.707 ± 0.012** |
| D = 20, weak | 114 | 0.971 ± 0.001 | **0.980 ± 0.001** | 0.652 ± 0.022 | **0.722 ± 0.041** |
| D = 30, distractors | 107 | 0.964 ± 0.002 | **0.979 ± 0.001** | 0.645 ± 0.040 | **0.701 ± 0.020** |
| D = 30, weak | 109 | 0.967 ± 0.002 | **0.979 ± 0.000** | 0.654 ± 0.038 | **0.709 ± 0.004** |

**Overall accuracy along the budget** (COMPASS / Sobol, mean of 3):

| configuration | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| D = 10 | .947 / .948 | .960 / .954 | .964 / .956 | .976 / .959 | .983 / .968 | .987 / .975 |
| D = 20, distractors | .926 / .933 | .951 / .945 | .959 / .956 | .967 / .960 | .972 / .965 | .980 / .967 |
| D = 20, weak | .924 / .935 | .954 / .945 | .956 / .953 | .965 / .957 | .976 / .964 | .980 / .971 |
| D = 30, distractors | .899 / .899 | .942 / .931 | .955 / .949 | .963 / .954 | .973 / .958 | .979 / .964 |
| D = 30, weak | .898 / .900 | .946 / .933 | .957 / .950 | .963 / .956 | .970 / .959 | .979 / .967 |

**Budget at which COMPASS reaches Sobol's 64k overall accuracy** (log-linear interpolation): 14.6k at D = 10, 15.8k and 17.6k with distractors at D = 20 and 30, and 22.3k and 23.1k with weak dimensions. That is about 4.4×, 4.0×, 3.6×, 2.9× and 2.8× fewer simulations. On the hard set, Sobol's curve is not monotone in some configurations, so the equivalent hard-set figures are not reliable.

**Compute per repetition** (simulation + selection): Sobol 42–43 s; COMPASS 114–121 s. The engine overhead is 72–78 s at every dimension, against 78 s at D = 10.

## Reading

1. **Both methods lose about one point of overall accuracy** going from 10 to 20–30 dimensions. That loss comes from the map-maker having to learn a wider input, not from the selection rule.
2. **COMPASS keeps its lead** in overall accuracy (+0.009 to +0.015, against +0.012 at D = 10). Its lead on the hard set shrinks from +0.11 to +0.06–0.09.
3. **The saving is smaller when importance is spread over many weak dimensions** (about 2.8× against 3.6–4.0×). Real environment maps will probably lie between the two cases. We plan with a saving of about 3×.
4. **At 2k (the Sobol start) COMPASS is slightly behind at D = 20** (0.926 vs 0.933): in one repetition the first map was poor. It catches up by 4k.
5. **We did not measure how many outcomes the weak dimensions actually change.** Life fraction and hard-set size barely move (0.161–0.162), so their effect is small, as intended.

## Consequences for stage 4

COMPASS is used unchanged for the stage 4 maps. With bottom-up simulations expected to cost 0.1–10 s each, a 3× saving at 20–30 dimensions means days rather than weeks of machine time per map. A screening step (importance on a cheap first map, then refinement in the important dimensions) is kept in reserve in case the real maps behave worse than this test.
