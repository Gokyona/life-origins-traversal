# Stage 3e: joint nested map and comparison with a standard active-learning rule (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Joint map run:** 19:26–19:38 Toronto time (12 min) · **Straddle comparison:** 19:38–20:32 Toronto time (55 min) · **Plan:** [stage3_plan.md](../docs/stage3_plan.md) · **Code:** [stage3/s3e_joint.py](../stage3/s3e_joint.py), [s3e_straddle.py](../stage3/s3e_straddle.py)

## Bottom line

Two engine questions were left open before the environmental search of Part II.

1. **Can one map cover all four levels of life without contradicting itself?** Yes. A single network whose outputs are nested by construction never predicts a higher level without the lower ones. Four separately trained maps do this at 1.1–1.4% of conditions. The joint map is as accurate (equal within noise at every level and budget, slightly better for L4) and trains about 3× faster. It becomes the default map-maker for multi-level maps.
2. **Is COMPASS's selection rule better than a standard active-learning rule for finding boundaries?** No, and not worse either. The straddle rule of Bryan et al. (2005), the standard heuristic for locating where a function crosses a threshold, gives the same accuracy as COMPASS's disagreement rule at every budget from 4k to 64k. All four active rules need about 4× fewer simulations than Sobol sampling to reach the same overall accuracy. COMPASS's advantage therefore comes from the loop as a whole (ensemble map-maker, rounds, exploration share, no repeats), not from its particular scoring formula. We say this plainly in the stage 3 paper.

## (a) Joint nested map

**Design.** One network (11 inputs → 256 → 256 → 4 outputs z1…z4) with a product link:

P(≥ L1) = s(z1), P(≥ L2) = P(≥ L1)·s(z2), P(≥ L3) = P(≥ L2)·s(z3), P(≥ L4) = P(≥ L3)·s(z4),

where s is the logistic function. Each probability is at most the one below it, so the map is nested whatever the weights are. The loss is the sum of the four binary cross-entropies, computed in log space, on all four labels of every simulation. Training uses the same early stopping as the separate maps (10% validation, patience 5, at most 10,000 steps). The comparison is four separate maps (the stage 1–3 map-maker, hidden width 256) trained on the same simulations.

**Setting.** The evolving protocell model of [stage 3c](stage3c_evolution_b70.md) (11 dimensions, L1–L4). Ground truth: 4,000 conditions × 64 seeds (life fractions L1 0.558, L2 0.338, L3 0.173, L4 0.112). Training data: Sobol points with one simulation each, 3 repetitions per budget.

**Results** (mean of 3 repetitions; overall accuracy / hard-set accuracy; hard sets L1 317, L2 118, L3 145, L4 453 points; "violations" = test conditions where a higher level is predicted without a lower one):

| budget | map | L1 | L2 | L3 | L4 | violations | training |
|---|---|---|---|---|---|---|---|
| 8k | separate | 0.937 / 0.68 | 0.959 / 0.60 | 0.947 / 0.63 | 0.953 / 0.76 | 1.13% | 17 s |
| | joint | 0.935 / 0.68 | 0.957 / 0.60 | 0.950 / 0.60 | 0.956 / 0.79 | **0** | 4 s |
| 16k | separate | 0.948 / 0.73 | 0.963 / 0.60 | 0.951 / 0.64 | 0.959 / 0.80 | 1.39% | 21 s |
| | joint | 0.945 / 0.71 | 0.962 / 0.62 | 0.956 / 0.63 | 0.963 / 0.83 | **0** | 6 s |
| 32k | separate | 0.954 / 0.74 | 0.968 / 0.66 | 0.962 / 0.70 | 0.962 / 0.82 | 1.39% | 30 s |
| | joint | 0.954 / 0.74 | 0.971 / 0.69 | 0.964 / 0.68 | 0.966 / 0.84 | **0** | 10 s |

**On COMPASS-chosen data.** In practice maps are trained on points the engine chose, not on Sobol points. With 32,000 points chosen by COMPASS for L4 (one repetition), separate / joint: L1 0.946 / 0.946, L2 0.978 / 0.972, L3 0.967 / **0.974**, L4 0.968 / **0.971**; violations 1.2% / 0. The joint map again removes all contradictions. It gains on the levels the points were chosen for (L3 hard 0.67 → 0.77) and loses a little on L2, which the points were not chosen for.

## (b) COMPASS against the straddle rule

**Design.** The same engine loop as the stage 1d/2 default (fast mode, 5-member bootstrap ensemble, 2,048 Sobol starting points, budget × 4/3 per round, 30% of each round Sobol exploration). Only the rule that picks the other 70% differs:

| method | selection rule |
|---|---|
| engine3 | the points where the ensemble disagrees most (top-k of sd) |
| engine4_soft (COMPASS default) | sample points with probability ∝ sd² |
| straddle | top-k of 1.96·sd − \|mean − 0.5\| (Bryan et al. 2005) |
| straddle_soft | sample with probability ∝ max(0, straddle)² |
| sobol_1seed | no selection: Sobol points, one simulation each |

The straddle score is high for points that are both uncertain and predicted near the boundary (probability 0.5). It is the standard heuristic for threshold (level-set) estimation and a natural baseline for COMPASS.

**Setting.** The coupled protocell model of [stage 2b](stage2b_protocell_benchmark_b70.md), level L3, with the same ground truth (4,000 conditions × 64 seeds) and seeds as the stage 2b benchmark; hard set 107 points. Budget 64,000, 5 repetitions.

**Overall accuracy** (mean ± sd over 5 repetitions):

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| sobol_1seed | 0.947 ± 0.002 | 0.953 ± 0.002 | 0.955 ± 0.001 | 0.961 ± 0.002 | 0.969 ± 0.002 | 0.976 ± 0.002 |
| engine3 | 0.946 ± 0.002 | 0.956 ± 0.002 | 0.966 ± 0.001 | 0.975 ± 0.001 | 0.981 ± 0.001 | 0.986 ± 0.001 |
| engine4_soft | 0.946 ± 0.002 | 0.959 ± 0.002 | 0.964 ± 0.002 | 0.977 ± 0.003 | 0.983 ± 0.001 | 0.987 ± 0.001 |
| straddle | 0.946 ± 0.002 | 0.957 ± 0.003 | 0.965 ± 0.002 | 0.975 ± 0.001 | 0.983 ± 0.001 | 0.987 ± 0.001 |
| straddle_soft | 0.946 ± 0.002 | 0.959 ± 0.001 | 0.966 ± 0.002 | 0.976 ± 0.000 | 0.983 ± 0.002 | 0.987 ± 0.001 |

**Hard-set accuracy** (107 points near the boundary):

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| sobol_1seed | 0.598 ± 0.026 | 0.622 ± 0.051 | 0.598 ± 0.017 | 0.626 ± 0.017 | 0.645 ± 0.031 | 0.699 ± 0.042 |
| engine3 | 0.574 ± 0.020 | 0.609 ± 0.032 | 0.626 ± 0.025 | 0.675 ± 0.021 | 0.731 ± 0.030 | 0.793 ± 0.023 |
| engine4_soft | 0.574 ± 0.020 | 0.626 ± 0.013 | 0.624 ± 0.011 | 0.692 ± 0.031 | 0.744 ± 0.017 | 0.772 ± 0.038 |
| straddle | 0.574 ± 0.020 | 0.641 ± 0.039 | 0.654 ± 0.016 | 0.664 ± 0.026 | 0.736 ± 0.023 | 0.776 ± 0.026 |
| straddle_soft | 0.574 ± 0.020 | 0.621 ± 0.019 | 0.639 ± 0.039 | 0.677 ± 0.024 | 0.748 ± 0.037 | 0.791 ± 0.025 |

**Compute per repetition** (simulation + selection, without scoring): Sobol 44 s; each of the four active rules 122–123 s, i.e. about 78 s of engine overhead. The selection rules themselves cost the same; the overhead is the ensemble training.

**Reading.** All four active rules reach at 16k the overall accuracy Sobol reaches at 64k (0.975–0.977 vs 0.976), about **4× fewer simulations**, as in stage 2b. Between the active rules the differences are at most 0.002 overall and within one standard deviation on the hard set at every budget (the hard set has only 107 points, so one point is 0.9 percentage points). Straddle is slightly ahead on the hard set at 4k–8k and the soft disagreement rule at 16k; neither lead is significant.

## What this means

1. **The joint nested map is adopted.** Maps for several levels of life are now made by one network that cannot contradict the definition of the levels (L4 implies L3 implies L2 implies L1), at no cost in accuracy and a third of the training time.
2. **The scoring formula is not COMPASS's contribution.** Disagreement and straddle are interchangeable here. We keep the soft disagreement rule as default (it needs no threshold and extends unchanged to probabilities that are not centred on 0.5), and we describe COMPASS as a combination of known components whose value is the whole loop and its validation, not a new acquisition function.
3. **Part I's engine questions are settled.** The remaining limits (late plateau near 0.96–0.99, small hard sets, models that presuppose cells) are written up in the stage 3 paper.

## Limits

- One model per question: the joint map was tested on the 3c model, the straddle comparison on the 2b model at one level.
- Three (a) and five (b) repetitions; hard sets of 107–453 points. Differences below about two points of hard accuracy cannot be resolved.
- Only one COMPASS-chosen data set for (a).
- Other active-learning methods (Gaussian-process level-set methods, BALD) were not compared; straddle was chosen because it is the established threshold heuristic and fits the same ensemble without a new model.

## Reference

Bryan, B., Schneider, J., Nichol, R. C., Miller, C. J., Genovese, C. R. & Wasserman, L. (2005). Active learning for identifying function threshold boundaries. *Advances in Neural Information Processing Systems 18 (NIPS 2005)*, 163–170. [proceedings](https://proceedings.neurips.cc/paper/2005/hash/8e930496927757aac0dbd2438cb3f4f6-Abstract.html) (straddle score as defined there: 1.96·σ̂ − |f̂ − t|).

## Raw data

`s3e_joint_results.json` and `s3e_straddle_results.json` (in [stage3/](../stage3/)) hold every repetition. Logs are on the benchmark machine and will be uploaded with the other raw results.
