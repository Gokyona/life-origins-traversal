# Stage 1c benchmark: re-checking stage 1b with a fair map-maker (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU, ~91–93 W, 2800 MHz) · **Wall time:** 1 h 38 min · **Code:** [`stage1b/`](../stage1b) at commit `f115c5a`

## Bottom line

**The engine's advantage is real but smaller than stage 1b suggested, and it is an advantage in cost, not in final accuracy.** With a map-maker that keeps learning as data grows, every method improves a lot. The engine (`engine3`) reaches with **32,000 runs** the accuracy that the best baseline reaches with **256,000 runs** (0.989 / 0.947 vs 0.990 / 0.948 on all / hard conditions), so it needs about **8× fewer simulations** for the same map. At very large budgets the baselines catch up and the methods nearly tie.

This corrects stage 1b: its large gap (+8 to +12 points on hard conditions at every budget) was partly an artefact of an under-trained map-maker.

## Why stage 1c

In stage 1b ([report](stage1b_benchmark_b70.md)) the baselines stayed flat from 2k to 256k runs, which pointed at the shared map-maker. Diagnosis confirmed it: the fixed 600-step network fit underfits. With longer training or early stopping, `random_1seed` at 64k already reaches about 0.98 overall and 0.87–0.93 on the hard set, instead of 0.945 / 0.70.

Stage 1c changes only the map-maker (`--surrogate es`, `fit_es` in `stage1b/surrogate.py`):

- Final maps for scoring: 2×256 network, up to 10,000 steps, 10% of the data held out for validation, best checkpoint kept, stop after 5 evaluations without improvement.
- Engine's bootstrap ensemble (used only to choose points): 2×128 network, up to 3,000 steps, same early stopping.

Everything else is identical to stage 1b: the same 10-condition model (life fraction 0.108), the same 8,000 held-out test conditions × 128 runs (666 hard conditions with true success probability 0.1–0.4 or 0.6–0.9), the same methods, 5 repetitions. Command: `python bench1b.py --device xpu --surrogate es --budget 256000 --reps 5 --n-test 8000 --gt-seeds 128`.

## Results

Accuracy on all 8,000 test conditions, by simulation runs spent (mean ± standard deviation over 5 repetitions):

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k | 256k |
|---|---|---|---|---|---|---|---|---|
| random_1seed | 0.948 ± 0.007 | 0.962 ± 0.003 | 0.969 ± 0.004 | 0.975 ± 0.002 | 0.978 ± 0.001 | 0.983 ± 0.001 | 0.986 ± 0.001 | 0.989 ± 0.001 |
| sobol_1seed | 0.953 ± 0.008 | 0.963 ± 0.002 | 0.969 ± 0.001 | 0.975 ± 0.002 | 0.980 ± 0.001 | 0.983 ± 0.001 | 0.987 ± 0.002 | 0.990 ± 0.001 |
| sobol_8seed | 0.936 ± 0.006 | 0.947 ± 0.006 | 0.958 ± 0.002 | 0.964 ± 0.003 | 0.973 ± 0.002 | 0.978 ± 0.002 | 0.982 ± 0.001 | 0.987 ± 0.001 |
| **engine3** | 0.953 ± 0.008 | **0.970 ± 0.002** | **0.978 ± 0.002** | **0.986 ± 0.001** | **0.989 ± 0.002** | **0.991 ± 0.001** | **0.991 ± 0.001** | **0.991 ± 0.000** |
| engine3_explore50 | 0.953 ± 0.008 | 0.968 ± 0.002 | 0.975 ± 0.001 | 0.983 ± 0.000 | 0.988 ± 0.001 | 0.989 ± 0.001 | **0.991 ± 0.001** | 0.990 ± 0.001 |

Accuracy on the 666 hard (stochastic) conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k | 256k |
|---|---|---|---|---|---|---|---|---|
| random_1seed | 0.734 ± 0.041 | 0.794 ± 0.030 | 0.815 ± 0.026 | 0.839 ± 0.012 | 0.870 ± 0.012 | 0.897 ± 0.008 | 0.922 ± 0.009 | 0.944 ± 0.002 |
| sobol_1seed | 0.753 ± 0.051 | 0.803 ± 0.019 | 0.822 ± 0.015 | 0.844 ± 0.019 | 0.871 ± 0.011 | 0.894 ± 0.011 | 0.924 ± 0.009 | 0.948 ± 0.005 |
| sobol_8seed | 0.700 ± 0.029 | 0.747 ± 0.043 | 0.798 ± 0.010 | 0.810 ± 0.023 | 0.848 ± 0.013 | 0.864 ± 0.014 | 0.893 ± 0.008 | 0.926 ± 0.010 |
| **engine3** | 0.753 ± 0.051 | 0.821 ± 0.020 | **0.859 ± 0.017** | **0.920 ± 0.006** | **0.947 ± 0.011** | **0.962 ± 0.004** | **0.962 ± 0.004** | **0.957 ± 0.006** |
| engine3_explore50 | 0.753 ± 0.051 | **0.824 ± 0.021** | 0.845 ± 0.010 | 0.899 ± 0.004 | 0.930 ± 0.005 | 0.946 ± 0.005 | 0.960 ± 0.007 | 0.955 ± 0.004 |

(At 2k all Sobol-based methods share the same first points, so they tie by construction.)

Engine accuracy at 256k per repetition (all / hard): `engine3` 0.991/0.965, 0.990/0.959, 0.991/0.950, 0.990/0.950, 0.991/0.959; `engine3_explore50` 0.989/0.953, 0.990/0.947, 0.991/0.955, 0.992/0.961, 0.990/0.956.

Time per repetition: baselines 41–78 s; engines 480–550 s (17 refinement rounds, five network fits per round). Ground truth: 3 s.

## Findings

1. **About 8× fewer simulations for the same map.** `engine3` at 32k matches `sobol_1seed` at 256k on both metrics. At 16k it already beats every baseline at 64k.
2. **The gap is largest at medium budgets.** On the hard set the engine leads `sobol_1seed` by +7.6 points at 16k, +7.6 at 32k and +6.8 at 64k, far above the spread across repetitions.
3. **At very large budgets the methods converge.** At 256k the engine is still ahead on the hard set (0.957 vs 0.948) but the overall accuracies tie (0.991 vs 0.990). With unlimited simulations, careful ordering matters little; with expensive simulations it matters a lot.
4. **The engine plateaus after 64k and slips slightly at 256k** (hard 0.962 → 0.957). The late rounds keep packing points into a thin band near the boundary. Annealing or capping refinement late in the run is the obvious fix.
5. **More exploration does not help.** The 50/50 variant is never clearly better than the 30/70 split.
6. **Correction to stage 1b.** The fair map-maker lifts every method: `random_1seed` at 256k goes from 0.944 / 0.702 to 0.989 / 0.944. Stage 1b's "baselines never reach the engine" was an artefact of underfitting; the correct statement is "the engine reaches the same accuracy with far fewer runs".
7. **Cost of the engine itself.** Choosing points costs about 8 minutes of network training per repetition here, versus about 1 minute for the baselines. For the real target, where one simulation of a spatial protocell may take minutes or hours, this overhead is negligible next to an 8× saving in simulations.

## Relation to existing methods

The engine's parts are established techniques: quasi-random Sobol designs; query-by-committee (Seung, Opper & Sompolinsky, 1992) and deep ensembles (Lakshminarayanan et al., 2017) for disagreement-based sampling; separating model uncertainty from intrinsic randomness (e.g. BALD, Houlsby et al., 2011); and active learning of thresholds and level sets (Bryan et al., 2005; Gotovos et al., 2013). Robotic, algorithm-driven exploration of chemical space is also used in origin-of-life chemistry (e.g. the Cronin group). This project does not claim a new algorithm. What it adds is applying this combination to map the conditions under which a stochastic origin-of-life model produces life, with a benchmark against known ground truth, repetitions and negative results reported. A formal literature search is still to be done.

## Next steps

- Anneal or cap refinement in late rounds to remove the plateau and dip after 64k.
- Apply the engine to established origin-of-life models (autocatalytic sets, hypercycles, compositional/GARD-type models) to test that the saving holds beyond our own toy model.
- Stage 2: a spatial protocell model (membranes, growth, division), where each run is far more expensive and the 8× saving matters most.

## Raw data

`results_1c.json` and `bench_1c.log` are on the benchmark machine and will be added to this folder.
