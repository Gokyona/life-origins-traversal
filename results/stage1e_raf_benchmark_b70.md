# Stage 1e benchmark: the engine on an established model, Kauffman autocatalytic sets (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Run:** 13:09–13:53 Toronto time (44 min) · **Code:** [`stage1e_raf/`](../stage1e_raf) at commit `38fab7d` with the stage 1b harness at `f2464a4`, run with `--fast`

## Bottom line

**The engine's advantage carries over to a model we did not design, but it is smaller: about 2× fewer simulations instead of 8×.** On the RAF model, the default engine (`engine4_soft`) reaches with 16k runs the accuracy that Sobol sampling reaches with 32k (hard 0.937 vs 0.939, all 0.966 vs 0.966), and with 32k what Sobol reaches with 64k (hard 0.955 vs 0.955, all 0.970 vs 0.970). How much the engine saves depends on the shape of the problem.

## The model

Kauffman's autocatalytic sets in the binary polymer formulation of Hordijk & Steel (2004). Molecules are binary strings up to length n, each split of a string is one reversible ligation/cleavage reaction, the food set is every string up to length t, and each molecule catalyses each reaction with probability p. "Life" here means a RAF set exists: a self-sustaining autocatalytic network, detected exactly with the Hordijk–Steel maxRAF reduction batched on the GPU. That corresponds to level L1 of our life detector.

Seven conditions: maximum length n (5–10), catalysis level f = mean reactions catalysed per molecule (0.1–10, log scale), food length t (1–3), fraction of possible reactions present (0.3–1), fraction of molecules able to catalyse (0.1–1), and two distractors.

Before the benchmark, the implementation was checked against theory. With t = 2 and all reactions present, P(RAF) = 0.5 at f ≈ 1.1 (n = 6), 1.25 (n = 8) and 1.3 (n = 10), and the transition sharpens with n. This matches the published threshold of 1–2 reactions per molecule. The design check and its caveats (the catalyst-fraction condition turned out to have no measurable effect, and reaction density partly overlaps with f) are in [`stage1e_raf/README.md`](../stage1e_raf/README.md).

One run costs about 1 ms on the B70, roughly 1,000× more than the stage 1b toy model, so the budget is 64k runs. Ground truth: 4,000 random conditions × 64 runs (159 s). Life fraction 0.418; hard set (true success probability 0.1–0.4 or 0.6–0.9) 932 conditions, 23% of the space versus 8% in stage 1b.

Command: `python bench_raf.py --device xpu --surrogate es --fast --budget 64000 --reps 5 --n-test 4000 --gt-seeds 64 --methods random_1seed,sobol_1seed,engine3,engine4_soft --out results_raf.json`.

## Results

Accuracy on all 4,000 test conditions (mean ± standard deviation over 5 repetitions):

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| random_1seed | 0.934 ± 0.009 | 0.941 ± 0.004 | 0.948 ± 0.002 | 0.958 ± 0.006 | 0.966 ± 0.002 | 0.970 ± 0.003 |
| sobol_1seed | 0.940 ± 0.004 | 0.945 ± 0.003 | 0.952 ± 0.004 | 0.957 ± 0.003 | 0.966 ± 0.002 | 0.970 ± 0.003 |
| engine3 | 0.942 ± 0.004 | 0.946 ± 0.005 | 0.956 ± 0.002 | **0.966 ± 0.001** | **0.970 ± 0.001** | **0.973 ± 0.003** |
| **engine4_soft** | 0.942 ± 0.004 | **0.949 ± 0.003** | **0.959 ± 0.003** | **0.966 ± 0.002** | **0.970 ± 0.002** | **0.973 ± 0.001** |

Accuracy on the 932 hard (stochastic) conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| random_1seed | 0.837 ± 0.030 | 0.859 ± 0.020 | 0.880 ± 0.008 | 0.914 ± 0.019 | 0.939 ± 0.007 | 0.952 ± 0.008 |
| sobol_1seed | 0.856 ± 0.011 | 0.869 ± 0.008 | 0.891 ± 0.012 | 0.912 ± 0.009 | 0.939 ± 0.009 | 0.955 ± 0.007 |
| engine3 | 0.856 ± 0.011 | 0.875 ± 0.018 | 0.907 ± 0.005 | **0.938 ± 0.002** | 0.952 ± 0.006 | 0.963 ± 0.007 |
| **engine4_soft** | 0.856 ± 0.011 | **0.876 ± 0.012** | **0.914 ± 0.010** | 0.937 ± 0.005 | **0.955 ± 0.006** | **0.964 ± 0.002** |

At 64k per repetition (all / hard): `engine4_soft` 0.972/0.960, 0.974/0.965, 0.972/0.963, 0.973/0.966, 0.973/0.963; `sobol_1seed` 0.968/0.950, 0.970/0.960, 0.966/0.945, 0.974/0.964, 0.971/0.954.

Time per repetition: baselines 74–96 s (almost all of it simulation and scoring); engines 162–184 s (12 refinement rounds).

## Findings

1. **The advantage transfers, at about 2×.** From 8k runs on, both engines beat both baselines on the hard set, and `engine4_soft` matches Sobol at twice the budget. On the toy model the saving was about 8×.
2. **Why the saving is smaller here.** The uncertain region is broad: 23% of the space is hard, and the life boundary is spread along f, t and reaction density rather than forming a thin shell. Concentrating tests near the boundary pays off most when the boundary is thin. Part of the ceiling near 0.97 overall also comes from the 64-run reference answers.
3. **Soft selection is again the best and most stable engine.** It is ahead of engine3 at 4k–8k (hard 0.914 vs 0.907 at 8k) and has the smallest spread at 64k (0.964 ± 0.002).
4. **The engine's own cost is now about the same size as the simulation cost.** Training adds about 90–100 s per repetition, against about 75 s for 64k simulations at 1 ms each. As simulations get more expensive (stage 2), the training overhead becomes relatively small and the saving in simulations dominates.

## What this means for the project

The ordered engine is not tied to our own toy model: on a published origin-of-life model with a known threshold, it maps the life boundary with fewer simulations than random or Sobol sampling. The size of the gain is problem-dependent, between about 2× and 8× in the two models tested so far, so it should be reported as a range, not a single number.

## Raw data

`results_raf.json` and `bench_raf.log` are on the benchmark machine and will be added to this folder.
