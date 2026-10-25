# Stage 1d benchmark: faster engine and late-round selection (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Wall time:** 44 min · **Code:** [`stage1b/`](../stage1b) at commit `f2464a4`, run with `--fast`

## Bottom line

1. **The engine is 4× cheaper to run.** It now takes about 2 minutes per repetition instead of 8.5, with unchanged accuracy.
2. **Soft selection is the new default.** It draws refinement points with probability proportional to disagreement squared, and it is the best variant: better at small budgets and the most stable at the largest.
3. **The plateau after 64k runs is not fixed.** Every engine variant tops out near 0.96 on the hard set. The limit now seems to lie in the map-maker or the reference answers, not in point selection.

## Speed-up (`--fast`)

Profiling on the B70 showed that engine time went almost entirely into training small networks, not into simulation: 256,000 runs of the model take under a second. Each training step cost 2.3–3.6 ms whatever the data size. That means the GPU was waiting on kernel launches from the CPU rather than computing. Larger batches cost little extra (16× the data for about 2× the time). Bfloat16 and `torch.compile` gave no gain.

The fix trains the engine's five bootstrap networks as one batched network (`surrogate.fit_es_multi`, stacked weights and batched matrix products) and scores all checkpoints of a run with one batched fit. Each member keeps its own early stopping and best weights. Results:

| | before | `--fast` | speed-up |
|---|---|---|---|
| engine, seconds per repetition | ~510 | 124–133 | ~4× |
| baseline (`sobol_1seed`), seconds per repetition | ~60 | 14–20 | ~3.5× |

A parity check (2 repetitions, `sobol_1seed` and `engine3`) matched stage 1c within the spread between repetitions, for example engine3 at 32k: 0.989 / 0.939 before vs 0.989 / 0.940 with `--fast` (all / hard). The fast path uses different random draws, so numbers are not bit-identical.

The engine is still slower in wall time than the baselines on this toy model, because it retrains its networks in each of 17 rounds. What it saves is simulations. When one simulation is expensive, as it will be from stage 2 on, that is the cost that matters. For example, at 1 s per simulation, 256k baseline runs take about 71 hours and 32k engine runs about 9, plus about 2 minutes of network training.

## Variants tested

In stage 1c `engine3` stopped improving after 64k runs. Two fixes were tested, alone and together:

- `engine4_anneal`: the explore share rises linearly from 30% at 32k runs to 70% at the full budget.
- `engine4_soft`: refinement points are drawn with probability proportional to disagreement squared instead of strictly taking the top of the list, which spreads them along the boundary.
- `engine4_both`: both.

Same model, test set (8,000 conditions × 128 runs, 666 hard) and 5 repetitions as stage 1c. Command: `python bench1b.py --device xpu --surrogate es --fast --budget 256000 --reps 5 --n-test 8000 --gt-seeds 128 --methods sobol_1seed,engine3,engine4_anneal,engine4_soft,engine4_both --out results_1d.json`.

## Results

Accuracy on all 8,000 test conditions (mean ± standard deviation over 5 repetitions):

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k | 256k |
|---|---|---|---|---|---|---|---|---|
| sobol_1seed | 0.956 ± 0.005 | 0.964 ± 0.004 | 0.971 ± 0.000 | 0.975 ± 0.002 | 0.978 ± 0.002 | 0.984 ± 0.001 | 0.986 ± 0.001 | 0.989 ± 0.001 |
| engine3 | 0.956 ± 0.005 | 0.967 ± 0.003 | 0.977 ± 0.002 | 0.985 ± 0.001 | **0.989 ± 0.001** | **0.991 ± 0.001** | **0.991 ± 0.001** | 0.990 ± 0.001 |
| engine4_anneal | 0.956 ± 0.005 | 0.967 ± 0.003 | 0.977 ± 0.002 | 0.985 ± 0.001 | **0.989 ± 0.001** | 0.990 ± 0.001 | **0.991 ± 0.001** | **0.991 ± 0.001** |
| **engine4_soft** | 0.956 ± 0.005 | **0.970 ± 0.001** | **0.980 ± 0.001** | **0.986 ± 0.001** | 0.988 ± 0.001 | **0.991 ± 0.001** | **0.991 ± 0.001** | **0.991 ± 0.001** |
| engine4_both | 0.956 ± 0.005 | **0.970 ± 0.001** | **0.980 ± 0.001** | **0.986 ± 0.001** | 0.988 ± 0.001 | 0.990 ± 0.001 | **0.991 ± 0.001** | **0.991 ± 0.001** |

Accuracy on the 666 hard (stochastic) conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k | 128k | 256k |
|---|---|---|---|---|---|---|---|---|
| sobol_1seed | 0.772 ± 0.027 | 0.807 ± 0.021 | 0.826 ± 0.010 | 0.843 ± 0.012 | 0.866 ± 0.014 | 0.903 ± 0.011 | 0.920 ± 0.010 | 0.942 ± 0.009 |
| engine3 | 0.772 ± 0.027 | 0.811 ± 0.018 | 0.849 ± 0.014 | 0.916 ± 0.010 | **0.942 ± 0.011** | **0.959 ± 0.004** | 0.961 ± 0.006 | 0.953 ± 0.004 |
| engine4_anneal | 0.772 ± 0.027 | 0.811 ± 0.018 | 0.849 ± 0.014 | 0.916 ± 0.010 | **0.942 ± 0.011** | 0.957 ± 0.004 | 0.960 ± 0.006 | 0.959 ± 0.011 |
| **engine4_soft** | 0.772 ± 0.027 | **0.830 ± 0.020** | **0.874 ± 0.012** | **0.918 ± 0.006** | 0.938 ± 0.008 | 0.956 ± 0.009 | 0.962 ± 0.006 | **0.960 ± 0.002** |
| engine4_both | 0.772 ± 0.027 | **0.830 ± 0.020** | **0.874 ± 0.012** | **0.918 ± 0.006** | 0.938 ± 0.008 | 0.953 ± 0.005 | **0.963 ± 0.006** | 0.958 ± 0.005 |

(Variants share their early points by construction: anneal equals engine3 up to 32k, and both equals soft up to 32k.)

At 256k, per repetition (all / hard): `engine4_soft` 0.992/0.957, 0.991/0.960, 0.990/0.957, 0.991/0.962, 0.990/0.959; `engine4_anneal` 0.992/0.972, 0.991/0.962, 0.988/0.941, 0.991/0.963, 0.990/0.953; `engine3` 0.989/0.945, 0.990/0.953, 0.990/0.951, 0.991/0.953, 0.990/0.957.

## Findings

1. **Soft selection helps most when data are scarce.** On the hard set it scores +1.9 points over engine3 at 4k and +2.5 at 8k. At 256k it has the highest mean among the single fixes and by far the smallest spread (0.960 ± 0.002). It becomes the default engine.
2. **Annealing has weak evidence.** It raises the 256k hard-set mean from 0.953 to 0.959, but its spread across repetitions (±0.011) is as large as the gain.
3. **The ceiling near 0.96 is not a selection problem.** All engines reach about 0.96 on the hard set by 64k–128k and stop there. Likely causes are the map-maker's capacity and early-stopping noise, plus some noise in the 128-run reference answers. They will be examined when it matters for a real model.
4. **The stage 1c conclusion stands under the fast map-maker:** engines at 32k (0.988–0.989 / 0.938–0.942) match `sobol_1seed` at 256k (0.989 / 0.942), about 8× fewer simulations.

## Next

Test the engine on an established origin-of-life model: Kauffman's autocatalytic sets in the binary polymer formulation of Hordijk and Steel. It has a known threshold, about 1–2 reactions catalysed per molecule, against which the engine's map can be checked.

## Raw data

`results_1d.json`, `bench_1d.log` and the interrupted slow run (`bench_1d_partial_slow.log`) are on the benchmark machine and will be added to this folder.
