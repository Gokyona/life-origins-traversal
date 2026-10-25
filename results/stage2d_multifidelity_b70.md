# Stage 2d: multi-fidelity COMPASS on coupled protocells (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Runs:** fidelity check; main benchmark 15:47–16:19 Toronto time; follow-up 16:26–17:07 · **Code:** [`stage2/fidelity_check.py`](../stage2/fidelity_check.py), [`stage2/mf_bench.py`](../stage2/mf_bench.py), `protocell.simulate(..., fidelity=0|1)`

## Bottom line

**A partial negative result.** Screening with cheap, short simulations helps a little when the budget is small, but at larger budgets no multi-fidelity variant beats COMPASS that uses full simulations only. The best variant (`engine5_decay`) is slightly ahead up to 8k and ties from 32k on. COMPASS on full simulations (`engine4_soft`) stays the default engine; multi-fidelity is kept as an option for small budgets or for models where the cheap version is much cheaper than here.

## The cheap simulator

Fidelity 0 runs 32 cell slots for 300 steps instead of 128 for 1,000; the generation requirement of L2/L3 is scaled with the number of steps. On 4,000 random conditions × 16 runs at each fidelity (level L3):

| measure | value |
|---|---|
| cost per simulation | 0.078 ms (cheap) vs 0.611 ms (full): **7.8× cheaper** |
| P(L3) over the space | 0.164 (cheap) vs 0.148 (full) |
| label agreement (p > 0.5) | **0.972**; confusion (rows full, columns cheap): no/no 3319, no/yes 89, yes/no 25, yes/yes 567 |
| correlation of probabilities | 0.934 |
| agreement on clear-cut conditions (full p ≤ 0.1 or ≥ 0.9; 97% of the space) | 0.985 |
| agreement near the boundary | 0.597 |
| bias near the boundary | where full p is in [0.1, 0.4), cheap p averages 0.45 vs 0.23: **the cheap model is optimistic exactly at the boundary** |

## Methods

Budgets are counted in full-simulation-equivalent cost (one cheap run = 1/7.8). All methods are scored by the same map-maker, predicted at full fidelity, against full-fidelity ground truth (4,000 conditions × 64 runs; life fraction 0.162; hard set 107 conditions, 2.7%).

- `sobol_1seed`: full simulations at Sobol points.
- `engine4_soft`: COMPASS on full simulations only (default since stage 1d).
- `engine5_mf`: every round spends 30% of its cost on cheap Sobol screening and the rest on full simulations chosen by COMPASS; the map-maker gets the fidelity as an extra input.
- `engine5_decay`: cheap share 50% in the first round, falling linearly to 0 at half the budget.
- `engine5_screen`: cheap runs only rule regions out: cheap "life" results are discarded, and exploration points the map calls clearly lifeless (mean p < 0.02 and ensemble spread < 0.02) are skipped. Discarded runs still count towards the cost.
- `engine5_both`: decay and screen together.

## Results

Accuracy on all 4,000 test conditions (mean ± standard deviation over 5 repetitions), by full-simulation-equivalent cost:

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| sobol_1seed | 0.945 ± 0.003 | 0.955 ± 0.002 | 0.957 ± 0.002 | 0.962 ± 0.002 | 0.969 ± 0.002 | 0.977 ± 0.002 |
| **engine4_soft** | 0.945 ± 0.004 | 0.959 ± 0.001 | 0.964 ± 0.001 | **0.977 ± 0.002** | **0.983 ± 0.002** | **0.986 ± 0.001** |
| engine5_mf | 0.951 ± 0.001 | 0.957 ± 0.004 | **0.969 ± 0.002** | 0.975 ± 0.002 | 0.982 ± 0.002 | 0.984 ± 0.002 |
| engine5_decay | 0.951 ± 0.002 | **0.961 ± 0.002** | 0.967 ± 0.002 | 0.974 ± 0.001 | 0.981 ± 0.001 | 0.985 ± 0.001 |
| engine5_screen | 0.935 ± 0.008 | 0.955 ± 0.002 | 0.962 ± 0.002 | 0.976 ± 0.001 | 0.980 ± 0.002 | 0.983 ± 0.002 |
| engine5_both | 0.924 ± 0.009 | 0.954 ± 0.003 | 0.960 ± 0.002 | 0.971 ± 0.003 | 0.979 ± 0.001 | 0.985 ± 0.001 |

Accuracy on the 107 hard conditions:

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| sobol_1seed | 0.581 ± 0.049 | 0.611 ± 0.056 | 0.619 ± 0.035 | 0.637 ± 0.024 | 0.641 ± 0.022 | 0.705 ± 0.024 |
| engine4_soft | 0.578 ± 0.026 | 0.641 ± 0.017 | 0.630 ± 0.024 | 0.693 ± 0.022 | 0.736 ± 0.049 | 0.774 ± 0.028 |
| engine5_mf | 0.585 ± 0.025 | 0.611 ± 0.046 | 0.662 ± 0.022 | 0.665 ± 0.040 | 0.742 ± 0.032 | 0.755 ± 0.035 |
| engine5_decay | 0.594 ± 0.025 | 0.624 ± 0.019 | 0.637 ± 0.029 | 0.621 ± 0.037 | 0.723 ± 0.013 | 0.774 ± 0.037 |
| engine5_screen | 0.594 ± 0.030 | 0.621 ± 0.017 | 0.647 ± 0.027 | 0.673 ± 0.042 | 0.733 ± 0.041 | 0.750 ± 0.030 |
| engine5_both | 0.596 ± 0.036 | 0.604 ± 0.040 | 0.641 ± 0.025 | 0.650 ± 0.016 | 0.725 ± 0.030 | 0.776 ± 0.031 |

Runs used at each checkpoint (mean over repetitions, cheap / full; for the screen variants the cheap count is the runs kept, all are counted in cost):

| method | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| engine5_mf | 4792/1385 | 9584/2771 | 22674/5093 | 40366/10824 | 82212/21460 | 149774/44798 |
| engine5_decay | 7987/976 | 15335/2033 | 31488/3963 | 51983/9335 | 72385/22719 | 72385/54719 |
| engine5_screen | 3937/1385 | 7877/2771 | 18676/5093 | 33247/10824 | 67675/21460 | 123303/44798 |
| engine5_both | 6565/976 | 12615/2033 | 25923/3963 | 42778/9335 | 59567/22719 | 59567/54719 |

The screen variants skipped about 33k (`screen`) and 43k (`both`) clearly lifeless exploration points in total over 5 repetitions. Time per repetition: Sobol 55–58 s, engine4_soft 139–145 s, multi-fidelity variants 155–165 s.

## Findings

1. **Cheap screening helps only while full simulations are scarce.** At 2k–8k the multi-fidelity variants are up to 0.006 ahead of engine4_soft (8k: 0.969 vs 0.964, about what engine4_soft reaches at 12k).
2. **At larger budgets it does not pay here.** From 16k on every variant is level with or slightly behind engine4_soft (64k: 0.983–0.985 vs 0.986); all stay clearly ahead of Sobol.
3. **Why.** The cheap model is only 7.8× cheaper and is biased exactly where accuracy is decided: near the boundary it is optimistic. Money spent on cheap runs is money not spent on full runs at the boundary.
4. **Decay is the safest variant**: never meaningfully worse, slightly better early. **Screen alone is worse early**: once cheap "life" results are thrown away, the cheap runs carry little information, and skipping exploration points thins the early map.

## Decision

Keep `engine4_soft` (COMPASS on full simulations) as the default. Use `engine5_decay` when the budget is small or when a cheap model is far cheaper (tens to hundreds of times) than the full one, as expected for spatial or more detailed models later. A correction map that learns the cheap model's bias near the boundary is a possible future improvement.

## Next

Stage 2e: maps of the conditions for life levels L1, L2 and L3 with COMPASS, and a comparison of which conditions each level needs.

## Raw data

`results_mf.json`, `bench_mf.log`, `results_mf2.json` and `bench_mf2.log` are on the benchmark machine and will be added to this folder.
