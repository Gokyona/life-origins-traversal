# Stage 2: dividing protocells

Plan: [`../docs/stage2_plan.md`](../docs/stage2_plan.md).

## 2a: stochastic corrector population (`scm.py`)

A well-mixed population of 128 protocells after Szathmáry & Demeter (1987). Each cell holds K mutually dependent gene types and a parasite; balanced, parasite-poor cells grow faster; cells divide at size S and share their molecules at random; copying errors turn genes into parasites. Eight conditions: division size S (8–512), parasite advantage σ (1–3), mutation rate μ (0.001–0.3), number of gene types K (1–3), degradation, food, and two distractors. Detector: L1 (some cell still functional) and L2 (at least half the cells functional after ≥ 15 divisions per cell slot). Deviations from the original model are listed in the module docstring.

Check against known behaviour (200 populations per point, P(L2)):

| sweep | result |
|---|---|
| division size S, σ = 2, μ = 0.01, K = 2 | 1.00 for S ≤ 32, 0.98 at 64, 0.11 at 128, 0 at ≥ 256: parasites win in large cells |
| same, σ = 1.2 | 1.00 up to S = 256, 0.25 at 512: weaker parasites allow larger cells |
| K = 3, σ = 1, μ = 0.001 | 0 at S = 8 (genes lost at division), 1.00 for S ≥ 16; K = 1 and 2 always 1.00 |
| mutation rate μ, S = 64, σ = 1.5 | 1.00 for μ ≤ 0.01, 0.56 at 0.03, 0 at ≥ 0.1: a critical mutation rate |
| food, S = 64, σ = 1.5, μ = 0.01 | 0 at food ≤ 0.35, 0.14 at 0.5, 1.00 at ≥ 0.7 |

About 0.8 ms per population on the Intel Arc Pro B70. Over the whole space P(L1) = 0.615 and P(L2) = 0.198. Because a population of 128 cells averages out much of the randomness, only about 2% of conditions are "hard" (success probability 0.1–0.4 or 0.6–0.9): the boundary is sharp.

```bash
python validate_scm.py --device xpu
python bench_scm.py --device xpu --surrogate es --fast --budget 64000 --reps 5 --n-test 4000 --gt-seeds 64 \
  --methods random_1seed,sobol_1seed,engine4_soft --out results_scm.json
```

## 2b: division driven by the cell's own chemistry (`protocell.py`)

Extends 2a with a membrane, after Gánti's chemoton idea that replication, metabolism and membrane growth are coupled. Gene type 0 is the membrane enzyme: it makes lipid from food. The membrane sets how much a cell can hold; crowding slows copying, and a cell holding more than twice its capacity bursts. **A cell divides only when its membrane has doubled**; there is no rule on its contents. Ten conditions: the eight of 2a (S is now the membrane division size) plus membrane synthesis rate km (0.1–30, log) and membrane turnover dm (0–0.02).

New detector level **L3 (coupled protocell)**: L2 holds, and the cells are in homeostasis: the mean content/capacity ratio of functional cells at the end is between 0.25 and 1.5, and bursts in the second half of the run are fewer than 10% of divisions.

Check against expected behaviour (200 populations per point; reference S = 64, σ = 1.5, μ = 0.01, K = 2, d = 0.005, food = 1, dm = 0.005):

| sweep | result |
|---|---|
| membrane rate km at the reference | 0.1: L2 = L3 = 0, 4.9 generations, bursts = 29% of divisions; 0.15–0.2: 0 (8–11 generations); **0.3–5: L3 = 1.00** (17–66 generations, ratio 1.30 → 0.33); 10 and 30: L2 = 1.00 but L3 = 0 (ratio 0.21, 0.10: under-filled vesicles) |
| same, food = 0.4 | L3 window narrows to km ≈ 2–5 |
| division size S, σ = 2 | L3 = 1.00 for S ≤ 128, 0.29 at 256, 0 at 512: parasites still win in large cells |
| K = 3, σ = 1, μ = 0.001 | 0 at S = 8, 1.00 for S ≥ 16 (genes lost at division in tiny cells); K = 1, 2 always 1.00 |
| mutation rate μ | L3 = 1.00 up to μ = 0.03, 0 at ≥ 0.1: critical mutation rate kept |

So division emerges from the internal chemistry in a reference setting (the 2b acceptance criterion in the plan), with a window in membrane synthesis rate: too slow gives crowding and bursting, too fast gives vesicles that split before their contents double. Contrary to our first expectation, the fast-membrane failure is under-filling rather than outright gene loss (at S ≥ 16).

About 0.9 ms per population on the Intel Arc Pro B70. Over the whole space P(L1) = 0.561, P(L2) = 0.306, P(L3) = 0.152.

```bash
python validate_protocell.py --device xpu
python bench_protocell.py --device xpu --surrogate es --fast --budget 64000 --reps 5 --n-test 4000 --gt-seeds 64 \
  --methods random_1seed,sobol_1seed,engine4_soft --out results_protocell.json
```
