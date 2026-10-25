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
