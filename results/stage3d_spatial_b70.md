# Stage 3d: cross-check on an independent spatial model (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Design and validation:** from about 18:50 Toronto time · **Map run:** 18:51–19:02 Toronto time (11 min) · **Plan:** [stage3_plan.md](../docs/stage3_plan.md) · **Code:** [stage3/spatial.py](../stage3/spatial.py), [validate_spatial.py](../stage3/validate_spatial.py), [maps_sp.py](../stage3/maps_sp.py)

## Bottom line

We rebuilt the evolving protocell model of [stage 3c](stage3c_evolution_b70.md) with space: cells sit on a grid, a daughter can only take a neighbouring place, and molecules may leak to neighbours. Every qualitative conclusion of stages 2 and 3a–3c holds in this second model. Space makes the boundaries tighter: with almost no leak, the L2 and L3 regions are about a quarter smaller than in the well-mixed model, and the largest cell size that still allows evolution is roughly halved under strong parasites. Space also brings a new, strongly limiting condition: **compartments must be nearly sealed**. At the reference, heredity collapses once more than about 0.8% of molecules leak per step, and evolution once more than about 0.7% do. COMPASS again maps the space more accurately than Sobol sampling.

## Model

`spatial.py` is `evolve.py` (stage 3c: levels L1–L4, heritable quality q, eight neutral markers) with three changes:

- **Grid.** The 128 cell slots form a 16 × 8 torus (8 × 4 at low fidelity); each slot has 8 neighbours.
- **Local placement.** The second daughter goes to a neighbouring slot: a dead or empty one first (unless the `no_selection` knockout is on), never one that is itself dividing, otherwise a random neighbour. If two dividing cells pick the same slot in the same step, a random one wins and the other's second daughter is lost (its parent still divides). Selection is therefore local: a fast lineage spreads only into its neighbourhood.
- **Leak.** Each molecule moves to a random neighbouring cell with probability lk per step (all molecules leaving one cell in a step go to the same random neighbour). lk is a new condition, 0 to 0.01.

The condition space now has twelve dimensions: the ten real conditions S, σ, μ, K, d, food, km, dm, η, lk, and the two distractors. Cost: 1.28 ms per simulation, against 1.03 ms for the well-mixed model.

## Validation (200 populations per setting; spatial | well-mixed; lk = 0 unless stated)

P(L1 L2 L3 L4) and, where useful, mean q at the end. With 200 populations the binomial standard deviation is at most 0.035.

| setting | spatial | well-mixed |
|---|---|---|
| reference (η = 0.1) | 1 / 1 / 1 / 0.81, q 3.27 | 1 / 1 / 1 / 0.88, q 3.62 |
| σ = 2, S = 16 | L4 0.81 | L4 0.83 |
| σ = 2, S = 64 | L4 0.70 | L4 0.78 |
| σ = 2, S = 128 | L3 0.97, L4 0.38 | L3 1.00, L4 0.49 |
| σ = 2, S = 256 | L3 0.02, L4 0.01 | L3 0.20, L4 0.04 |
| σ = 2, S = 512 | 0 | 0 |
| μ = 0.003 / 0.01 / 0.03 | L4 0.92 / 0.82 / 0.44 | L4 0.92 / 0.82 / 0.53 |
| μ = 0.1 | L2 0 | L2 0 |
| km = 0.1 / 0.3 / 2 / 10 | L2 0; L3 0.98 & L4 0.37; L4 0.86; L3 0 | L2 0; L3 1.00 & L4 0.53; L4 0.82; L3 0 |
| K = 3, S = 8 (tiny cells) | 0 | 0 |

**Knockouts in space.** `equal_split` (no random sharing at division): L4 0.69 at the reference (well-mixed 0.75); in the corner S = 256, σ = 2, L3 falls to 0.00 (well-mixed 0.01; baselines 0.02 and 0.24). `no_selection` (no preference for dead slots): L4 0.81 at the reference in both; in the corner L3 0.04 (well-mixed 0.18). The 3a conclusions therefore hold in space: random sharing matters where cells are large and parasites strong, and dropping the dead-slot preference changes little at the reference but costs in that corner.

**Leak at the reference** (400 populations per value, binomial sd ≤ 0.025), P(L2) / P(L4):

| lk | 0–0.0045 | 0.005 | 0.0055 | 0.006 | 0.0065 | 0.007 | 0.0075 | 0.008 | 0.0085 | 0.009 | 0.0095 | 0.01 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L2 | 1.00 | 1.00 | 0.99 | 0.97 | 0.90 | 0.78 | 0.58 | 0.42 | 0.26 | 0.16 | 0.08 | 0.02 |
| L4 | 0.71–0.85 | 0.68 | 0.62 | 0.60 | 0.53 | 0.38 | 0.32 | 0.22 | 0.12 | 0.08 | 0.06 | 0.02 |

The critical leak (first grid point below 0.5) is 0.008 for heredity (L2) and 0.007 for evolution (L4). L1 stays at 1.00 throughout: leaky cells keep growing, but they no longer keep their genes together. In the corner S = 256, σ = 2, any leak of 0.001 or more removes life entirely. Leaked parasites infect neighbouring cells, which undoes the protection that compartments give.

**Whole space** (4,096 random conditions): spatial P(L1) 0.539, P(L2) 0.223, P(L3) 0.102, P(L4) 0.062; well-mixed 0.569, 0.345, 0.183, 0.114.

## Map (COMPASS vs Sobol, 32,000 simulations each, for L3 and L4; ground truth 4,000 points × 64 seeds for L1–L4)

| | L1 | L2 | L3 | L4 |
|---|---|---|---|---|
| spatial, all conditions | 0.527 | 0.223 | 0.102 | 0.071 |
| spatial, lk < 0.001 only (n = 401, ± about 0.02) | 0.531 | 0.259 | 0.135 | 0.100 |
| well-mixed ([stage 3c](stage3c_evolution_b70.md), no lk) | 0.558 | 0.338 | 0.173 | 0.112 |

Part of the smaller spatial region comes from leak, part from space itself: even with almost no leak, the L2 and L3 regions are about a quarter smaller than well-mixed, while L4 is only slightly smaller. The four levels are perfectly nested in the ground truth, and 70% of L3 conditions reach L4 (64% well-mixed).

**Accuracy** (overall / hard points near the boundary):

| level | sims | 2k | 4k | 8k | 16k | 32k |
|---|---|---|---|---|---|---|
| L3 (hard set 100) | COMPASS | 0.952 / 0.57 | 0.960 / 0.57 | 0.968 / 0.62 | 0.976 / 0.66 | **0.983 / 0.72** |
| | Sobol | 0.948 / 0.56 | 0.951 / 0.56 | 0.959 / 0.61 | 0.966 / 0.63 | 0.967 / 0.64 |
| L4 (hard set 248) | COMPASS | 0.956 / 0.621 | 0.967 / 0.669 | 0.971 / 0.730 | 0.975 / 0.754 | **0.977 / 0.786** |
| | Sobol | 0.954 / 0.605 | 0.960 / 0.593 | 0.965 / 0.702 | 0.970 / 0.698 | 0.969 / 0.681 |

At 32k COMPASS has about half the overall errors of Sobol for L3 (1.7% vs 3.3%) and 10.5 points higher hard accuracy for L4. The hard sets here are small (100 and 248 points), so the hard accuracies carry several points of sampling error.

**Importance** (drop in accuracy when a condition is shuffled):
- L3: km 0.115, S 0.083, K 0.065, μ 0.048, food 0.044, σ 0.041, d 0.019, lk 0.017, dm 0.012, η 0.002, distractors ≤ 0.
- L4: km 0.072, S 0.053, K 0.045, μ 0.031, σ 0.029, food 0.026, η 0.025, lk 0.009, d 0.005, dm 0.003, distractors 0.

lk ranks low in importance even though it is decisive at the reference. Most of its range (lk < 0.006) is harmless at small cell sizes, and it bites mainly together with large S (see below).

**Largest division size S that still gives life, at lk = 0** (map read-off; "–" = no limit within S ≤ 512):

| σ | 1.0 | 1.2 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|---|
| L3 spatial | – | – | 338 | 239 | 194 | 169 |
| L3 well-mixed ([2e](stage2e_maps_b70.md)) | – | – | 388 | 256 | 223 | 223 |
| L4 spatial | – | 416 | 194 | 111 | 84 | 74 |
| L4 well-mixed ([3c](stage3c_evolution_b70.md)) | – | – | 256 | 181 | 147 | 169 |

**Leak and cell size interact.** In the lk × S slice at the reference, the largest cell size for L3 falls from about 340 without leak to about 55 at lk = 0.01. The more leaky the compartments, the smaller cells must be.

Figures: [accuracy](figures_3d/accuracy.png), [importance](figures_3d/importance.png), partial dependence [L3](figures_3d/pd_L3.png) and [L4](figures_3d/pd_L4.png), slices [lk × S](figures_3d/slice_lk_x_S.png), [lk × σ](figures_3d/slice_lk_x_sigma.png), [η × food](figures_3d/slice_eta_x_food.png), [S × σ](figures_3d/slice_S_x_sigma.png).

## What this means

1. **The stage 3 conclusions are not an artefact of one model design.** Compartments, small cells, a low mutation rate, a window of membrane synthesis, and heritable variation for evolution are all required in the spatial model as well.
2. **Compartments must be nearly sealed.** Leak above about 0.7% of molecules per step destroys heredity and evolution at the reference, and much less leak is tolerated by large cells.
3. **Local competition narrows the room for evolution.** Under strong parasites the largest cell size that still evolves is roughly halved. A possible reason is that good lineages spread more slowly when they can only take neighbouring places, while parasites spread within neighbourhoods. This is a hypothesis we have not tested.

## Limits

- One grid size (16 × 8) and one neighbourhood (8 neighbours). Larger grids may allow spatial patterns (waves, clusters) that this small torus cannot show.
- The leak model sends all molecules leaving a cell in one step to the same neighbour. This is a simplification.
- As in 3c, L4 is measured over a fixed run length.

## Raw data

`maps_sp_compact.json` (in [stage3/](../stage3/)) holds curves, importance, partial dependence, slices, the lk sweep and the read-offs. Full logs and result files: [raw/stage3/](../raw/stage3/).
