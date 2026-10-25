# Combination list over the system properties of the spatial model (stage 3 side result)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 · **Run:** about 9 min · **Summary data:** [stage3/combos_sp_summary.json](../stage3/combos_sp_summary.json)

**Read this first.** The ten conditions here are the conditions of the [stage 3d](stage3d_spatial_b70.md) spatial model. Most of them (division size S, parasite advantage σ, gene types K, error rate μ, membrane synthesis km, heritable variation η) describe the chemical system itself, not the environment. Only food, degradation and leak correspond roughly to the project's original environmental conditions. So this list answers "which chemical systems of this model become alive", not "which environments give rise to life". The environment-driven version is planned for stage 4 ([sketch](../docs/stage4_sketch.md)).

## Method

- **Grid.** Each of the ten conditions takes 4 representative values (at u = 0.125, 0.375, 0.625, 0.875 of its range). K takes 1, 2 or 3. The distractors are fixed at 0.5. That gives 3 × 4⁹ = 786,432 combinations.
- **Maps.** Each combination is scored with the COMPASS maps of the spatial model (32,000 simulations per level, refitted with the same seeds as in 3d). For each combination we record P(L1)…P(L4) and the highest level with P > 0.5.
- **Check.** 2,000 combinations, 400 per predicted highest level, were simulated directly with 16 seeds each.

## Results

| highest level | combinations | share |
|---|---|---|
| none | 364,805 | 46.4% |
| L1 self-maintenance | 245,154 | 31.2% |
| L2 heredity | 94,752 | 12.0% |
| L3 coupled protocell | 32,746 | 4.2% |
| L4 evolution | 48,975 | 6.2% |

**Accuracy against direct simulation.** Agreement on the yes/no label per level is 97.8% (L1), 97.6% (L2), 93.8% (L3) and 93.1% (L4). The highest level agrees in 85% of cases overall. For combinations predicted to stop at exactly L3 it agrees in only 63%, because they sit between the L2 and L4 boundaries.

**Nesting.** The levels are mapped separately, and 7,466 combinations (0.9%) break the nesting, for example L4 above 0.5 while L3 is just below. This is why more combinations end at L4 than at L3. The joint nested map of stage 3e is meant to remove this.

**Typical evolving systems.** In the 20 most robust L4 patterns, all or almost all of the 256 settings of the remaining conditions evolve. These patterns share:
- moderate membrane synthesis (km ≈ 0.85);
- small cells (S ≈ 13–38);
- one gene type;
- low error (μ ≤ 0.0085);
- weak parasites (σ ≤ 1.75);
- plenty of food (0.84).

Degradation, membrane turnover, leak and heritable variation hardly matter within these patterns.

**Single conditions.** No single grid value rules out evolution completely, but the share of L4 combinations varies strongly. Examples:

| condition | value | share of combinations reaching L4 |
|---|---|---|
| membrane synthesis km | 0.2 | 0.15% |
| | 3.5 | 16% |
| | 14.7 | 0.36% |
| division size S | 13 | 11.6% |
| | 304 | 1.1% |
| heritable variation η | 0.025 | 0.7% |
| | 0.175 | 9.4% |

## Files

The full table (`combos_sp.csv.gz`, 3.4 MB, one row per combination) and the reduced table with only L3/L4 rows (`combos_sp_L3L4.csv.gz`, 0.5 MB, 81,721 rows) are on the benchmark machine. They will be uploaded with the other raw results.
