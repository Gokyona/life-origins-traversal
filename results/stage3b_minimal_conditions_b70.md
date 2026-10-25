# Stage 3b: minimal sufficient conditions for each level of life (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Run:** 18:20–18:28 Toronto time (8 min) · **Plan:** [stage3_plan.md](../docs/stage3_plan.md)

## Bottom line

For each life level of the stage 2b protocell model we found a small set of conditions which, when held in the given ranges, give life with about 95% probability whatever the other conditions do. Each set was confirmed with fresh direct simulations. The higher the level, the more conditions must be held and the narrower the ranges: the boxes cover 17%, 6.9% and 1.9% of the condition space for L1, L2 and L3.

## Method

1. Start with all ten conditions free over their full ranges.
2. Using the COMPASS map of each level (32k simulations, as in [stage 2e](stage2e_maps_b70.md)), estimate the share of random points inside the current box that the map calls "life".
3. At each step, try every contiguous interval (on a 10-bin grid) on every condition, and constrain the one that buys the largest gain for the smallest loss of box volume.
4. Stop when at least 95% of the box is life (a 90% variant is also reported).
5. Confirm the final box by direct simulation: 2,000 random points inside it × 16 simulations each.

## Results (95% target)

| level | conditions that must be held | free conditions | box volume | direct check: share of points with P > 0.5 / mean P |
|---|---|---|---|---|
| L1 self-maintenance | mutation rate μ ≤ 0.011; membrane synthesis km ≥ 0.67; division size S ≥ 32; degradation d ≤ 0.018 | σ, K, food, dm, z1, z2 | 17.0% | 0.960 / 0.958 |
| L2 heredity | km ≥ 2.8; μ ≤ 0.017; food ≥ 0.40; S ≥ 23; gene types K ≤ 2 | σ, d, dm, z1, z2 | 6.9% | 0.946 / 0.945 |
| L3 coupled protocell | **km 1.08–2.79**; μ ≤ 0.017; food ≥ 0.50; **S 11–181**; K ≤ 2 | σ, d, dm, z1, z2 | 1.9% | 0.950 / 0.949 |

With a 90% target the boxes grow to 25% (L1), 10.8% (L2) and 2.1% (L3) of the space; direct checks give 0.914, 0.914 and 0.933.

## Findings

1. **Each level adds requirements.** L1 needs low mutation, enough membrane, cells that are not tiny, and moderate degradation. L2 adds enough food (to complete many generations) and at most two gene types. L3 narrows membrane synthesis to a window and division size to a range with an upper limit.
2. **Parasite advantage never needs to be held.** Every box also holds the mutation rate low (μ ≤ 0.011 for L1, ≤ 0.017 for L2 and L3). Parasites arise only by mutation, so they stay rare and their replication advantage matters little. The direct confirmation points show this: in the corner with large cells and strong parasites (S > 256, σ > 2), life is still found at 91% of points in the L1 box (257 points) and 91% in the L2 box (217 points), only about 4 points below the box averages; the L3 box excludes that corner through its upper limit on cell size. Degradation and membrane turnover can also be left free for L2 and L3.
3. **The distractors are never chosen**, as expected.
4. **Sufficient, not unique.** These boxes are one simple description of a region that guarantees life. Life also exists outside them, in regions that a single box cannot describe.

## Raw data

`s3b.json` and `s3b.log` (full greedy step sequences) are on the benchmark machine and will be added with the search code.
