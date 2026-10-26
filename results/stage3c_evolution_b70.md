# Stage 3c: Darwinian evolution (L4) in the protocell model (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Design and validation:** from about 18:30 Toronto time · **Map run:** 18:36–18:42 Toronto time (6 min) · **Plan:** [stage3_plan.md](../docs/stage3_plan.md) · **Code:** [stage3/evolve.py](../stage3/evolve.py), [validate_evolve.py](../stage3/validate_evolve.py), [maps_l4.py](../stage3/maps_l4.py)

## Bottom line

We added a fourth level of life, **L4: Darwinian evolution**, to the stage 2b protocell model. A population reaches L4 when it is a working protocell population (L3) *and* a heritable quality of its machinery rises clearly above what genetic drift alone produces. L4 occupies 11% of the condition space. Every L4 point is an L3 point, and 64% of L3 conditions also reach L4. Evolution needs three things beyond L3: some heritable variation (η ≥ about 0.02 per division), a mutation rate for parasites below about 0.04, and enough food to run many generations. Less food can be offset by more heritable variation. COMPASS again maps the level better than plain Sobol sampling with the same budget.

## Design

- **Heritable quality.** Every cell carries a number q (start 1, clipped to [0, 10]). It speeds up the cell's whole metabolism with diminishing returns, h(q) = 2q/(1+q): h(1) = 1 and h never exceeds 2, so quality can at most double the rate. h multiplies gene, parasite and membrane synthesis alike, so faster cells do not break homeostasis merely by growing faster.
- **Inheritance.** At division each daughter receives its parent's q plus a random step N(0, η). The size of the step, **η** (0 to 0.2), is a new condition, so the condition space now has eleven dimensions (nine real ones plus the two distractors).
- **Drift control inside every simulation.** Each cell also carries **8 neutral markers**, inherited, mutated and clipped exactly like q but with no effect. They share the cells' genealogy, so without selection the mean of q behaves like one more marker.
- **Detector.** L4 = L3 **and** dq = mean q − mean marker > 0.15 **and** z = dq / (spread of the 8 marker means) > 3, all measured over the functional cells at the end of the run. The marker spread measures how much drift this particular population had.
- **Why 8 markers.** A first version compared q with a single paired neutral marker at a fixed threshold. It gave 39% false positives at the reference and 26% over the space, because the mean of one marker in a single population can drift by 0.4–0.6 from genealogy alone. Using several markers measures that drift instead of assuming it.

## Validation (200 populations per setting; reference = stage 2b reference with η = 0.1)

**Adaptive rise at the reference.** Mean q (functional cells) climbs steadily, while the neutral markers barely move:

| step | 0 | 100 | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900 | end |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mean q | 1.00 | 1.11 | 1.37 | 1.71 | 2.04 | 2.32 | 2.59 | 2.87 | 3.15 | 3.40 | 3.65 |
| mean marker | 1.00 | 1.00 | 1.00 | 1.01 | 1.03 | 1.05 | 1.08 | 1.11 | 1.14 | 1.17 | 1.21 |

P(L3) = 1.00 and P(L4) = 0.86. The slight rise of the markers comes from clipping at 0, which affects q in the same way and therefore cancels in dq.

**Heritable variation.** With η = 0 nothing can evolve, and the detector correctly gives no L4:

| η | 0 | 0.02 | 0.05 | 0.1 | 0.2 |
|---|---|---|---|---|---|
| P(L4) | 0.00 | 0.26 | 0.68 | 0.85 | 0.87 |
| mean dq | 0 | +0.27 | +1.12 | +2.43 | +4.48 |
| mean z | – | 2.4 | 3.7 | 4.9 | 5.8 |

**Controls.**
- *No selection on q (double-neutral run: q has no effect either):* P(L4) = 0.01 at the reference (dq +0.03 ± 0.44, mean z 0.2). Over 4,096 random conditions it gives P(L4) = 0.002, which is a false-positive rate of 1.5% among L3 populations. This is the false-positive rate of the detector.
- *Where L3 fails* (S = 512 with σ = 2, or km = 0.1): no L4, by construction.
- *The `no_selection` knockout of stage 3a is not a "no selection on q" control.* That knockout only removes the preference for replacing dead cells. Cells with higher q still divide faster and so leave more descendants, which is selection on q. Accordingly it still gives P(L4) = 0.81 (dq +2.39). The proper control for selection on q is the double-neutral run above.
- *equal_split knockout:* P(L4) = 0.75. Evolution of q does not need random sharing at division, at least where L3 itself survives.

**Whole space** (4,096 random conditions; 1.02 ms per simulation): P(L1) 0.567, P(L2) 0.336, P(L3) 0.172, P(L4) 0.104. Among L3 populations, dq is +1.43 ± 1.21 with selection against −0.00 ± 0.33 in the double-neutral run (99th percentile 0.93), and P(L4 | L3) = 0.60.

L4 is limited by statistical power when η is small or the run has few generations, so near its boundary the outcome is genuinely stochastic. COMPASS treats this like any other noisy boundary.

## Map (COMPASS vs Sobol, 32,000 simulations each; ground truth 4,000 points × 64 seeds)

| level | life fraction (ground truth) |
|---|---|
| L1 | 0.558 |
| L2 | 0.338 |
| L3 | 0.173 |
| L4 | 0.112 |

The levels are perfectly nested in the ground truth (no point with L2 but not L1, L3 but not L2, or L4 but not L3). 64% of L3 conditions reach L4.

**Accuracy of the L4 map** (overall / on the 453 hard points near the boundary):

| simulations | 2k | 4k | 8k | 16k | 32k |
|---|---|---|---|---|---|
| COMPASS | 0.935 / 0.653 | 0.953 / 0.775 | 0.958 / 0.808 | 0.968 / 0.848 | **0.973 / 0.876** |
| Sobol | 0.942 / 0.691 | 0.949 / 0.751 | 0.954 / 0.786 | 0.961 / 0.817 | 0.961 / 0.815 |

COMPASS is behind at 2k (it has not yet found the boundary) and ahead from 4k on. At 32k it makes 6.1 percentage points fewer errors on hard points; Sobol stops improving between 16k and 32k. Engine overhead: 106 s for COMPASS against 37 s for Sobol over the whole run. The final map (an ensemble trained on all COMPASS data) has accuracy 0.969, hard accuracy 0.863, and a life fraction of 0.106.

**Importance** (permutation, drop in accuracy): km 0.111, μ 0.072, K 0.057, S 0.050, η 0.045, food 0.044, σ 0.013, d 0.012, dm 0.011; distractors 0.000.

**Read-offs from the map**
- *Heritable variation:* L4 needs η ≥ about 0.02. Averaged over the other conditions, P(L4) rises from 0 at η = 0 to 0.06 at η = 0.03, 0.12 at η = 0.1 and 0.145 at η = 0.2 (it levels off because the other conditions must also allow L3).
- *Mutation rate:* the critical parasite mutation rate for L4 at the reference is μ ≈ 0.041, lower than for L3 alone.
- *Food:* L4 needs food ≥ about 0.32 at the reference.
- *Cell size and parasites:* the largest division size that still gives L4 falls as parasites get stronger: S ≤ 256, 181, 147 and 169 for σ = 1.5, 2, 2.5 and 3 (no upper limit inside the range for σ ≤ 1.2). This is the stochastic corrector's limit, as in L2 and L3.
- *Food against heritable variation (empirical trade-off):* with less food, populations run fewer generations, so more variation per division is needed to show a clear rise. In the η × food slice at the reference, P(L4) passes 0.5 at η ≈ 0.02 for food 1.0, at η ≈ 0.03 for food 0.5, at η ≈ 0.07 for food 0.33 and only at η ≈ 0.13 for food 0.31; at food 0.29 it stays below 0.35 even at η = 0.2. We report this as an empirical shape of the map and do not fit a formula to it.

Figures: [accuracy](figures_3c/accuracy.png), [importance](figures_3c/importance.png), [partial dependence](figures_3c/pd_L4.png), slices [η × μ](figures_3c/slice_eta_x_mu.png), [η × S](figures_3c/slice_eta_x_S.png), [η × food](figures_3c/slice_eta_x_food.png), [S × σ](figures_3c/slice_S_x_sigma.png).

## What this means

1. **Evolution sits inside self-maintenance, heredity and coupling.** In this model there is no Darwinian evolution without a working protocell population, and about two thirds of working protocell conditions also evolve.
2. **The extra requirements are about information and time:** heritable variation must exist, copying errors that create parasites must stay rare (μ below about 0.04), and the population needs enough resources to run many generations.
3. **Detecting evolution needs an internal drift control.** A single neutral comparison is fooled by genealogy; eight neutral markers bring the false-positive rate to about 1%.

## Limits

- q affects all synthesis equally and saturates at h = 2, a deliberately simple choice. Other forms (a trade-off between speed and accuracy, or an effect on the membrane only) could change where L4 sits.
- L4 is measured over a fixed run length, so "no L4" near the boundary often means "not yet detectable" rather than "impossible".
- One model. Stage 3d checks the main conclusions on a spatial version.

## Raw data

`maps_l4_compact.json` (in [stage3/](../stage3/)) holds curves, importance, partial dependence and slices. Full logs and result files: [raw/stage3/](../raw/stage3/).
