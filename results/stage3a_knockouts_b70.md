# Stage 3a: mechanism knockouts of the protocell model (Intel Arc Pro B70)

**Date:** 2026-10-25 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Run:** 17:41–18:17 Toronto time (36 min) · **Code:** [`stage3/knockouts.py`](../stage3/knockouts.py), [`stage3/validate_knockouts.py`](../stage3/validate_knockouts.py), [`stage3/s3a_run.py`](../stage3/s3a_run.py) · **Plan:** [stage3_plan.md](../docs/stage3_plan.md)

## Bottom line

In the stage 2b protocell model, **compartments** and **selection between cells** are the mechanisms life cannot do without. Without compartments, life disappears from most of the condition space and the coupled protocell (L3) is impossible; without selection between cells, heredity (L2) is lost in 11% of the space and nothing is gained. Random sharing of molecules at division has almost no net effect over the whole space, but it is decisive in one dangerous region (large cells, strong parasites). Coupling division to the membrane is not needed for heredity but defines the protocell level. Degradation and extra gene types are pure burdens.

## Knockouts

Each knockout removes one mechanism from [`protocell.py`](../stage2/protocell.py); everything else, including the ten conditions, is unchanged.

| knockout | what is removed |
|---|---|
| `equal_split` | random sharing at division: every molecule type is split exactly in half (odd counts rounded at random) |
| `no_compartments` | cells: one well-mixed soup per population with the same total starting material, crowding against a fixed capacity, and a chemostat outflow that keeps molecules turning over; no membrane, lysis or division; L3 impossible |
| `no_selection` | selection between cells: a daughter replaces a random cell, whether it is dead or alive |
| `no_coupling` | the membrane: division at content size S as in stage 2a; L3 undefined |
| `no_decay` | losses: degradation d = 0 and membrane turnover dm = 0 |
| `K1` | complexity: a single gene type does all jobs |

The soup in `no_compartments` is a modelling choice: with a different outflow rule its numbers would change, so we read it qualitatively.

## Validation (direct simulation, 200 populations per point)

P(L3) along division size S at parasite advantage σ = 2 (other conditions at the 2b reference); for `no_coupling` P(L2), for `no_compartments` P(L2) of the soup:

| knockout | S = 8 | 16 | 64 | 128 | 256 | 512 |
|---|---|---|---|---|---|---|
| none | 1.00 | 1.00 | 1.00 | 1.00 | 0.25 | 0 |
| equal_split | 1.00 | 1.00 | 1.00 | 0.94 | **0.01** | 0 |
| no_compartments | 0 | 0 | 0 | 0 | 0 | 0 |
| no_selection | **0** | 0.81 | 1.00 | 1.00 | 0.21 | 0 |
| no_coupling (L2) | 1.00 | 1.00 | 1.00 | 0.00 | 0 | 0 |
| no_decay | 1.00 | 1.00 | 1.00 | 1.00 | 0.20 | 0 |
| K1 | 1.00 | 1.00 | 1.00 | 1.00 | **0.91** | 0 |

The soup survives only when parasites have little or no advantage (σ ≤ 1.2; L1 = L2 = 1.00 at σ = 1 for K = 1–3). Without selection, small cells fail: cells that lose genes at division are no longer replaced preferentially.

## Life regions

Ground truth: the same 4,000 random conditions × 64 simulations for every variant. Life fraction per level, with the share of the whole space gained / lost relative to the unchanged model:

| knockout | L1 self-maintenance | L2 heredity | L3 coupled protocell |
|---|---|---|---|
| none | 0.564 | 0.310 | 0.162 |
| equal_split | 0.575 (+0.026 / −0.015) | 0.325 (+0.026 / −0.012) | 0.167 (+0.019 / −0.014) |
| no_compartments | 0.159 (+0.039 / **−0.445**) | 0.092 (+0.045 / **−0.264**) | impossible |
| no_selection | 0.538 (+0.002 / −0.028) | 0.201 (+0.000 / **−0.109**) | 0.112 (+0.000 / −0.050) |
| no_coupling | 0.601 (+0.164 / −0.127) | 0.442 (+0.232 / −0.100) | undefined |
| no_decay | 0.773 (+0.219 / −0.009) | 0.463 (+0.155 / −0.002) | 0.288 (+0.136 / −0.009) |
| K1 | 0.783 (+0.224 / −0.005) | 0.488 (+0.179 / −0.001) | 0.185 (+0.070 / −0.046) |

COMPASS maps of every variant (16k simulations per level) reach 0.96–0.995 accuracy on all test conditions (hard sets 0.55–0.87), so the engine stays reliable when the mechanisms change.

| knockout | L1 all / hard | L2 all / hard | L3 all / hard |
|---|---|---|---|
| none | 0.966 / 0.775 | 0.981 / 0.663 | 0.979 / 0.757 |
| equal_split | 0.965 / 0.764 | 0.981 / 0.688 | 0.971 / 0.549 |
| no_compartments | 0.990 / 0.840 | 0.988 / 0.759 | – |
| no_selection | 0.962 / 0.774 | 0.985 / 0.800 | 0.984 / 0.753 |
| no_coupling | 0.976 / 0.848 | 0.984 / 0.867 | – |
| no_decay | 0.976 / 0.786 | 0.983 / 0.808 | 0.970 / 0.755 |
| K1 | 0.977 / 0.819 | 0.995 / 0.855 | 0.993 / 0.817 |

## Does random sharing at division matter?

Between-cell variance of the gene-0 fraction among functional cells at the end (mean over 200 populations):

| setting | unchanged | equal_split |
|---|---|---|
| reference (S = 64, σ = 1.5) | 0.0507 | 0.0431 |
| S = 256, σ = 2 | 0.0174 (P(L3) 0.32) | 0.0012 (P(L3) 0.01) |

At the reference, copying and mutation noise alone produce most of the variation between cells, so removing random sharing changes little. In large cells with strong parasites, random sharing is what creates the variation that selection acts on; without it, variation nearly vanishes and so does life. (The equal_split value at S = 256 rests on few functional cells.) Over the whole space the gains and losses roughly cancel (+1.9% / −1.4% at L3), because random sharing helps against parasites in large cells and costs genes in very small ones. In this model, then, the core of the stochastic corrector is **compartments plus selection between them**; random sharing is an extra source of variation that becomes essential where other variation is too weak.

## Findings

1. **Compartments are essential.** Without them L1 loses 45% of the space and L2 26%, and L3 cannot exist. A shared soup is taken over by parasites unless they have almost no advantage.
2. **Selection between cells is the second essential mechanism.** Removing it loses 11% (L2) and 5% (L3) of the space and gains nothing.
3. **Random sharing at division** has little net effect but is decisive for large cells with strong parasites.
4. **Membrane coupling** is not needed for heredity (removing it even enlarges the L2 region, by removing the km constraint) but is what makes a coupled protocell (L3) possible.
5. **Degradation and extra gene types** are burdens: removing them enlarges every life region.

## Raw data

`s3a.json` and `s3a.log`: [raw/stage3/](../raw/stage3/).
