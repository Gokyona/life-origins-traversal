# Stage 3 plan: what is necessary, what is sufficient, and when does evolution start

**Status:** draft for review, 25 October 2026 · **Engine:** COMPASS (see the [stage 1](../paper/COMPASS_stage1.md) and [stage 2](../paper/COMPASS_stage2.md) papers)

## 1. Goal

Stage 2 produced maps: for a protocell model with ten conditions, where in the condition space do self-maintenance (L1), heredity through division (L2) and a coupled protocell (L3) appear. A map says *where* life appears, but not *why*, nor what is the least that must be true. Stage 3 asks three questions:

1. **Necessary mechanisms.** Which ingredients of the model are needed at all? If we switch off a mechanism (compartments, random sharing at division, competition between cells, the membrane coupling) does the life region shrink, move, or disappear?
2. **Minimal sufficient conditions.** What is the smallest set of conditions that must be controlled, and within which ranges, so that life appears with high probability whatever the other conditions do? This is the "knockout → minimal sufficient condition set" step of the original proposal.
3. **The start of evolution (L4).** Under which conditions does a protocell population not only persist but improve, i.e. undergo Darwinian evolution?

Two extensions support these questions: a spatial version of the model (local competition), and engine work so that nested levels are mapped consistently.

## 2. Steps

**3a. Mechanism knockouts (on the stage 2b model).** Build switch-off variants of `protocell.py` and map each with COMPASS at L1–L3, using the same ground-truth protocol:
- *no compartments*: all molecules in one large well-mixed pool (one cell of the whole population size);
- *no random sharing*: division splits every molecule type exactly in half (removes the chance variation the stochastic corrector relies on);
- *no selection between cells*: the daughter replaces a random cell regardless of whether cells are alive or growing;
- *no membrane coupling*: division by content size as in 2a (the L3 coupling switched off).
Output: life fraction per level for each variant, and maps showing where life is lost or gained. Expected: no compartments and no random sharing should remove protection against parasites (life only at σ ≈ 1 or μ ≈ 0); the test is whether the maps show this.

**3b. Minimal sufficient conditions.** From the stage 2 maps, and confirmed by direct simulation:
- for each level, find the largest box in condition space (ranges on as few conditions as possible) inside which P(life) ≥ 0.95, with every other condition free over its full range;
- greedy search: start with all conditions free, constrain the condition that raises the guaranteed probability most, repeat until 0.95 is reached; record the order and the ranges;
- confirm each box with fresh direct simulations at random points inside it;
- compare the sets for L1, L2 and L3: which extra conditions does each level need?
Output: a short, checkable statement per level, for example "L3 is guaranteed if km is in [a, b], μ < c and S in [d, e], whatever the other conditions".

**3c. Darwinian evolution (L4).** Give every gene copy a heritable quality q (for example its copying speed or catalytic efficiency) that mutates by small random steps at copying. Over a long run, test whether mean q in functional cells rises beyond what drift alone produces (compare against a control run in which selection on q is switched off). L4 = L3 plus a significant adaptive increase in q. Map L4 with COMPASS. This also tests where heredity is accurate enough for selection to act, which connects to the critical mutation rate found in stage 2.

**3d. Spatial version.** Cells on a two-dimensional grid; a daughter replaces a neighbour instead of a random cell; optionally molecules leak between neighbours. Compare well-mixed and spatial maps (compare Hogeweg & Takeuchi 2003). Simulations will cost more, so this is also where multi-fidelity screening (stage 2d) gets a second chance.

**3e. Engine.** A joint map for nested levels that enforces L3 ⊂ L2 ⊂ L1 (removes the 0.4% inconsistency seen in stage 2), and a comparison of COMPASS with one standard active-learning method (BALD or level-set estimation) on the stage 2 benchmarks.

**3f. Stage 3 paper and homepage.** References verified against originals.

## 3. Milestones and acceptance

| step | output | accept when | rough effort |
|---|---|---|---|
| 3a | knockout variants, maps L1–L3 | each variant validated by direct sweeps; effect on life fraction measured with confidence intervals | 1–2 days |
| 3b | minimal sufficient condition sets for L1–L3 | each box confirmed by fresh direct simulations (P ≥ 0.95 within error) | 1–2 days |
| 3c | evolving protocells, L4 detector, L4 map | adaptive increase clearly above the no-selection control in a reference setting | about 1 week |
| 3d | spatial model, spatial maps | spatial model validated; differences from well-mixed maps quantified | about 1 week |
| 3e | joint nested map; BALD or level-set comparison | nesting violations eliminated without loss of accuracy | a few days |
| 3f | stage 3 paper | references verified | a few days |

## 4. Decisions needed

1. **Order.** Start with 3a + 3b on the existing model (recommended: cheap, builds directly on the stage 2 maps, answers "what is necessary and sufficient" first), or start with L4 evolution, or with the spatial model.
2. **Scope.** Include all of 3a–3e in stage 3 (recommended), or stop after 3a–3c and move spatial structure to stage 4.
3. **Volunteer computing.** Keep it for after stage 3 (recommended), or design it now.

## References

- Hogeweg, P. & Takeuchi, N. (2003). Multilevel selection in models of prebiotic evolution: compartments and spatial self-organization. *Origins of Life and Evolution of the Biosphere* 33, 375–403.

Further references (for L4 criteria and active-learning baselines) will be verified against the originals before they are used.
