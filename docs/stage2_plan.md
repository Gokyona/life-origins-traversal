# Stage 2 plan: from self-sustaining chemistry to dividing protocells

**Status:** draft for review, 25 October 2026 · **Engine:** COMPASS (see the [stage 1 paper](../paper/COMPASS_stage1.md))

## 1. Goal

Stage 1 showed that COMPASS maps the boundary between "life" and "no life" with fewer simulations than random or Sobol sampling, and that its map recovers a known threshold (autocatalytic sets, level L1). Stage 2 moves one level up the life detector. The question is: **under which conditions does a population of protocells keep growing, dividing and passing on its working molecules, instead of being destroyed by parasites, dilution or errors?**

The output is a set of condition maps for levels L1–L3, checked against what theory already says about compartment models, and a measured saving from COMPASS on a model where each simulation is expensive.

## 2. Life detector for stage 2

Each level is a yes/no test computed automatically at the end of a simulation:

| level | meaning | test (first version) |
|---|---|---|
| L0 | ordinary chemistry | none of the below |
| L1 | self-maintaining chemistry | the essential molecules persist in at least one cell to the end of the run |
| L2 | heredity through division | the population survives G generations of division and the essential gene types are kept in a majority of cells (parasites do not take over) |
| L3 | coupled protocell | as L2, and cell growth and division are driven by the cell's own chemistry (membrane made inside the cell), not imposed from outside |
| L4 | Darwinian evolution | deferred to stage 3 |

## 3. Model, built in steps

**2a. Stochastic corrector population (well mixed).** Start from the stochastic corrector model (Szathmáry & Demeter 1987; re-examined by Grey, Hutson & Szathmáry 1995). Each protocell contains a few types of mutually dependent template replicators and a parasite (a replicator that copies faster but contributes nothing). Templates copy using food, with errors that turn genes into parasites. A cell divides when its content reaches a size threshold, and molecules are shared at random between the two daughters. Cell growth rate depends on having a balanced set of gene types. The population has a fixed capacity. One simulation is one population followed for many generations, batched on the GPU like stage 1.

*Acceptance:* reproduce known qualitative behaviour. Group selection keeps parasites in check only when cells are small enough; with large cells parasites win; with very small cells genes are lost by chance at division (see e.g. Szathmáry & Demeter 1987; the trade-off is also discussed in later work on compartment models such as Hogeweg & Takeuchi 2003 and Fontanari, Santos & Szathmáry 2006). There is a mutation rate above which the population collapses.

**2b. Coupling chemistry, replication and membrane.** Add a membrane component made inside the cell by a reaction that one gene type catalyses, so that division follows from the cell's own growth (the idea behind Gánti's chemoton). Optionally use the stage 1e RAF network as the cell's metabolism, linking L1 and L3.

**2c. Spatial version (later).** Cells on a grid with local competition and diffusion, to test how spatial structure changes the maps (compare Hogeweg & Takeuchi 2003). Planned after the well-mixed version works.

## 4. Condition space (8–12 conditions)

Candidates: food supply, dilution or death rate, copying fidelity, parasite advantage, division size, number of essential gene types (1–3), membrane synthesis rate, wet–dry cycling or temperature, catalysis level f (if the RAF core is used), plus two distractors with no effect, as in stage 1.

## 5. Engine work

- Default engine: COMPASS with soft selection and batched map training (`engine4_soft --fast`).
- **Multi-fidelity.** Screen conditions with short runs on small populations, then confirm near the boundary with long runs. This is the main new engine feature for stage 2 and is where expensive models should benefit most.
- **Ground truth without brute force everywhere.** 128 runs per test condition will be too expensive. Plan: a smaller held-out set (about 2,000 conditions × 32 runs), plus slice checks against known behaviour, as in the stage 1e threshold analysis.
- Optional: compare COMPASS with one standard active-learning method (BALD or level-set estimation) on the same benchmark.

## 6. Milestones and acceptance criteria

| step | output | accept when | rough effort |
|---|---|---|---|
| 2a | stochastic corrector population model + L1/L2 detector | reproduces the qualitative results in section 3; cost per simulation measured | about 1 week |
| 2b | membrane coupling, L3 detector | division emerges from internal chemistry in a reference setting | about 1 week |
| 2d | condition space and multi-fidelity COMPASS | engine runs end to end on the new model | a few days |
| 2e | maps for L1, L2, L3 | COMPASS saves at least 2× simulations vs Sobol at equal accuracy; maps agree with section 3 checks | 1–2 days of GPU time |
| 2f | stage 2 paper, homepage update | references verified against originals | a few days |

Total: about 4–6 weeks at a part-time pace. The largest uncertainties are the coupling design (2b) and the cost of one simulation (measured in 2a).

## 7. Decisions needed

1. **Starting model:** build on the stochastic corrector model (recommended, results can be checked against the literature) or design a new protocell model from scratch.
2. **Scope:** stop at L3 (recommended) or include L4, Darwinian evolution.
3. **Space:** well-mixed population first, spatial grid later (recommended), or spatial from the start.

## References

- Fontanari, J. F., Santos, M. & Szathmáry, E. (2006). Coexistence and error propagation in pre-biotic vesicle models: a group selection approach. *Journal of Theoretical Biology* 239, 247–256.
- Grey, D., Hutson, V. & Szathmáry, E. (1995). A re-examination of the stochastic corrector model. *Proceedings of the Royal Society of London B* 262, 29–35.
- Hogeweg, P. & Takeuchi, N. (2003). Multilevel selection in models of prebiotic evolution: compartments and spatial self-organization. *Origins of Life and Evolution of the Biosphere* 33, 375–403.
- Szathmáry, E. & Demeter, L. (1987). Group selection of early replicators and the origin of life. *Journal of Theoretical Biology* 128, 463–486.

Citation details above were checked against Szathmáry's publication list and the Utrecht University research portal. The specific claims attributed to each paper will be checked against the full texts before the stage 2 paper.
