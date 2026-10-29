# Paper outline: the COMPASS method (paper 1)

**Status:** outline for review, 29 October 2026, Toronto time · Gokyona Wahiro (pen name) · Sources: [stage 1](COMPASS_stage1.md), [stage 2](COMPASS_stage2.md), [stage 3](COMPASS_stage3.md) papers, [stage 4a](../results/stage4a_scaling_b70.md), the stage 4 reports, and the companion [Stage 4 paper](COMPASS_stage4.md)

**Targets** (decided by the owner, 29 October 2026): an arXiv preprint (full version) and an ALIFE 2027 full paper (3–8 pages, a shortened version).

## Working title

**COMPASS: mapping where life-like behaviour appears in stochastic models, with an ordered, non-repeating active search**

## Main argument

1. **The problem.** Origin-of-life and artificial-life models are stochastic, slow and high-dimensional. Asking "under which conditions does a level of life appear?" means finding a boundary in condition space. Grid or random sampling wastes most simulations far from that boundary.
2. **The method.** COMPASS combines known parts:
   - Sobol coverage (Sobol' 1967; Owen 1998);
   - one simulation per condition, never repeated;
   - a shared neural map of P(level | conditions);
   - selection by disagreement between bootstrap ensemble members (query by committee, Seung et al. 1992; deep ensembles, Lakshminarayanan et al. 2017).

   We do not claim the parts as new. Our contribution is the combination for this question, a benchmark protocol (every method gets the same map-maker and is scored against brute-force ground truth), and an honest record of what did not work.
3. **Result: map accuracy per simulation.** At equal accuracy, COMPASS needs fewer simulations than Sobol or random sampling in every model tested:

   | Model | Dimensions | Saving against Sobol or random | Source |
   |---|---|---|---|
   | Toy protocell (chemostat) | 10 | about 8× | stage 1 |
   | RAF sets, binary polymer model | 7 | about 2× | stage 1e |
   | Stochastic corrector, L2 | 8 | about 8× | stage 2a |
   | Coupled protocell, L3 | 10 | about 4× | stage 2b |
   | Coupled protocell padded with extra dimensions | 20–30 | 2.8–4.0× | stage 4a |
   | Bottom-up chemistry, E1 / E2 | 29 | 2× / 1.6× to Sobol's final accuracy; 2.4× / 2.9× to 95% of the achievable gain | stage 4c-1, curve analysis |

   The saving is largest where the boundary is thin and smallest where most of the space is plainly negative.
4. **Result: maps agree with theory.**
   - The RAF threshold of Hordijk, Kauffman & Steel (2011) is reproduced to within 2%.
   - The threshold read off a COMPASS map with 32,000 simulations lies within 7% of direct simulation.
   - Thresholds in the stage 2 maps agree with validation runs.
5. **Result: what matters in the engine** (stage 3e, 2d, 4a).
   - A standard boundary rule (straddle, Bryan et al. 2005) does as well as COMPASS's own scoring. The saving comes from the loop as a whole, not from the formula.
   - A joint, nested map of all levels removes contradictions between levels at no cost in accuracy, and runs about 3× faster.
   - Multi-fidelity screening is a partial negative result: it helps only at small budgets, because the cheap model is biased at the boundary.
   - The engine's cost does not grow with dimension, at about 75 s per 64,000-simulation run.
6. **Result: rare levels in a 54-dimensional bottom-up model** (companion paper).
   - COMPASS found, and confirmed with new seeds, environments for levels that the ground truth never hit: long-copy heredity, and heredity inside vesicles.
   - At an equal budget this is suggestive, not significant (4 against 0 and 1 against 0, p ≈ 0.06).
   - **COMPASS does not give a uniform gain in hits per run.** Compared with a single-seed random sampler, it hits more for some levels (heredity in map4c3, 32 against an expected 15), fewer for others (catalytic networks in a focused box, 4–6 against about 20). Most of its raw advantage over the ground truth comes from the ground truth's 16 seeds per point. The engine targets uncertainty at the boundary, not hits.
7. **Limits.**
   - Savings depend on boundary thickness and are smaller in mostly negative spaces.
   - When simulations are cheap, the engine's own cost matters.
   - Only Sobol and random sampling are compared at full scale; straddle only in stage 3e.
   - Most comparisons have 1–3 repetitions.

## Structure (ALIFE version, about 8 pages)

1. Introduction: maps instead of single histories; why the boundary matters.
2. Related work: active learning and level-set estimation (Settles 2009; Bryan et al. 2005; Gotovos et al. 2013; Houlsby et al. 2011), quasi-random designs, and models of the origin of life (autocatalytic sets, stochastic corrector).
3. The COMPASS loop: algorithm box, defaults, the joint nested map.
4. Benchmark protocol: identical map-maker, brute-force ground truth, hard sets near the boundary, compute overhead.
5. Results: the savings table; scaling with dimension; agreement with theory; the engine ablations.
6. Use on a bottom-up model: rare levels, and the hit-rate analysis with its design effect.
7. Limits and negative results.
8. Code, data and AI use.

The arXiv version adds the per-model sections, the negative first benchmark of stage 1, the full multi-fidelity result and all accuracy curves.

## Figures (planned)

1. The COMPASS loop (diagram).
2. Accuracy against budget, COMPASS against Sobol, one panel per model.
3. Saving against dimension (stage 4a).
4. A map with its known threshold overlaid (RAF).
5. Ablation: COMPASS scoring, straddle, Sobol.

## Checks before writing

- Take every number from the stage reports and re-check it against the result files, as for paper 2.
- Upload lot-s4/upload_curves (the 2.4× / 2.9× analysis) so that the number can be reproduced.
- Confirm the ALIFE 2027 call (dates and page limit) once it is posted.
- arXiv needs an endorser in the category: find a domain contact.
- Disclose AI use.
