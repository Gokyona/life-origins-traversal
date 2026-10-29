# Paper outline: Stage 4, the environment search

**Status:** outline for review, 29 October 2026, Toronto time · Gokyona Wahiro (pen name) · Earlier reports: [stage 1](COMPASS_stage1.md), [stage 2](COMPASS_stage2.md), [stage 3](COMPASS_stage3.md)

## Working title

**How far up the ladder of life can environments carry a bottom-up chemistry? An ordered search over 13 environmental variables**

## Main argument

1. **Method.** An ordered, non-repeating active search (COMPASS) can map a fully bottom-up origin-of-life chemistry over up to 54 dimensions:
   - 13 environmental variables;
   - up to 41 unmeasured constants, kept as explicit unknowns.

   COMPASS finds rare levels that random sampling of the same region does not find. This claim needs an equal-budget check (see "Checks before writing").
2. **Result: the ladder.** In this model, environments can carry chemistry from activated monomers to coupled protocells (E1 → E2 → L2 → E3 → L3). Each step narrows the set of environments:

   | Level | Share of environments where it occurs |
   |---|---|
   | Activated feedstock | about 45% of the space |
   | Long chains | about 10% |
   | Self-sustaining catalytic networks | about 0.3%, and only with highly catalytic short chains |
   | Heredity (template copying) | a narrow class; 5 confirmed environments in the 4c-3 maps |
   | Heredity inside stable vesicles | 4 environments found by the search itself; 0 in random samples |
   | Coupled protocells | wherever heredity and vesicles coexist (298 of 435 coexistence runs) |

3. **Result: the environment.** Every environment that reaches the top is a cool to mild pool (0–25 °C) with these features:
   - temperature cycles of days to weeks across the melting point of short duplexes;
   - slow washout, with molecules living long compared with a cycle;
   - accurate copying (fewer than about 0.3 errors per copy);
   - near-neutral to slightly acidic water (pH 4.4–7.6), with membranes stable there.

   The hardest single step is not catalysis or heredity alone. It is **heredity and stable membranes in the same environment**: half of the heredity environments found are acidic or warm, where membranes in the model cannot form or rarely hold.
4. **Result: the ceiling.** The ladder stops before Darwinian evolution.
   - Selection among vesicles happens as one fast sweep.
   - After it, copying refills every vesicle to the same load, so no heritable variation remains.
   - This is a property of the model (sequences have no phenotype), not of environments. It gives a sharp prediction: **sustained evolution in such a system needs a sequence-dependent trait that copying does not erase.**
5. **What is built in and what is derived.** Every result is classified as built in, partly built in or derived (as in each stage report). The paper is honest about which conclusions follow from assumptions.

## Structure

1. **Introduction.** The phase diagram of life, rather than one lucky history. Why environments, not system properties (the limit of Part I). The ladder E1–E3, L1–L4.
2. **Model.** The bottom-up chemistry, layer by layer (4c-1 to 4c-4), with each rate law and its source:
   - feedstock and activation;
   - condensation and hydrolysis;
   - catalysis;
   - template copying with errors and duplex melting;
   - amphiphiles, the membrane window and vesicles.

   Also covered: layer A (13 variables, ranges with sources, from [stage4a_ranges.md](../docs/stage4a_ranges.md)) and layer C (unmeasured constants as dimensions).
3. **Search method.** COMPASS in brief (cite the stage 1–3 reports), scaling to 54 dimensions (4a), focused maps with masking, confirmation with 16 runs, and paired random-copy and β = 0 controls.
4. **Validation.** One table of all pre-registered tests per layer. It includes the failures (4c-3 test 9 at #1558; 4c-5 tests 3 and 6) and the design changes made before results.
5. **Results.**
   - 5.1 The ladder: fractions per level, full and focused maps.
   - 5.2 The environments at the top: the 4 search-found coexistence points and the 4c-3 heredity points, in a table of environmental conditions.
   - 5.3 The environmental sequence (4d-seq): how each condition must change from level to level, with conflicts marked.
   - 5.4 Real settings (4e): which known site types match the environments found.
   - 5.5 Coupled protocells and the ceiling at L4.
   - 5.6 Engine performance on the bottom-up model.
6. **What follows from the assumptions:** the combined built-in or derived table.
7. **Discussion.**
   - Comparison with the literature: hot-spring wet–dry cycling (Damer & Deamer 2020), alkaline vents (Martin et al. 2008; work by Lane and colleagues, references to be chosen and verified), cold or icy settings (eutectic freezing), thermal gradients in pores (Baaske et al. 2007), and fatty-acid membranes and Mg²⁺ (Mansy & Szostak 2008; Adamala & Szostak 2013; Jordan et al. 2019).
   - Where our model agrees, disagrees, or cannot say.
8. **Limits.**
   - A 2-letter alphabet and a 64-letter cap.
   - Coarse membrane laws.
   - Copying does not depend on Mg²⁺.
   - No sequence phenotype.
   - Layer C ranges decide much (for example pKa_app near 6 at three of the four top points).
   - The focused regions were chosen from earlier maps.
   - No full joint map (4d) and no environment-change simulation.
9. **Predictions and next steps.** Testable predictions for experiments, a new layer with sequence-dependent traits, and other planetary settings (stage 5).

## Figures (planned)

1. The ladder: the fraction of the space per level, full against focused maps.
2. Layer A conditions of the top environments: radar or parallel-coordinates plot.
3. The environmental sequence (4d-seq): per-variable ranges by level.
4. Heredity inside vesicles, with a snapshot of the lineage enclosed share and the vesicle family tree.
5. The L4 ceiling: test 6b (founder sweep) and the equalised load of test 6c.
6. COMPASS against ground truth: hits found per budget in the focused regions.

## Checks before writing

- **Equal-budget comparison.** "Random sampling finds none" compares unequal budgets: 8,000 ground-truth runs (500 points × 16 seeds) against 16,000 COMPASS runs. We need an equal-budget statement, including a Poisson bound on the rate the ground truth could have missed, or the claim must be softened.
- **Re-check every number in the paper** against the result files (md5 listed in each stage report).
- **Verify every citation against the original**, especially those not yet opened: Millero et al. 2008 ion values, and Patel et al. 2015.
- **Upload code and raw data** from the benchmark machine (lot-s4/upload_* folders with MD5SUMS) before submission.
- **Disclose the use of AI** in the research and writing.

## Venue (to decide with the owner)

The earlier plan was two papers:
- paper 1, the COMPASS method (Part I), for arXiv and ALIFE 2027;
- paper 2, Stage 4, for a journal such as *Origins of Life and Evolution of Biospheres*, *Artificial Life* or *Life*.

This outline is paper 2. It can also stand alone, with Part I summarised in the method section.
