# Paper outline: Stage 4, the environment search

**Status:** outline for review, 29 October 2026, Toronto time · Gokyona Wahiro (pen name) · Earlier reports: [stage 1](COMPASS_stage1.md), [stage 2](COMPASS_stage2.md), [stage 3](COMPASS_stage3.md)

## Working title

**How far up the ladder of life can environments carry a bottom-up chemistry? An ordered search over 13 environmental variables**

## Main argument

1. **Method.** An ordered, non-repeating active search (COMPASS) can map a fully bottom-up origin-of-life chemistry over up to 54 dimensions:
   - 13 environmental variables;
   - up to 41 unmeasured constants, kept as explicit unknowns.

   COMPASS found and confirmed specific environments for levels so rare that random samples of the same regions contained none. A first equal-budget check shows that this does **not** prove COMPASS finds rare points more often than random sampling would (see "Checks before writing"). Its proven advantage is map accuracy per simulation.
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

   The hardest single step is not catalysis or heredity alone. It is **heredity and stable membranes in the same environment**: 5 of the 7 confirmed heredity environments are acidic (pH 2.6–4.7). At pH below 4 no membrane in the model can form; at the other two, only 0.4–1.4% of the membrane-constant space allows one.
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
  - *First check (29 October), done before writing.* Hits per simulated (feasible) run:
    - box3 v2: the ground truth has 0 joint hits in about 2,130 feasible runs, so its 95% upper bound is about 1.4 × 10⁻³ per run. COMPASS has 5 in about 6,460, or 7.7 × 10⁻⁴.
    - box3: 0 in about 210 feasible ground-truth runs, against 2 in about 2,300 for COMPASS.

    COMPASS's hit rate is **inside** the range the ground truth allows, and the ground-truth runs come in groups of 16 seeds per point, which makes its bound even wider. **So the data do not show that COMPASS finds rare levels more often than random sampling would at equal budget.** The paper will say instead:
    - these levels are rare (below about 10⁻³ per run even in the focused regions);
    - COMPASS found and confirmed specific environments;
    - COMPASS's proven advantage is map accuracy per simulation (Part I, 4a, 4c-1), not a proven advantage in finding rare points.

  - *Pooled equal-budget test (29 October, hitrate4c.json, md5 0deabe99).* COMPASS is truncated to the same number of feasible runs as the ground truth in each map (4c-3 full, box2 rep 0, 4c-4 full, box3, box3 v2; 38,672 feasible runs each).
    - **L2_long:** COMPASS 4 runs, ground truth 0. If both had the same rate, all 4 hits would fall in COMPASS with probability 1/16, so p ≈ 0.06, one-sided.
    - **Joint:** 1 against 0.
    - **E3 strict in box3 v2:** 220 against 208 runs, the same per run.
    - **L2 strict:** COMPASS hits more distinct points, but most full-map L2 strict hits are noise-level short copies.
    - The ground truth tests far fewer distinct environments, about 2,400 points against about 38,700 for COMPASS, because it repeats each point with 16 seeds.

    **Wording for the paper:** "At equal budget, COMPASS found heredity with long copies 4 times and the joint level once, where the ground truth found neither (p ≈ 0.06). This is suggestive but not significant. Where a level is common enough to compare (E3 strict), hit rates per run are equal." The proven advantage remains map accuracy per simulation.
- **Re-check every number in the paper** against the result files (md5 listed in each stage report).
- **Verify every citation against the original**, especially those not yet opened: Millero et al. 2008 ion values, and Patel et al. 2015.
- **Upload code and raw data** from the benchmark machine (lot-s4/upload_* folders with MD5SUMS) before submission.
- **Disclose the use of AI** in the research and writing.

## Venue

**Decided by the owner, 29 October 2026: two papers.**
- **Paper 1, the COMPASS method:** Part I, plus engine scaling (4a) and engine performance on the bottom-up model. Targets: arXiv and ALIFE 2027.
- **Paper 2, Stage 4 (this outline):** the environment search. Target: a journal such as *Origins of Life and Evolution of Biospheres*, *Artificial Life* or *Life*.

Section 5.6 (engine performance) moves to paper 1 and is only summarised here.
