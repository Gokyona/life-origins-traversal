# How far up the ladder of life can environments carry a bottom-up chemistry? An ordered search over 13 environmental variables

**Draft 1, for review** · Gokyona Wahiro (pen name), independent researcher, Toronto · 29 October 2026

Code and data: <https://github.com/Gokyona/life-origins-traversal> · Text: CC BY 4.0 · Code: MIT · Method: [stage 1](COMPASS_stage1.md), [stage 2](COMPASS_stage2.md), [stage 3](COMPASS_stage3.md) reports · Outline and checks: [stage4_outline.md](stage4_outline.md)

> **Status of this draft.** All numbers come from the stage reports listed in each section. Before submission, every number will be checked again against the result files (md5 sums are given in the reports), and every citation against its original. Places that still need work are marked **[TODO]**.

---

## Abstract

Which environments can turn non-living chemistry into something that inherits, encloses itself and evolves? We built a fully bottom-up artificial chemistry with nothing living put in. It runs from feedstock through activated monomers, chains, catalysts and template copying to amphiphiles and vesicles. Its only inputs are 13 environmental variables, encoded as 19 dimensions with ranges taken from the literature. Rate constants that no measurement fixes are kept as explicit unknown dimensions, up to 54 dimensions in all. We searched this space layer by layer with COMPASS, an ordered, non-repeating active-learning engine, and detected seven levels automatically against paired controls.

The model climbs the ladder up to coupled protocells, and each step narrows the environments sharply:
- activated feedstock appears in about 45% of the space;
- persistent chains in about 10–16%;
- self-sustaining catalytic networks only when short chains are highly catalytic (0.3% of a focused region);
- heredity by template copying in a narrow class of environments (5 confirmed);
- heredity inside stable vesicles in 4 environments found by the search;
- coupled protocells wherever heredity and vesicles coexist (298 of 435 coexistence runs).

All environments at the top are cool to mild (0–25 °C), slowly flushed, enclosed water bodies with temperature cycles of days to weeks across the melting point of short duplexes, and near-neutral to slightly acidic water. The ladder shows no conflict between consecutive levels: one stable environment of this kind could carry the chemistry through all levels reached. The hardest step is having heredity and stable membranes in the same place. The ladder stops before Darwinian evolution. Copying families sweep once, after which copying erases all heritable variation in vesicle content, because sequences in the model have no phenotype. We therefore predict that sustained evolution in such systems needs a sequence-dependent trait that copying does not erase. For every result we state whether it follows from the model's assumptions or is derived.

## 1. Introduction

Most simulations of the origin of life choose one set of conditions and wait for something life-like to appear. A negative result then says little: the conditions may have been wrong, or the run too short. We ask the opposite question: **in which environments does a given chemistry climb the ladder of life, and where does it stall?** We answer it with a map rather than a single history.

Part I of this project ([stage 1–3 reports](COMPASS_stage3.md)) built the mapping engine, COMPASS, and validated it on test models. There it needed 2–8 times fewer simulations than random or Sobol sampling for the same map accuracy. Part I had one fundamental limit: its conditions described the chemical system (cell size, error rate, parasite advantage), not the environment. This paper reports Part II. Its inputs are environmental variables only, and its chemistry is built from the bottom up, so that every level of life must emerge.

We use an extended ladder of seven levels, each scored automatically:
- **E1:** activated feedstock reaches a working concentration.
- **E2:** chains of at least 8 units persist.
- **L1:** a self-sustaining catalytic network (a persistent RAF set, Hordijk & Steel 2004) that raises its own products above a paired run without catalysis.
- **L2:** heredity. A sequence lineage outlives its molecules because it is copied, beating two paired random-copy controls, with copies of at least 6 letters (L2_long).
- **E3:** persistent vesicles that hold chains.
- **L3:** a coupled protocell. Vesicle growth depends on content, and the content is passed on by copying.
- **L4:** Darwinian evolution. A heritable vesicle trait rises beyond drift.

**[TODO]** Short review of scenario literature (hot-spring pools, vents, ice, lakes) and of earlier artificial-chemistry maps.

## 2. Model

The simulator was built in five layers. Each layer was validated before the next was added, and none changes the layers below. Full specifications: [4c-1](../docs/stage4c1_spec.md), [4c-2](../docs/stage4c2_spec.md), [4c-3](../docs/stage4c3_spec.md), [4c-4](../docs/stage4c4_spec.md), [4c-5](../docs/stage4c5_spec.md).

**Environment (layer A).** There are 13 variables, encoded as 19 dimensions:
- temperature, pressure, pH, H₂ availability and UV flux;
- salinity and Mg²⁺;
- feedstock supply and composition;
- residence time and mixing between sites;
- evaporation and thermal-gradient concentration;
- wet–dry and hot–cold cycling (depth and period);
- mineral surfaces (amount and type);
- available time.

Ranges and their sources are in [stage4a_ranges.md](../docs/stage4a_ranges.md). They run from acid springs to alkaline vents and from flowing channels to closed pools.

**Unknown constants (layer C).** Constants that no measurement fixes, such as condensation rate, catalytic probability, copy speed-up, error rate, duplex melting temperature and membrane stability limits, are dimensions of the map with literature-based ranges. Their number grows from 10 (4c-1) to 41 (4c-4).

**Layers:**
1. **Polymer chemistry (4c-1).**
   - Feedstock becomes monomers, which are activated by UV, H₂ or drying.
   - Activated monomers start and lengthen chains, which hydrolyse, wash out, adsorb to minerals and break under UV.
   - Hydrolysis follows Li & Breaker (1999) for nucleotide-like chains and Radzicka & Wolfenden (1996) for peptide-like chains.
   - Each run has 64 reaction sites with 512 chain slots each.
2. **Catalysis (4c-2).** Any chain of length ≥ L_cat may, by a fixed random rule, catalyse the joining of two other molecules.
3. **Template copying (4c-3).**
   - Nucleotide-like chains of a 2-letter alphabet (0/1) are copied letter by letter into their complement, with an error rate ε.
   - Copies stay paired until the duplex melts, with a two-state duplex energy set by the melting temperature of an 8-pair duplex, T_m8.
   - Chains are capped at 64 letters.
4. **Amphiphiles and vesicles (4c-4).**
   - Amphiphiles from the feedstock assemble above a critical concentration.
   - They form vesicles only inside a stability window: |pH − pKa_app| ≤ w_pH, Mg²⁺ and salinity below limits, and the hot phase below T_mem.
   - Vesicles enclose nearby chains at random at formation and at rehydration.
   - They grow faster when filled, by osmotic pressure through a factor β (Chen, Roberts & Szostak 2004).
   - They divide under shear or cycling (Zhu & Szostak 2009), sharing contents at random.
   - Chains of 4 or more letters cannot leave a vesicle (Mansy & Szostak 2008).
5. **Coupling (4c-5).** No new rules. Vesicle lineages are recorded and L3 and L4 are only detected.

## 3. Search and detection

**COMPASS.**
- An ensemble of neural classifiers is trained on all runs so far.
- New environments are drawn from a Sobol sequence and selected where the ensemble disagrees most. Each environment is tested once.
- At 20–30 dimensions, the engine still needed about 3–4× fewer simulations than Sobol sampling for the same accuracy (stage 4a).

**Maps per layer.**
- **Full map:** 1,000 random environments × 16 seeds as ground truth, and 32,000 COMPASS runs.
- **Focused maps:** regions chosen from the full map. Environments outside a region's constraints are masked, not simulated. Each has 500 × 16 ground truth and 16,000 COMPASS runs.

**Confirmation.** Every candidate for a rare level was re-run with 16 new seeds.

**Controls.**
- L1: a paired run without catalysis.
- L2: two paired runs where copying makes chains of the same length with random sequence.
- L3: paired runs with β = 0 (no osmotic coupling) and with random-sequence copying.

**Pre-registration.** Every layer had validation tests with criteria written before they ran. Design changes made before any result are recorded in each specification. Failed tests are reported as failed (Section 4).

## 4. Validation

**[TODO]** One combined table. The key points:
- **4c-1:** chain length distributions follow the analytic equilibrium of reversible polymerisation; wet–dry cycling lengthens chains.
- **4c-2:** the RAF threshold of Hordijk, Kauffman & Steel (2011) is reproduced (P = 0.5 at 1.05 and 1.27 catalysed reactions per molecule, against 1.06 and 1.25 in stage 1e). No networks appear without catalysis.
- **4c-3:**
  - copy fidelity is as specified;
  - Eigen's error threshold emerges (heredity lost at 0.5–1 errors per copy);
  - short parasites take over;
  - one pre-registered test failed at one of two points (test 9 at #1558); a finer grid, fixed before running, passed.
- **4c-4:** all nine tests passed. Without amphiphiles the model is bit-identical to 4c-3.
- **4c-5:**
  - Six of eight pre-registered tests passed.
  - **Test 3 failed** on its effect-size threshold: coupling was detected in 16 of 16 runs, but with a median ρ of 0.27 against a guessed 0.3.
  - **Test 6 failed**, the L4 positive control (3 of 16): its observation window opened after selection had finished. A revised test 6b, fixed before it ran, passed (16 of 16, with 1 of 16 false positives).
  - The point study of Section 5.5 ran after these failures as a recorded deviation from the stop rule.

## 5. Results

### 5.1 The ladder narrows at every step

| Level | Where it occurs | Source |
|---|---|---|
| E1, activated feedstock | 44–45% of environments | 4c-1, 4c-2, 4c-3, 4c-4 full maps |
| E2, persistent chains | 9.5–16% | 4c-1 (16%), later layers 8.5–10.7% |
| L1, catalytic network | 0 in the full space; 0.3% of a focused region with highly catalytic short chains | 4c-2 |
| L2_long, heredity | 0 in the full space; 5 confirmed environments from focused maps, plus 2 from a repetition | 4c-3 |
| E3, filled vesicles | 0.2% of the full space | 4c-4 |
| Heredity inside vesicles | 0 in all ground truth; 4 confirmed environments found by COMPASS | 4c-4 |
| L3, coupled protocell | at all 10 known coexistence environments, 298 of 435 coexistence runs | 4c-5 |
| L4, Darwinian evolution | none | 4c-5 |

Two observations break the expected order:
- **Catalytic networks do not need long chains.** 69% of L1 runs have no persistent chains of 8 or more, and the networks are made of short oligomers.
- **Vesicles alone do not bring heredity.** In the full map, the rare filled vesicles carry copies of only about 3 letters.

**On rarity and the search.** At equal numbers of simulated runs, COMPASS found heredity with long copies 4 times and the joint level once, where the ground truth found neither (one-sided p ≈ 0.06). This is suggestive but not significant. Where a level is common enough to compare per run (filled vesicles in a focused map), the rates are equal. The engine's proven advantage remains map accuracy per simulation, shown in Part I and in 4a and 4c-1 (about 2× for E1).

### 5.2 The environments at the top

| Point | Temperature (cycling) | Cycle period | pH | Salinity | Mg²⁺ | Residence time | Heredity and vesicles / 16 |
|---|---|---|---|---|---|---|---|
| engine12538 | 4.5 °C | 63 h | 7.58 | 0.74 M | 1.0 mM | 11,178 h | 16 |
| engine11230 | 7.1 °C (0–17 °C) | 484 h | 5.85 | 0.76 M | 38.7 mM | 63,326 h | 16 |
| engine15233 | 25.2 °C (19–31 °C) | 242 h | 4.39 | 0.83 M | 0.35 mM | 29,686 h | 16 |
| engine15539 | 0.0 °C (0–10 °C) | 992 h | 5.84 | 0.39 M | 0.55 mM | 20,978 h | 13 |

Common features at all four points:
- the cool phase lies below and the hot phase above T_m8, so duplexes pair, copy and separate once per cycle;
- molecules live far longer than a cycle;
- copying makes fewer than about 0.3 errors per 60-letter copy;
- inherited lineages are 86–100% enclosed in vesicles.

The 4c-3 heredity environments without membranes are more varied (8–39 °C, pH 2.6–8.0, cycles from 1 h to 26 days). **The hardest step is having heredity and stable membranes in the same place.** Five of the seven confirmed heredity environments are acidic (pH 2.6–4.7): at pH below 4 no membrane in the model can form, and at the other two only 0.4–1.4% of the membrane-constant space allows one.

### 5.3 The environmental sequence

We pooled all map runs and read off the conditions under which each level is reached ([4d-seq](../results/stage4dseq_sequence.md)).
- **No conflicts.** For every pair of consecutive levels and every variable, the 10–90% ranges overlap (114 comparisons).
- **The largest step is from catalysis to heredity.** It needs cold water (median 59 → 9 °C), deep cycling, strong thermal gradients (median factor 3.6×10³ → 2.5×10⁶), little UV and residence times of years.
- **Mixing between sites reverses.** It is high for vesicles alone but very low for heredity inside vesicles.

So the model needs no change of environment between levels, though it does not rule one out. Whether products survive such a change was not simulated. For the rare levels, the numbers describe the points found, not whole regions.

### 5.4 Real settings

[4e](../docs/stage4e_environment_classes.md).
- **Closest match: a cool, closed-basin lake or pond** with seasonal or multi-day temperature cycles and slow flushing. Carbonate-rich closed-basin lakes (Toner & Catling 2020) match residence time, near-neutral pH and cycling.
- **Volcanic hot-spring pools** (Damer & Deamer 2020) are too hot for copying in this model, and fresher than the points found.
- **Alkaline vents** (Kelley et al. 2005; Martin et al. 2008) do not match: no heredity point lies above pH 8, and vent flow implies short residence.
- **The cold preference agrees in direction with ice** as a medium for RNA replication (Attwater et al. 2010), but the model has no freezing.
- **Two conflicts:**
  - Salt. The top points are as salty as sea water, where measured fatty-acid vesicles fail (Milshteyn et al. 2018).
  - Acidity. Three of the four need a membrane with pKa_app near 6, below measured fatty-acid and alcohol mixtures, which form vesicles from about pH 6.5 (Jordan et al. 2019).

### 5.5 Coupled protocells, and the ceiling at L4

At the 10 known coexistence environments (16 seeds per point and membrane set, 608 runs), heredity and filled vesicles coexisted in 435 runs.

| Measure | Main runs | β = 0 control | Random-copy control |
|---|---|---|---|
| L3 loose (content helps growth) | 349 | 0 of 608 | 515 of 608 |
| L3 strict (and passed on by copying) | 298 | 0 | 0 |
| L4 | 0 | 0 | 0 |

- **The coupling is physical.** Filled vesicles grow faster because of osmotic pressure, and the coupling vanishes without it.
- **What makes a filled vesicle a coupled protocell is heredity.** With random-sequence copying, growth coupling remains, but nothing is passed on.
- **A typical coupled protocell** holds 3–12 chains of 64 letters, divides at 1.8–4.1 times the membrane of a newborn vesicle, and runs through about 290–300 vesicle generations.

**L4 does not occur.**
- A founder test shows that selection does happen: copying vesicle families leave 100% of the descendants against a drift null (16 of 16 runs).
- But the sweep is over within a few hundred hours. Copying then refills every surviving vesicle to the same load, so no heritable variation in content remains. When founders with 16 and with 1–4 templates were compared, their surviving families carried the same load per volume at the first checkpoint.
- At the real points, the correlation between an ancestor's load and its number of descendants is even slightly negative (−0.07 to −0.32).

## 6. What follows from the assumptions

| Conclusion | Verdict |
|---|---|
| Heredity needs temperature cycling across T_m8 | **partly built in**: strand separation by melting is a rule. Which cycles work is derived |
| Molecules must outlive a cycle | **partly built in**, a consequence of the copying rule |
| Error threshold near 0.5–1 errors per copy | **derived** |
| Short parasites win | **derived** |
| Acidic environments (pH < 4) cannot hold membranes | **built in** by the layer C ranges, which are already more permissive than measured fatty-acid mixtures |
| Heredity ends up inside vesicles | **partly built in** (the enclosure rules); that heredity survives being split over compartments is derived |
| The top environments are cool, slowly flushed, cycling pools | **derived**, within the layer C ranges |
| No conflicts along the ladder | **derived**, from the pooled maps |
| Content helps growth | **built in** (β); its presence at every coexistence point is a check |
| Contents are inherited by copying | **derived** |
| No sustained Darwinian evolution | **built in by omission**: sequences have no phenotype |
| Inherited sequences are about 60 letters | **built in** by the 64-letter cap |

## 7. Discussion

**[TODO]** Expand. Points to make:
- The model favours cool, enclosed, slowly flushed water with long temperature cycles. This lies between the hot-spring and closed-lake scenarios, and away from vents.
- **The coexistence of heredity and membranes** is a sharper bottleneck than catalysis or heredity alone. This echoes the experimental tension between Mg²⁺-dependent copying and Mg²⁺-sensitive membranes (Adamala & Szostak 2013), which our model only half represents.
- **The L4 ceiling** is a statement about what heredity alone can do. Copying keeps information, but without a phenotype that information gives no lasting advantage. This gives a testable prediction for both models and experiments: protocell populations evolve beyond the first sweep only if sequence affects copying, catalysis or the membrane.
- **COMPASS's role.** It made a 54-dimensional bottom-up map feasible on one GPU and found confirmed points for the rarest levels. Its advantage in finding rare points is suggestive but not proven.

## 8. Limits

- A 2-letter alphabet and a 64-letter cap. The cap sets how long inherited sequences are, not whether heredity happens.
- Coarse membrane laws (a stability window and an exponential critical concentration). The salt and acidity conflicts in Section 5.4 depend on them.
- Copying does not depend on Mg²⁺.
- Sequences have no phenotype beyond length and mass.
- The focused regions were chosen from earlier maps. No full joint map over all levels was run, and the extent of the L3 region is not mapped.
- No time-varying environment was simulated, and carry-over between levels is open.
- The layer C ranges decide much. Rare-level results are based on few points.

## 9. Predictions

1. Environments that can carry chemistry to coupled protocells share these features:
   - cool water (around 0–25 °C);
   - temperature cycles of days to weeks crossing the melting point of short duplexes;
   - residence times of years;
   - little UV and poor mixing between micro-sites.
2. Heredity and membranes compete for the same environmental variables (pH, salt, Mg²⁺). Experiments that combine nonenzymatic copying with fatty-acid vesicles should find their joint window narrower than either alone.
3. Without a sequence-dependent phenotype, compartmentalised replicators show one selective sweep and then stasis.

## References

**[TODO]** Merge from the stage reports; every entry must carry its check note.
- Adamala & Szostak (2013), *Science* 342, 1098–1100.
- Attwater et al. (2010), *Nat. Commun.* 1, 76.
- Chen, Roberts & Szostak (2004), *Science* 305, 1474–1476.
- Damer & Deamer (2020), *Astrobiology* 20, 429–452.
- Eigen (1971), *Naturwissenschaften* 58, 465–523 **[check pages]**.
- Hordijk & Steel (2004), *J. Theor. Biol.* 227, 451–461 **[check]**.
- Hordijk, Kauffman & Steel (2011), *Int. J. Mol. Sci.* 12, 3085–3101 **[check]**.
- Jordan et al. (2019), *Nat. Ecol. Evol.* 3, 1705–1714.
- Kelley et al. (2005), *Science* 307, 1428–1434.
- Li & Breaker (1999), *JACS* 121, 5364–5372 **[check]**.
- Mansy & Szostak (2008), *PNAS* 105, 13351–13355.
- Martin et al. (2008), *Nat. Rev. Microbiol.* 6, 805–814.
- Milshteyn et al. (2018), *Life* 8, 11.
- Radzicka & Wolfenden (1996), *JACS* 118, 6105–6109 **[check]**.
- Toner & Catling (2020), *PNAS*, doi:10.1073/pnas.1916109117 **[volume and pages]**.
- Zhu & Szostak (2009), *JACS* 131, 5705–5713.

## Use of AI

**[TODO]** State how AI assistants were used in the design, coding, analysis and writing, and that the author checked and is responsible for every claim.
