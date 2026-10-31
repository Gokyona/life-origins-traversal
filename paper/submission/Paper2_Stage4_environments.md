# How far up the ladder of life can environments carry a bottom-up chemistry? An ordered search over 13 environmental variables

**Gokyona Wahiro** (ORCID [0009-0008-1514-8856](https://orcid.org/0009-0008-1514-8856)) · Independent researcher, Toronto, Canada · 30 October 2026

Code, raw data and detailed stage reports: <https://github.com/Gokyona/life-origins-traversal> · Text: CC BY 4.0 · Code: MIT

Method paper: *COMPASS: mapping where life-like behaviour appears in stochastic models, with an ordered, non-repeating active search*, Zenodo preprint, doi:[10.5281/zenodo.23825371](https://doi.org/10.5281/zenodo.23825371)

---

## Abstract

Which environments can turn non-living chemistry into something that inherits, encloses itself and evolves? We built a fully bottom-up artificial chemistry with nothing living put in. It runs from feedstock through activated monomers, chains, catalysts and template copying to amphiphiles and vesicles. Its only inputs are 13 environmental variables, encoded as 19 dimensions with ranges taken from the literature. Rate constants that no measurement fixes are kept as explicit unknown dimensions, up to 54 dimensions in all. We searched this space layer by layer with COMPASS, an ordered, non-repeating active-learning engine, and detected seven levels automatically against paired controls.

The model climbs the ladder up to coupled protocells, and each step narrows the environments sharply:
- activated feedstock appears in about 45% of the space;
- persistent chains in about 8.5–16%;
- self-sustaining catalytic networks only when short chains are highly catalytic (0.3% of a focused region);
- heredity by template copying in a narrow class of environments (5 confirmed);
- heredity inside stable vesicles in 4 environments found by the search;
- coupled protocells wherever heredity and vesicles coexist (298 of 435 coexistence runs).

All environments at the top are cool to mild (0–25 °C), slowly flushed, enclosed water bodies with temperature cycles of days to weeks across the melting point of short duplexes, and near-neutral to slightly acidic water. The ladder shows no conflict between consecutive levels: one stable environment of this kind could carry the chemistry through all levels reached. The hardest step is having heredity and stable membranes in the same place. The ladder stops before Darwinian evolution. Copying families sweep once, after which copying erases all heritable variation in vesicle content, because sequences in the model have no phenotype. We therefore predict that sustained evolution in such systems needs a sequence-dependent trait that copying does not erase. For every result we state whether it follows from the model's assumptions or is derived. A supplementary study (*added 31 October 2026*) shows that heredity at the confirmed points depends in part on the order in which the copy step fills scarce chain slots, a numerical rule rather than chemistry; with fair allocation it survives at only some of them (Section 7).

## 1. Introduction

Most simulations of the origin of life choose one set of conditions and wait for something life-like to appear. A negative result then says little: the conditions may have been wrong, or the run too short. We ask the opposite question: **in which environments does a given chemistry climb the ladder of life, and where does it stall?** We answer it with a map rather than a single history.

Part I of this project ([stage 1–3 reports](https://github.com/Gokyona/life-origins-traversal/blob/main/paper/COMPASS_stage3.md)) built the mapping engine, COMPASS, and validated it on test models. There it needed 2–8 times fewer simulations than random or Sobol sampling for the same map accuracy. Part I had one fundamental limit: its conditions described the chemical system (cell size, error rate, parasite advantage), not the environment. This paper reports Part II. Its inputs are environmental variables only, and its chemistry is built from the bottom up, so that every level of life must emerge.

We use an extended ladder of seven levels, each scored automatically:
- **E1:** activated feedstock reaches a working concentration.
- **E2:** chains of at least 8 units persist.
- **L1:** a self-sustaining catalytic network (a persistent RAF set, Hordijk & Steel 2004) that raises its own products above a paired run without catalysis.
- **L2:** heredity. A sequence lineage outlives its molecules because it is copied, beating two paired random-copy controls, with copies of at least 6 letters (L2_long).
- **E3:** persistent vesicles that hold chains.
- **L3:** a coupled protocell. Vesicle growth depends on content, and the content is passed on by copying.
- **L4:** Darwinian evolution. A heritable vesicle trait rises beyond drift.

**Where life could start.** The main scenarios each rest on a physical argument.
- *Volcanic hot-spring pools* add wet–dry cycling, which drives polymerisation and encapsulates polymers in membranes (Damer & Deamer 2020). Fatty-acid vesicles form in water from such pools (Milshteyn et al. 2018).
- *Alkaline hydrothermal vents* supply H₂ and natural gradients (Martin et al. 2008). Lost City fluids are 40–90 °C at pH 9–11 (Kelley et al. 2005), and thermal gradients in pores can concentrate nucleotides by more than 10⁸ (Baaske et al. 2007).
- *Ice* concentrates RNA in the eutectic phase and protects it from hydrolysis, so that a polymerase ribozyme makes longer products (Attwater et al. 2010).
- *Closed-basin carbonate-rich lakes* can accumulate phosphate to high levels (Toner & Catling 2020).

Each scenario is usually argued for one or two steps at a time. Here we ask how the whole sequence of steps constrains the environment.

**Models of the steps.** Theory gives sharp results for single steps:
- the error threshold of copying (Eigen 1971);
- the level of catalysis at which autocatalytic sets appear (Hordijk & Steel 2004; Hordijk, Kauffman & Steel 2011);
- the trade-off between replicators inside vesicles (Fontanari, Santos & Szathmáry 2006).

These models are defined over system parameters, such as catalytic probability or error rate. Only some of these map to an environment. We keep the same mechanisms but make the environment the input, and treat every parameter that no measurement fixes as an explicit unknown dimension.

## 2. Model

The simulator was built in five layers. Each layer was validated before the next was added, and none changes the layers below. Full specifications: [4c-1](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/stage4c1_spec.md), [4c-2](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/stage4c2_spec.md), [4c-3](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/stage4c3_spec.md), [4c-4](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/stage4c4_spec.md), [4c-5](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/stage4c5_spec.md).

**Environment (layer A).** There are 13 variables, encoded as 19 dimensions:
- temperature, pressure, pH, H₂ availability and UV flux;
- salinity and Mg²⁺;
- feedstock supply and composition;
- residence time and mixing between sites;
- evaporation and thermal-gradient concentration;
- wet–dry and hot–cold cycling (depth and period);
- mineral surfaces (amount and type);
- available time.

Ranges and their sources are in [stage4a_ranges.md](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/stage4a_ranges.md). They run from acid springs to alkaline vents and from flowing channels to closed pools.

**Unknown constants (layer C).** Constants that no measurement fixes, such as condensation rate, catalytic probability, copy speed-up, error rate, duplex melting temperature and membrane stability limits, are dimensions of the map with literature-based ranges. Their number grows from 10 (4c-1) to 35 (4c-4), so the map grows from 29 to 54 dimensions.

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
- At 20–30 dimensions, the engine still needed about 3× fewer simulations than Sobol sampling for the same accuracy (2.8–4.0× depending on the type of added dimension; stage 4a).

**Maps per layer.**
- **Full map:** 1,000 random environments × 16 seeds as ground truth, and 32,000 COMPASS runs (4c-1: 2,000 × 32 ground truth and 64,000 runs each for COMPASS and Sobol sampling).
- **Focused maps:** regions chosen from the full map. Environments outside a region's constraints are masked, not simulated. Each has 500 × 16 ground truth and 16,000 COMPASS runs (the 4c-2 focused map: 1,000 × 16 and 32,000).

**Confirmation.** From 4c-3 on, every candidate for a rare level was re-run with 16 new seeds. The 4c-2 catalytic networks were not confirmed this way.

**Controls.**
- L1: a paired run without catalysis.
- L2: two paired runs where copying makes chains of the same length with random sequence.
- L3: paired runs with β = 0 (no osmotic coupling) and with random-sequence copying.

**Pre-registration.** Every layer had validation tests with criteria written before they ran. Design changes made before any result are recorded in each specification. Failed tests are reported as failed (Section 4).

## 4. Validation

Test criteria were written into the layer's specification before the tests ran; where a test design itself was faulty, the redesign is recorded in the validation report. Every layer passed a kill-and-resume test with bit-identical output (4c-2 through a real power cut), and layers 4c-3 to 4c-5 passed a memory pre-flight at worst-case settings.

| Layer | Tests | Passed | What the tests show | Failures and changes |
|---|---|---|---|---|
| 4c-1 polymers | 6 | 6 | Chain lengths match the exact steady state (total-variation distance ≤ 0.003); the hydrolysis law matches its source (slope within 0.2%); wet–dry cycling lengthens chains (mean 3.12 against 2.12 letters at 80 °C); washout or no feed removes chains | Test (a) criterion corrected before running, from the Flory distribution to the exact steady state |
| 4c-2 catalysis | 4 | 4 | The RAF threshold is reproduced: P = 0.5 at 1.05 and 1.27 catalysed reactions per molecule, against 1.06 and 1.25 in stage 1e (Hordijk, Kauffman & Steel 2011). No networks without catalysis (0 of 72) or with catalysis reshuffled each step | 4 changes to the specification after the first validation run, before any map |
| 4c-3 copying | 9, plus 9b | 8, plus 9b | Measured copy error 0.0494 against ε = 0.05; strand separation needs cycling across T_m8; short parasites take over (short share 0.994); null false positives 1 of 64 and 0 of 64 | Test 4 passed on its main criteria but missed a secondary one (0.083 against 0.073), after two earlier designs failed for reasons of test design. Tests 2 and 3 ran on earlier code. **Test 9 failed** at #1558 (error threshold): heredity was already lost at the first grid point. Test 9b, with a finer grid and its criterion fixed before running, passed at both points; heredity is lost near ε·L_copy ≈ 1 (Eigen 1971) |
| 4c-4 vesicles | 9 | 9 | Without amphiphiles bit-identical to 4c-3; vesicles only inside the stability window; binomial encapsulation and division; filled vesicles gain membrane only with β > 0 | 3 rules changed after the memory pre-flight, before validation |
| 4c-5 coupling | 8, plus 6b | 6, plus 6b | Bit-identical to 4c-4; records consistent over about 200,000 vesicle ids; inheritance 16 of 16 with copying and 0 of 16 without; 0 of 80 false positives at real points | The inheritance criterion was revised three times and one fixture fixed, all before any result. **Test 3 failed** on its effect size: coupling was detected in 16 of 16 runs (p = 0.001 each), but ρ reached the guessed 0.3 in only 3 of 16 (median 0.27). **Test 6 failed**, the L4 positive control (3 of 16): its window opened after selection had finished. Test 6b, fixed before it ran, passed (16 of 16; 1 of 16 false positives) |

After the two 4c-5 failures, the stop rule forbade maps. The point study of Section 5.5 was run anyway as a recorded deviation; no 4c-5 map was run.

## 5. Results

### 5.1 The ladder narrows at every step

| Level | Where it occurs | Source |
|---|---|---|
| E1, activated feedstock | 44–47% of environments | 4c-1, 4c-2, 4c-3, 4c-4 full maps |
| E2, persistent chains | 8.5–16% | 4c-1 (16%), later layers 8.5–10.7% |
| L1, catalytic network | 0 of 16,000 ground-truth runs in the full space (39 COMPASS runs near the boundary); 0.3% of a focused region with highly catalytic short chains | 4c-2 |
| L2_long, heredity | 0 in the full space; 5 confirmed environments from a targeted scan and a focused map, plus 2 from a repetition | 4c-3 |
| E3, filled vesicles | 0.2% of the full space | 4c-4 |
| Heredity inside vesicles | 0 in all ground truth; 4 confirmed environments found by COMPASS | 4c-4 |
| L3, coupled protocell | at all 10 known coexistence environments, 298 of 435 coexistence runs | 4c-5 |
| L4, Darwinian evolution | none | 4c-5 |

Two observations break the expected order:
- **Catalytic networks do not need long chains.** 69% of L1 runs have no persistent chains of 8 or more, and the networks are made of short oligomers.
- **Vesicles alone do not bring heredity.** In the full map, the rare filled vesicles carry copies of only about 3 letters.

**On rarity and the search.** At equal numbers of simulated runs, COMPASS found heredity with long copies 4 times and the joint level once, where the ground truth found neither (one-sided p ≈ 0.06). This is suggestive but not significant. Where a level is common enough to compare per run (filled vesicles in a focused map), the rates are equal. The engine's proven advantage remains map accuracy per simulation, shown in Part I, 4a and 4c-1. In 4c-1, the only layer with a matched Sobol baseline, COMPASS reaches Sobol's final accuracy with about 2× (E1) and 1.6× (E2) fewer runs, and 95% of the achievable gain in accuracy over a trivial classifier with about 2.4× and 2.9× fewer runs (analysis of the stored accuracy curves, to be added to the repository).

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
- copying is accurate (about 0.08 errors per copy at engine12538; the 4c-3 heredity points need fewer than about 0.3);
- inherited lineages are 86–100% enclosed in vesicles.

The 4c-3 heredity environments without membranes are more varied (6–39 °C, pH 2.6–8.0, cycles from 1 h to 26 days). **The hardest step is having heredity and stable membranes in the same place.** Five of the seven confirmed heredity environments are acidic (pH 2.6–4.7): at pH below 4 no membrane in the model can form, and at the other two only 0.4–1.4% of the membrane-constant space allows one.

### 5.3 The environmental sequence

We pooled all map runs and read off the conditions under which each level is reached ([4d-seq](https://github.com/Gokyona/life-origins-traversal/blob/main/results/stage4dseq_sequence.md)).
- **No conflicts.** For every pair of consecutive levels and every variable, the 10–90% ranges overlap (114 comparisons).
- **The largest step is from catalysis to heredity.** It needs cold water (median 59 → 9 °C), deep cycling, strong thermal gradients (median factor 3.6×10³ → 2.5×10⁶), little UV and residence times of years.
- **Mixing between sites reverses.** It is high for vesicles alone but very low for heredity inside vesicles.

So the model needs no change of environment between levels, though it does not rule one out. Whether products survive such a change was not simulated. For the rare levels, the numbers describe the points found, not whole regions.

### 5.4 Real settings

[4e](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/stage4e_environment_classes.md).
- **Closest match: a cool, closed-basin lake or pond** with seasonal or multi-day temperature cycles and slow flushing. Carbonate-rich closed-basin lakes (Toner & Catling 2020) match residence time, near-neutral to slightly acidic pH and cycling.
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
| Heredity arises at all | **built in** by an extension rate about 10⁹–10¹¹ times faster than measured chemistry (Section 7) |
| Heredity (L2) at the confirmed points | **partly a numerical artefact**: it depends on the order in which the copy step fills scarce chain slots; with fair allocation it is lost at 2–6 of the 7 points (Section 7, *added 31 October 2026*) |
| Heredity points pair weakly yet copy accurately | **built in** by independent layer C ranges for pairing strength, accuracy and speed; no measured base pair combines these values (Section 7) |
| Inherited sequences are about 60 letters | **built in** by the 64-letter cap |

## 7. Discussion

**Which scenario the model favours.** The environments at the top are cool, enclosed, slowly flushed water bodies with temperature cycles of days to weeks. Among the scenarios in Section 1, this is closest to closed-basin lakes (Toner & Catling 2020). It shares cycling with the hot-spring scenario (Damer & Deamer 2020), but in this model hot water above the duplex melting point blocks copying. It points away from vents (Martin et al. 2008): no heredity point lies above pH 8, and through-flow gives residence times far shorter than the years the top points need. The cold preference agrees in direction with ice as a medium for replication (Attwater et al. 2010), by a different mechanism: in the model, cold water lengthens molecular lifetimes, while in ice it is concentration in the eutectic phase. It also agrees in direction with nonenzymatic copying driven by in-place activation, which works far better at −3 °C than at room temperature (Aitken et al. 2023). There, too, the mechanism differs: with pre-activated building blocks, cold slows copying (42% of primers extended at room temperature against 8% at +3 °C), so the benefit comes from the activation chemistry, which loses less reagent to side reactions in the cold, not from longer-lived molecules or duplexes (*revised 30 October 2026*; [part3_cold_gap.md](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/part3_cold_gap.md)).

These are statements about the model's chemistry, not about the history of the Earth. Where the model is too coarse to decide, we say so. The two clearest cases are the salt conflict and the acid-tolerant membranes needed at three of the four top points (Section 5.4).

**The bottleneck is coexistence.** No single step is the hardest. Catalysis needs highly catalytic short chains, and heredity needs cycling across T_m8 and long lifetimes, but each happens in its own region. The sharpest cut comes when heredity and membranes must occur in the same place. Most heredity environments are acidic, where no membrane in the model forms. This echoes an experimental tension: nonenzymatic copying needs Mg²⁺, which breaks fatty-acid membranes unless citrate chelates it (Adamala & Szostak 2013). Our model represents only the membrane side of that tension, because copying does not depend on Mg²⁺. Adding that dependence would most likely narrow the joint window further.

**The L4 ceiling.** The model shows what heredity alone can do inside compartments. Copying families win one fast sweep, as competition between protocells predicts (Chen, Roberts & Szostak 2004). After the sweep, copying refills every vesicle to the same load, so no heritable variation in content remains. Information is kept, but it gives no lasting advantage, because sequence affects nothing but length. This failure is informative. It turns into a prediction that models and experiments can test: protocell populations evolve beyond the first sweep only if sequence affects copying, catalysis or the membrane.

A later check, made while building an experimental platform on the same engine, shows a second, technical reason. A template can start at most one copy per time step, with probability 1 − exp(−λ·dt). At the ten coexistence environments, re-run with the same engine at the paper's step size, 99.1–100% of template-steps have λ·dt > 3, with a median λ·dt of at least 10⁶. Every template therefore copies once per step whatever its rate: copy throughput is set by the step size and the number of templates, not by any property of the template. A trait that changes the copy rate could therefore not be selected in this discretisation. The L4 ceiling is thus a joint result of the missing phenotype and of this cap on copy starts. Neither depends on the environment, but the second is a property of the numerical scheme, not of chemistry. The same check shows how lineages end: 93–99.9% of losses of inherited chains happen when their vesicle washes out, which is independent of content, and only 0.06–7% through bond hydrolysis. A trait acting on chain stability would therefore also have little to select on at these points.

**Copying is far faster than measured chemistry.** A further check, made on the same engine after the maps, puts all heredity results here on a conditional footing. The model's template-directed extension rate is k_on·[M*], where k_on multiplies the uncatalysed condensation constant by the template enhancement k_t and by the concentration factors of layer A (evaporation, thermal-gradient accumulation and the dry phase). At the coexistence point engine12538 this gives about 1.7 × 10⁶ per molecule per hour, or about 5 × 10¹⁰ letters per hour at the monomer levels reached. In experiments, nonenzymatic copying with activated nucleotides at 20–50 mM and 100–200 mM Mg²⁺ takes from 30 minutes to about a day per added letter, even with activated helper oligonucleotides (Prywes et al. 2016). The model's copying is therefore about 10⁹–10¹¹ times faster than measured chemistry, and a copy finishes long before its duplex melts.

When we set extension to a rate anchored in that literature (about one letter per hour) and apply the concentration factors only to the monomer concentration, template copying stops entirely. Short duplexes dissociate about 10⁶ times faster than a letter is added. Two further routes also failed to start heredity from scratch within the layer C ranges: templated ligation of chance oligomers (a piece must be 17–39 letters long to stay bound, and a random pool holds almost no such complements), and chance long primers at the most stable pairing allowed. A semi-analytic estimate from the engine's own formulas suggests what would be needed: 8-mer duplexes melting above about 50 °C, untemplated growth reduced 7–27-fold, and a supply of 6–8-letter oligomers of 1–10% of the monomer flux. This is the known problem of nonenzymatic replication, reproduced by the model. The heredity levels reported here (L2, E3, L3) should therefore be read as conditional on fast copying chemistry, an assumption the model did not make explicitly.

The same check found two smaller numerical effects. Untemplated growth is capped at 4 letters per chain per step, and this cap binds in about 90% of growable chain-steps at two of the three coexistence points checked. That leaves more activated monomer for copying than mass action would. However, 61–82% of chains already sit at the 64-letter cap, so the level detectors are probably affected only modestly. Also, the dissociation rate in the cycle-averaged step is an arithmetic mean over the wet and dry phases. This favours melting, whereas real copying would happen in the cold phase.

**No real base pair fits the heredity constants** (*added 30 October 2026*). The two letters of the model do not stand for a named base pair: duplex stability (ΔH_bp, T_m8), copying accuracy (ε) and copying speed are separate layer C constants. A later check translated them into real chemistry with nearest-neighbour parameters for RNA (Turner 2004 set, Xia et al. 1998). An 8-letter duplex of G·C pairs would have a model melting point T_m8 of about 78–92 °C, and one of A·U pairs about 2–25 °C. The confirmed heredity points have T_m8 of 11.9–30.1 °C, so they pair like A·U. Yet they also need fewer than about 0.0013–0.0074 errors per letter and fast copying, and measured nonenzymatic copying is accurate and fast only for G·C-rich templates; A and U are copied at least 1,000 times more slowly without helper oligonucleotides (Prywes et al. 2016). A pair of intermediate strength, such as A with 2-thiouridine, copies fast (Fang et al. 2025) but is about 70% as stable as G·C, which puts T_m8 near 65–90 °C. The stable end of the maps does not offer a way out: in the Stage 4 runs with T_m8 of 70–90 °C, none of about 30,700 runs reached L2_long, where about 9 were expected at the rate seen at 10–30 °C. The cause is bond hydrolysis during the hot phase needed to melt such duplexes. The strongest pairing enthalpies allowed (ΔH_bp of 53–56) were never sampled. This gap is independent of the copying-speed gap above: even at the model's copying speed, the heredity points combine constants that no measured base pair has. Both limits must be removed before the heredity results can describe real chemistry; this is the starting point of follow-up work ([part3_alphabet_transfer.md](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/part3_alphabet_transfer.md)).

**Heredity depends on how scarce chain slots are shared out** (*added 31 October 2026*). A supplementary study re-ran the 7 confirmed L2 points (16 seeds each, with both random-copy controls) under changed rules ([Stage 4 supplementary study](https://github.com/Gokyona/life-origins-traversal/blob/main/stage4/supplement_robustness)). The unchanged model reproduces the published rates (Fisher exact p ≥ 0.46 at every point). Several simplifications do not matter for L2:
- capping untemplated growth at 16 instead of 4 letters per step;
- washing out chains of 4 or more letters at half the rate;
- an exact Poisson sampler;
- resolving duplex melting by phase instead of averaging it over the cycle (tested at the 3 long-period points).

In every one of these variants, L2 holds at every point tested. Slowing hydrolysis of paired strands changes nothing, because no paired strand breaks in any run.

What does matter is the copy step. At these points copying is limited by free chain slots: 88–99% of copy attempts find no free slot. The engine fills the free slots in slot-index order, so the same long-lived templates are served step after step, and the top 10% of template lineages receive 81–91% of all copies. With the same number of copies but the slots filled in random order, L2 is lost at 2 of the 7 points (15 → 1 and 9 → 0 of 16 runs) and its share of chains halves at the other 5. With several copy starts per template and step, shared out in proportion to each template's copying rate, L2 survives at none of the 7 points when the wait for duplex melting is enforced, and fully at only 1 point (partly at 1–2 others) without that wait, whatever the slot order. Heredity at the confirmed points is therefore produced largely by a priority rule of the numerical scheme in a slot-limited regime, not by copying rates. Whether the same holds for the vesicle and protocell levels (E3, joint, L3) was not tested.

The same study found that the GPU's built-in Poisson sampler (torch.poisson on the Intel XPU), used for copy initiation and similar counts in all Stage 4 runs, is biased at small means (−1.1% at a mean of 0.02, with excess spread between seeds). Replacing it with an exact sampler leaves L2 unchanged at the 7 points; its effect at the coexistence points was not tested.

**The role of the search.** COMPASS made a map over up to 54 dimensions feasible on one consumer GPU. It found, and confirmed with new seeds, specific environments for the rarest levels. At equal budget it found long-copy heredity 4 times and the joint level once, where the ground truth found neither (p ≈ 0.06). This advantage in finding rare points is suggestive but not proven. Its proven advantage is map accuracy per simulation, about 2.4–2.9× fewer runs to reach 95% of the achievable accuracy gain in the one layer with a matched baseline. Much of its apparent advantage in raw hit counts comes from the ground truth's design of 16 seeds per point, not from targeting. The method is the subject of a separate paper.

## 8. Limits

- A 2-letter alphabet and a 64-letter cap. The cap sets how long inherited sequences are, not whether heredity happens.
- Coarse membrane laws (a stability window and an exponential critical concentration). The salt and acidity conflicts in Section 5.4 depend on them.
- Copying does not depend on Mg²⁺.
- Sequences have no phenotype beyond length and mass.
- Copy starts are capped at one per template per time step, and this cap binds at the top environments (Section 7). Results that depend on copy rates, including the L4 ceiling, may change with a finer step or with several starts per step. A later test confirms this for heredity: with several starts per step, shared out by rate, L2 is lost at most of the 7 confirmed points (Section 7, *added 31 October 2026*).
- Copying is limited by free chain slots at the heredity points, and the slots are filled in a fixed index order that favours the same templates. Filling them in random order loses L2 at 2 of 7 points and halves it at the others, so heredity depends in part on this numerical rule (Section 7, *added 31 October 2026*). The vesicle and protocell levels were not re-tested.
- Template-directed extension is about 10⁹–10¹¹ times faster than measured nonenzymatic copying (Section 7). At a literature rate, no heredity arises from scratch in this chemistry, so all heredity levels are conditional on fast copying.
- The two letters are not a named base pair. The heredity points combine A·U-like pairing strength with G·C-like copying accuracy and speed, which no measured pair offers; at G·C-like stability, no long-copy heredity was found in the maps (Section 7). This limit is independent of the copying-speed limit above.
- Untemplated growth is capped at 4 letters per chain per step, and the cycle-averaged step averages duplex dissociation over the wet and dry phases (Section 7).
- Vesicle counts include empty membrane remnants: about 1–2% of vesicles at the coexistence points. They hold no chains, so no level detector is affected.
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

Each entry was checked against the original or a publisher, library or author record; the date and source of each check are in the stage reports and the [condition review](https://github.com/Gokyona/life-origins-traversal/blob/main/docs/condition_review.md).

- Adamala, K. & Szostak, J. W. (2013). Nonenzymatic template-directed RNA synthesis inside model protocells. *Science* 342(6162), 1098–1100.
- Aitken, H. R. M., Wright, T. H., Radakovic, A. & Szostak, J. W. (2023). Small-molecule organocatalysis facilitates in situ nucleotide activation and RNA copying. *Journal of the American Chemical Society* 145(29), 16142–16149. doi:10.1021/jacs.3c04635. Read in full on 30 October 2026; title and authors checked on Crossref.
- Attwater, J., Wochner, A., Pinheiro, V. B., Coulson, A. & Holliger, P. (2010). Ice as a protocellular medium for RNA replication. *Nature Communications* 1, 76. doi:10.1038/ncomms1076.
- Baaske, P., Weinert, F. M., Duhr, S., Lemke, K. H., Russell, M. J. & Braun, D. (2007). Extreme accumulation of nucleotides in simulated hydrothermal pore systems. *PNAS* 104(22), 9346–9351.
- Chen, I. A., Roberts, R. W. & Szostak, J. W. (2004). The emergence of competition between model protocells. *Science* 305(5689), 1474–1476.
- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452.
- Eigen, M. (1971). Selforganization of matter and the evolution of biological macromolecules. *Die Naturwissenschaften* 58(10), 465–522.
- Fang, Z., Jia, X., Xing, Y. & Szostak, J. W. (2025). Nonenzymatic RNA copying with a potentially primordial genetic alphabet. *PNAS* 122(21), e2505720122. Read in full on 30 October 2026; title and authors checked on Crossref.
- Fontanari, J. F., Santos, M. & Szathmáry, E. (2006). Coexistence and error propagation in pre-biotic vesicle models: a group selection approach. *Journal of Theoretical Biology* 239, 247–256.
- Hordijk, W. & Steel, M. (2004). Detecting autocatalytic, self-sustaining sets in chemical reaction systems. *Journal of Theoretical Biology* 227(4), 451–461.
- Hordijk, W., Kauffman, S. A. & Steel, M. (2011). Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. *International Journal of Molecular Sciences* 12(5), 3085–3101.
- Jordan, S. F., Rammu, H., Zheludev, I. N., Hartley, A. M., Maréchal, A. & Lane, N. (2019). Promotion of protocell self-assembly from mixed amphiphiles at the origin of life. *Nature Ecology & Evolution* 3, 1705–1714. doi:10.1038/s41559-019-1015-y.
- Kelley, D. S. et al. (2005). A serpentinite-hosted ecosystem: the Lost City hydrothermal field. *Science* 307(5714), 1428–1434.
- Li, Y. & Breaker, R. R. (1999). Kinetics of RNA degradation by specific base catalysis of transesterification involving the 2′-hydroxyl group. *Journal of the American Chemical Society* 121(23), 5364–5372.
- Mansy, S. S. & Szostak, J. W. (2008). Thermostability of model protocell membranes. *PNAS* 105(36), 13351–13355.
- Martin, W., Baross, J., Kelley, D. & Russell, M. J. (2008). Hydrothermal vents and the origin of life. *Nature Reviews Microbiology* 6, 805–814.
- Milshteyn, D., Damer, B., Havig, J. & Deamer, D. (2018). Amphiphilic compounds assemble into membranous vesicles in hydrothermal hot spring water but not in seawater. *Life* 8(2), 11. doi:10.3390/life8020011.
- Prywes, N., Blain, J. C., Del Frate, F. & Szostak, J. W. (2016). Nonenzymatic copying of RNA templates containing all four letters is catalyzed by activated oligonucleotides. *eLife* 5, e17756. doi:10.7554/eLife.17756. Checked on 30 October 2026 on the publisher's page.
- Radzicka, A. & Wolfenden, R. (1996). Rates of uncatalyzed peptide bond hydrolysis in neutral solution and the transition state affinities of proteases. *Journal of the American Chemical Society* 118, 6105–6109. Values as tabulated in Wolfenden & Snider (2001), *Accounts of Chemical Research* 34(12), 938–945.
- Toner, J. D. & Catling, D. C. (2020). A carbonate-rich lake solution to the phosphate problem of the origin of life. *PNAS* 117(2), 883–888. doi:10.1073/pnas.1916109117. Volume, issue and pages checked on 29 October 2026 in the Virtual Planetary Laboratory publication record.
- Xia, T. et al. (1998). Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson–Crick base pairs. *Biochemistry* 37(42), 14719–14735. Values as tabulated in the Nearest Neighbor Database (Turner 2004 set), checked on 30 October 2026.
- Zhu, T. F. & Szostak, J. W. (2009). Coupled growth and division of model protocell membranes. *Journal of the American Chemical Society* 131(15), 5705–5713.

## Use of AI

The author used an AI assistant, Claude (Anthropic), throughout this project:
- **Design:** it proposed and drafted the model layers, the detection criteria and the pre-registered validation tests, which the author reviewed and approved.
- **Code and runs:** it wrote all simulation and analysis code and ran it on the author's computer, under rules set by the author.
- **Literature:** it found and checked the sources for the variable ranges and the comparisons.
- **Writing:** it drafted the stage reports and this paper.

The author chose the research question, the scope and the environmental framing. The author decided the main choices: a fully bottom-up chemistry, a 2-letter alphabet, which steps to run and when to stop. The author reviewed the results and is responsible for every claim. All code, specifications, results and failed tests are public in the repository.
