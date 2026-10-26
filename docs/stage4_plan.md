# Stage 4 plan: which environments give rise to life

**Status:** draft for review, 25 October 2026 · **Part II, first stage** · **Builds on:** [condition review](condition_review.md), [stage 4 sketch](stage4_sketch.md), [stage 3 paper](../paper/COMPASS_stage3.md)

## 1. Goal

Part I built and validated COMPASS on models whose conditions describe the chemical system (cell size, error rate, parasite advantage). Stage 4 asks the project's real question:

> **Given a local environment, how likely is it that chemistry there reaches each level of life?**

The only inputs are the 13 approved environmental variables (layer A). Everything else is computed from them, or is an explicitly uncertain constant that is shown on the map.

## 2. Architecture

```
layer A: 13 environmental variables (inputs, set independently)
   │   sub-models: known laws with measured parameters, each checked against experiments
   ▼
layer B: derived chemistry (computed, never set)
   concentration factor · polymer length distribution · polymerisation/hydrolysis balance ·
   copying error per position · catalytic probability · vesicle formation and stability ·
   membrane permeability · degradation rates · effective food supply
   │
   ▼
life model: gates E1–E3, then the stage 3 protocell simulation for L1–L4
   │
   ▼
COMPASS maps over layer A (+ layer C uncertain constants)
```

**Layer C: uncertain constants.** Where a sub-model needs a constant that has no measured value, or a published range, the constant becomes an extra map dimension over that range, instead of being fixed by us. Importance analysis then shows whether the conclusion depends on it. This is the rule that keeps the maps honest.

### 2.1 Extended life ladder

The stage 3 model starts with vesicles that already contain replicators. Stage 4 adds three earlier gates, so the map also says where chemistry stalls before cells exist:

| level | meaning | decided by |
|---|---|---|
| E1 building blocks | monomers (nucleotides or amino acids, and amphiphiles) are made or supplied, and concentrated above a working threshold | feedstock, UV, redox, concentration mechanisms, dilution, degradation |
| E2 polymers | chains long enough to carry information or catalyse (length ≥ L_min) form faster than they hydrolyse | polymerisation/hydrolysis balance, temperature, water activity, cycling, minerals |
| E3 compartments | vesicles form and stay intact | amphiphile concentration, pH, salinity, Mg²⁺, temperature |
| L1–L4 | as in Part I: self-maintenance, heredity, coupled protocell, Darwinian evolution | the stage 3 protocell simulation, with its parameters computed from layer B |

The levels stay nested: L1 requires E1–E3. The joint nested map of stage 3e handles all seven levels.

### 2.2 How the protocell parameters are derived

| Part I parameter | derived from (layer B) |
|---|---|
| copying error μ | error per position (temperature, Mg²⁺, pH) × template length (E2) |
| food | effective monomer supply after concentration and dilution |
| degradation d | hydrolysis and UV damage rates (temperature, pH, UV), plus washout (flow) |
| leak lk | membrane permeability (membrane composition, temperature, ions) |
| membrane synthesis km, division size S | amphiphile supply and vesicle size (E3 sub-model) |
| catalytic support of metabolism | catalytic probability → autocatalytic-set probability (the stage 1e RAF model) |
| heritable variation η | error rate on the functional part of the template |
| spatial structure | mixing and geometry (well-mixed vs grid; exchange between cells) |
| cycles | periodic modulation of rates, as in stage 1b |
| available time | run length in generations |

Parasite advantage σ and the number of gene types K are not environmental. They stay as layer C dimensions with their Part I ranges.

## 3. Ranges of the layer A variables (provisional)

The ranges cover the stated scope: liquid water, carbon chemistry, the surface or sea floor of an Earth-like planet. They will be fixed from published sources in step 4a, and every value will be checked against the original source.

| variable | provisional range or encoding |
|---|---|
| temperature | 0–150 °C (above 100 °C only where pressure keeps water liquid) |
| pressure | 1–1,000 bar |
| pH | 2–12 |
| redox | reducing to oxidising, as H₂ availability (log scale) |
| UV | none to full unshielded surface flux (log scale) |
| ion composition | salinity 0 to above sea water; Mg²⁺ and Ca²⁺ concentrations (log scale) |
| feedstock composition and supply | supply rate (log scale) plus composition fractions: nucleotide precursors, amino-acid precursors, amphiphiles |
| dilution and flow | residence time (log scale) |
| mixing and geometry | from well-mixed to pore-confined (Péclet number) |
| concentration mechanisms | evaporation factor and thermal-gradient factor (log scale) |
| cycling | amplitude and period of wet–dry and hot–cold cycles |
| mineral surfaces | catalytic surface availability (log scale) and type (none, clay, metal sulphide) |
| available time | 1 to 10⁶ years, translated into generations |

Some variables need more than one number (cycling, composition, ions), so layer A has about 18–20 dimensions in practice, plus about 5–10 layer C constants. That is about twice the largest dimension COMPASS has handled so far (12), hence the scaling test in 4a.

## 4. Steps

**4a. Framework and engine scaling.**
- Fix the ranges and encodings of layer A from sources.
- Define the sub-model interfaces and the E1–E3 gates.
- Test COMPASS at 20–30 dimensions: pad the stage 3 models with distractors and correlated dimensions, and measure the accuracy and saving against Sobol sampling. If the saving collapses, a screening step (importance on a cheap first map, then refinement in the important dimensions) is added before 4d.

**4b. Sub-models, each checked against experiments before use.**
The order goes from clearest experimental data to least:
1. Vesicle formation and stability (E3): the window in pH, salinity, divalent ions and temperature. Target: high Mg²⁺ breaks fatty-acid vesicles, and citrate protects them while still allowing RNA copying (Adamala & Szostak 2013).
2. Concentration: thermal gradients in pores concentrate nucleotides more than 10⁸-fold (Baaske et al. 2007); evaporation in wet–dry cycles.
3. Polymerisation/hydrolysis balance (E2): the length distribution under constant water and under wet–dry cycling. Target: cycles drive polymerisation and encapsulate polymers (Damer & Deamer 2020).
4. Copying error and rate per position: dependence on temperature and Mg²⁺.
5. Building-block supply (E1): UV-driven precursor synthesis against UV damage (Patel et al. 2015); the redox-driven supply of the vent scenario (Martin et al. 2008).
6. Catalytic probability: chain length and mineral surfaces → reactions catalysed per molecule → the RAF probability of stage 1e (Hordijk, Kauffman & Steel 2011).
Each sub-model gets a short report: the law used, the parameter sources, the experiment reproduced, and what is layer C. A sub-model that cannot reproduce its target is not used.

**4c. Coupled model.**
Chain layer A → layer B → gates → protocell simulation. Validate it on hand-picked environments: does it behave sensibly in clear cases (for example boiling, strongly oxidising water with no feedstock gives nothing)? Measure the cost per simulation.

**4d. Maps.**
- COMPASS maps for E1–E3 and L1–L4 over layer A + layer C, with the joint nested map.
- Importance of each environmental variable.
- Minimal sufficient environments: boxes that give each level with 95% probability.
- Environmental knockouts: no cycling, no minerals, no concentration mechanism.

**4e. Real environment classes.**
Report the maps for real settings, each defined as a region of layer A from the literature: deep-sea alkaline vents, land hot springs with wet–dry pools, shallow lagoons and tidal pools, cold settings (ice, eutectic brines), and mineral pores. The question per class: which levels are reachable, and which variable limits it?

**4f. Robustness and paper.**
- Vary detector thresholds (three settings per level).
- Check how conclusions depend on layer C.
- List testable predictions: environments the map calls favourable or impossible that have not been tested in a laboratory.
- Stage 4 paper; homepage and log.

## 5. Milestones and acceptance

| step | output | accept when |
|---|---|---|
| 4a | ranges with sources; interfaces; scaling test | COMPASS still beats Sobol at 20+ dimensions, or a screening step that restores it |
| 4b | six sub-models with reports | each reproduces its experimental target within stated error |
| 4c | coupled model | clear cases behave correctly; cost per simulation known |
| 4d | maps E1–L4 | accuracy on held-out ground truth reported; importance stable across repetitions |
| 4e | per-class report | each class traced to literature ranges |
| 4f | paper | references verified; predictions stated |

## 6. Risks

- **Too many unknown constants.** If layer C dominates importance, the honest result is "the answer depends on X, which has not been measured". That is still a useful result and a pointer for experiments.
- **Dimension.** See 4a.
- **The protocell core remains a model assumption.** Stage 4 computes its parameters from the environment, but does not show how the first replicators assemble from polymers. A bottom-up replacement for that step is a candidate for stage 5.
- **Correlated environments.** Real settings correlate variables. Maps are drawn over the whole space first and then reported per class (4e), so that unlikely combinations do not mislead.

## 7. Decisions needed

1. **Depth of the model.** Hybrid (recommended): sub-models compute layer B from the environment, gates E1–E3, then the validated protocell simulation for L1–L4. The alternative is a full bottom-up artificial chemistry from monomers, which is closer to reality but much slower and may not converge in reasonable time.
2. **Extended ladder E1–E3.** Add it (recommended), or keep L1–L4 only.
3. **Order of sub-models.** As listed in 4b (recommended: clearest experiments first), or start with a particular scenario such as alkaline vents or hot springs.

## References

All checked on 25 October 2026 for the [condition review](condition_review.md); see the review for links and notes.

- Adamala, K. & Szostak, J. W. (2013). Nonenzymatic template-directed RNA synthesis inside model protocells. *Science* 342(6162), 1098.
- Baaske, P. et al. (2007). Extreme accumulation of nucleotides in simulated hydrothermal pore systems. *PNAS* 104(22), 9346–9351.
- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452.
- Hordijk, W., Kauffman, S. A. & Steel, M. (2011). Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. *International Journal of Molecular Sciences* 12(5), 3085–3101.
- Martin, W., Baross, J., Kelley, D. & Russell, M. J. (2008). Hydrothermal vents and the origin of life. *Nature Reviews Microbiology* 6, 805–814.
- Patel, B. et al. (2015). Common origins of RNA, protein and lipid precursors in a cyanosulfidic protometabolism. *Nature Chemistry*. Volume and pages not yet confirmed from the original.

Further sources for parameter values will be verified against the originals before they are used.
