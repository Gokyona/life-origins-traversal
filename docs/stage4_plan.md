# Stage 4 plan: which environments give rise to life (bottom-up)

**Status:** revised 25 October 2026 after the owner chose a fully bottom-up model (decision at 23:54 Toronto time) · **Part II, first stage** · **Builds on:** [condition review](condition_review.md), [stage 4 sketch](stage4_sketch.md), [stage 3 paper](../paper/COMPASS_stage3.md)

## 1. Goal

Part I built and validated COMPASS on models whose conditions describe the chemical system (cell size, error rate, parasite advantage). Stage 4 asks the project's real question:

> **Given a local environment, how likely is it that chemistry there reaches each level of life?**

Stage 4 delivers two maps (owner's request, 26 October 2026):

1. **Environment map:** for every combination of environmental conditions, the probability of reaching each level of life.
2. **Composition map:** for the environments where life arises, what the first living systems are made of, and how much this varies between runs.

The only inputs are the 13 approved environmental variables (layer A). **Nothing living is put into the model.** Monomers, polymers, catalysts, templates, membranes and cells must all arise in the simulation from simple starting materials. The Part I protocell model is no longer the core of the simulation. It serves only as a reference: where the bottom-up chemistry happens to produce protocells, their behaviour should agree with Part I.

## 2. Architecture

```
layer A: 13 environmental variables (inputs, set independently)
   │   rate laws: known physical chemistry with measured parameters,
   │   each checked against experiments
   ▼
rate constants for every reaction class
   (activation, polymerisation, hydrolysis, catalysis, templated ligation,
    amphiphile assembly, permeation, degradation, transport)
   │
   ▼
bottom-up artificial chemistry on the GPU
   monomers → polymers → catalysts / templates → amphiphile aggregates and vesicles → protocells
   │
   ▼
automatic detector: E1–E3, L1–L4 measured on whatever emerges
   │
   ▼
COMPASS maps over layer A (+ layer C uncertain constants)
```

**Layer C: uncertain constants.** Where a rate law needs a constant with no measured value, the constant becomes an extra map dimension over a published or argued range, instead of being fixed by us. Importance analysis shows whether a conclusion depends on it.

### 2.1 The simulated chemistry

The simulation tracks molecules explicitly and stochastically (tau-leap, as in Part I), in many small reaction sites. A site is a pore, a droplet, or a patch of surface, coupled to its neighbours by mixing and flow.

| species class | what it is in the model | made by | lost by |
|---|---|---|---|
| feedstock | small precursors (abstracted by type, not by exact chemistry) | supply (feedstock, flow), UV- or redox-driven synthesis | flow, degradation |
| activated monomers | monomers of an alphabet of size A (A = 2–4), able to condense | activation (energy from UV, redox or dehydration) | deactivation, hydrolysis |
| polymers | sequences over the alphabet, explicit up to a maximum length, with length classes beyond | condensation; templated ligation | hydrolysis (rate depends on temperature, pH and water activity) |
| catalysts | polymers whose sequence catalyses a reaction, assigned at random with a probability that rises with length (RAF-style), and mineral surfaces | arise by chance among polymers | as polymers |
| templates | polymers that speed up the ligation of their complement, with a sequence-dependent error per position | arise by chance; copying | as polymers |
| amphiphiles | chains above a length threshold that self-assemble | synthesis from feedstock, or supply | degradation, flow |
| vesicles | aggregates of amphiphiles above a critical concentration; take up molecules at formation; grow by taking up amphiphiles; divide under shear or cycling | self-assembly | breakdown (ions, pH, temperature) |

All rates come from layer A through the rate laws. No rule says "replicate", "divide when doubled" or "be a cell". Copying is templated ligation, division is a physical effect of growth and shear, and compartments are vesicles that formed by chance around whatever was nearby.

### 2.2 Life ladder, measured on emergent structures

| level | meaning | detector |
|---|---|---|
| E1 building blocks | activated monomers above a working concentration | concentration time series |
| E2 polymers | chains of length ≥ L_min form faster than they hydrolyse | length distribution |
| E3 compartments | vesicles form and persist | vesicle census |
| L1 self-maintenance | a set of catalysts that collectively catalyse their own formation from feedstock persists (an RAF set in the reaction network actually realised) | RAF detection (Hordijk & Steel) on the active network, plus persistence over turnover times |
| L2 heredity | sequence information is copied: lineages of sequences persist longer than the molecules carrying them | sequence lineage tracking, against shuffled controls |
| L3 coupled protocell | vesicles whose contents help them grow and divide, and which pass on that content | correlation of content with vesicle growth and division, across generations of vesicles |
| L4 Darwinian evolution | a heritable property rises beyond drift | neutral-marker control as in stage 3c (neutral sequence positions as markers) |

The joint nested map of stage 3e handles all seven levels.

### 2.3 Composition record

Because every molecule is simulated, each run also reports what emerged, at almost no extra cost:

- **Monomers:** which alphabet dominates (nucleotide-like, amino-acid-like, mixed) and in what shares.
- **Polymers:** length distribution, sequence diversity, the share of chains that act as templates or catalysts.
- **Catalytic network:** whether an autocatalytic set exists, its size, which reactions it covers, and whether minerals take part.
- **Membranes:** amphiphile chain lengths, vesicle sizes, how much of the population is enclosed.
- **Order of emergence:** which came first in that run: polymers, catalytic network, templates or membranes ("metabolism first", "replicator first", "compartment first", or together).
- **Content of the first protocells (L3):** what a typical successful vesicle contains.

The composition map is built from these records. Multi-output maps are trained over layer A in the regions where life arises. The successful environments are then clustered into types (for example "hot, cycling pool type" or "alkaline vent type"), each with its typical composition. Outcomes are stochastic, so each point of the map gives a distribution of compositions, not a single answer. That distribution is itself a result: it separates compositions that are forced by the environment from those that are matters of chance.

## 3. Ranges of the layer A variables (provisional)

The ranges cover the stated scope: liquid water, carbon chemistry, the surface or sea floor of an Earth-like planet. Stage 5 will drop this scope and look at other solvents and environments, including chemistries not based on carbon. Step 4a fixes them from published sources, and every value is checked against the original source.

| variable | provisional range or encoding |
|---|---|
| temperature | 0–150 °C (above 100 °C only where pressure keeps water liquid) |
| pressure | 1–1,100 bar (surface to the deepest ocean) |
| pH | 2–12 |
| redox | reducing to oxidising, as H₂ availability (log scale) |
| UV | none to full unshielded surface flux (log scale) |
| ion composition | salinity 0 to above sea water; Mg²⁺ and Ca²⁺ concentrations (log scale) |
| feedstock composition and supply | supply rate (log scale) plus composition fractions: nucleotide-like precursors, amino-acid-like precursors, amphiphile precursors |
| dilution and flow | residence time (log scale) |
| mixing and geometry | from well-mixed to pore-confined (Péclet number) |
| concentration mechanisms | evaporation factor and thermal-gradient factor (log scale) |
| cycling | amplitude and period of wet–dry and hot–cold cycles |
| mineral surfaces | catalytic surface availability (log scale) and type (none, clay, metal sulphide) |
| available time | translated into simulated turnover times |

In practice that is about 18–20 dimensions, plus about 5–10 layer C constants. This is roughly twice the largest dimension COMPASS has handled so far (12).

## 4. Steps

The simulator is built in layers. Each layer is validated against known results before the next is added, and each layer gives its own partial map. A partial map is already a result: it shows where chemistry stalls.

**4a. Framework and engine scaling.**
- Fix the ranges and encodings of layer A from sources.
- Test COMPASS at 20–30 dimensions: the stage 2b model padded with distractors and with many weak real dimensions. This test is running.
- Plan the computational budget. Bottom-up runs will cost seconds to minutes each instead of about 1 ms. That makes COMPASS's saving matter much more, and gives multi-fidelity screening (stage 2d) a real chance: short or small runs as the cheap model.

**4b. Rate laws (environment → rate constants), each checked against experiments.**
- **Temperature:** Arrhenius dependence of all rates.
- **Hydrolysis:** dependence on pH, temperature and water activity.
- **Wet–dry cycling:** polymerisation and encapsulation under cycling (Damer & Deamer 2020).
- **Membranes:** vesicle stability in salt, Mg²⁺ and pH. High Mg²⁺ breaks fatty-acid vesicles, and citrate protects them while still allowing RNA copying (Adamala & Szostak 2013).
- **Concentration:** thermal gradients concentrate molecules (Baaske et al. 2007).
- **Feedstock:** UV-driven precursor synthesis against UV damage (Patel et al. 2015), and redox-driven supply (Martin et al. 2008).
- **Catalysis:** catalytic probability against length and minerals.

Each law gets a short report giving the law, the parameter sources, the experiment reproduced, and which constants are layer C.

**4c. Bottom-up simulator, layer by layer.**
1. *Polymer chemistry (E1, E2).* Feedstock, activation, condensation, hydrolysis. Validate: length distributions against the analytic equilibrium of reversible polymerisation; longer chains under wet–dry cycling.
2. *Catalysis (L1).* Random catalytic assignment by sequence. Validate: the emergence of autocatalytic sets near the threshold known from theory (Hordijk, Kauffman & Steel 2011), as in stage 1e, but now within the simulated chemistry.
3. *Templating (L2).* Templated ligation with errors. Validate: an error threshold for information maintenance (Eigen 1971), and the takeover by short, fast parasites known from experiments and theory.
4. *Amphiphiles and vesicles (E3).* Self-assembly above a critical concentration, uptake at formation, growth, division by shear or cycling, breakdown by ions. Validate against the membrane rate laws of 4b.
5. *Coupling (L3, L4).* Nothing new is added; L3 and L4 are only detected. Validate: where protocells emerge, compare their behaviour with Part I, including the compartment and selection effects of stage 3a.

At each layer: measure the cost per simulation, run COMPASS on that layer's levels, and publish a short report.

**4d. Full maps.** COMPASS maps for E1–E3 and L1–L4 over layer A + layer C with the joint nested map. Composition maps over the life regions, with clustering into environment types and the spread of compositions per point. Also: importance of each environmental variable; minimal sufficient environments; environmental knockouts (no cycling, no minerals, no concentration mechanism).

**4e. Real environment classes.** Report the maps for real settings, each defined as a region of layer A from the literature:
- deep-sea alkaline vents;
- land hot springs with wet–dry pools;
- shallow lagoons and tidal pools;
- cold settings (ice, eutectic brines);
- mineral pores.

For each, we ask which levels are reachable and which variable limits it.

**4f. Robustness and paper.**
- Vary the detector thresholds.
- Check how much the conclusions depend on layer C.
- State testable predictions.
- Write the stage 4 paper and update the homepage and log.

## 5. Milestones and acceptance

| step | output | accept when |
|---|---|---|
| 4a | ranges with sources; scaling test; compute plan | COMPASS still beats Sobol at 20+ dimensions, or a screening step restores it |
| 4b | rate laws with reports | each reproduces its experimental target within stated error |
| 4c-1 | polymer chemistry, E1–E2 map | analytic length distribution reproduced |
| 4c-2 | catalysis, L1 map | autocatalytic-set threshold reproduced inside the simulated chemistry |
| 4c-3 | templating, L2 map | error threshold and parasite takeover reproduced |
| 4c-4 | vesicles, E3 map | membrane stability window matches 4b |
| 4c-5 | coupling, L3–L4 detection | emergent protocells behave as in Part I where they appear |
| 4d | full maps | accuracy on held-out ground truth reported; importance stable across repetitions |
| 4e | per-class report | each class traced to literature ranges |
| 4f | paper | references verified; predictions stated |

## 6. Risks

- **L3 and L4 may never emerge** within feasible run times. The honest result is then a map of how far chemistry gets (E1–L2) under each environment, and of the bottleneck. This is still new information. Part I's results tell us what the missing step would have to achieve.
- **Cost.** A bottom-up run is orders of magnitude more expensive than a Part I run. Mitigations:
  - COMPASS's saving, which grows when the boundary is thin;
  - multi-fidelity screening;
  - layer-by-layer maps that rule out regions early. For example, where E2 is impossible, L1–L4 need not be simulated;
  - volunteer computing later.
- **Abstraction.** Molecules are represented by type and sequence, not by exact chemistry. The rate laws carry the environment's effect. Where a law is uncertain, it becomes layer C.
- **Too many unknown constants.** If layer C dominates importance, the result is "the answer depends on X, which has not been measured": a pointer for experiments.
- **Correlated environments.** Maps are drawn over the whole space first and then reported per real class (4e).

## 7. Decisions

1. **Depth of the model: fully bottom-up** (owner's choice, 25 October 2026, 23:54). No living structure is put in; everything must emerge.
2. **Extended ladder E1–E3:** included.
3. **Two deliverables:** environment map and composition map (owner's request, 26 October 2026, 00:30).
4. **Order:** layer by layer as in 4c, starting with polymer chemistry, unless the owner prefers to start from a particular scenario.

## References

All were checked on 25 October 2026. For Adamala & Szostak, Baaske et al., Damer & Deamer, Martin et al. and Patel et al., the [condition review](condition_review.md) gives links and notes. Eigen and Hordijk et al. were checked for the stage 1 report.

- Adamala, K. & Szostak, J. W. (2013). Nonenzymatic template-directed RNA synthesis inside model protocells. *Science* 342(6162), 1098.
- Baaske, P. et al. (2007). Extreme accumulation of nucleotides in simulated hydrothermal pore systems. *PNAS* 104(22), 9346–9351.
- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452.
- Eigen, M. (1971). Selforganization of matter and the evolution of biological macromolecules. *Die Naturwissenschaften* 58(10), 465–522.
- Hordijk, W., Kauffman, S. A. & Steel, M. (2011). Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. *International Journal of Molecular Sciences* 12(5), 3085–3101.
- Martin, W., Baross, J., Kelley, D. & Russell, M. J. (2008). Hydrothermal vents and the origin of life. *Nature Reviews Microbiology* 6, 805–814.
- Patel, B. et al. (2015). Common origins of RNA, protein and lipid precursors in a cyanosulfidic protometabolism. *Nature Chemistry*. Volume and pages not yet confirmed from the original.

Further sources for parameter values will be verified against the originals before they are used.
