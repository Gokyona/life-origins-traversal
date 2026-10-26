# Review of the condition list before stage 4

**Date:** 25 October 2026 · **Status:** reviewed; the owner approved the 13 layer-A variables on 25 October 2026 (23:19 Toronto time) · **Next:** [stage 4 sketch](stage4_sketch.md)

Before the environment search starts, we checked the original list of 12 conditions (P1–P12) against two questions. Is the list sufficient and independent? And how is the link from conditions to life established? The answer to both was "not yet". The list has been revised.

## 1. Sufficiency

The list is not sufficient. It leaves out several variables that major origin-of-life scenarios treat as important, or covers them only indirectly:

| missing condition | why it matters | source |
|---|---|---|
| pH and pH gradients | Decides whether fatty-acid vesicles form and whether RNA is stable. The proton gradient between alkaline vent fluid and sea water has been proposed as the first energy source. | Martin et al. 2008 |
| redox conditions (e.g. H₂ and CO₂) | The core energy source of the hydrothermal-vent hypothesis; "energy flux" alone does not say which kind of energy. | Martin et al. 2008 |
| ion composition: salinity, Mg²⁺ and other divalent ions | Non-enzymatic RNA copying needs Mg²⁺, but high Mg²⁺ breaks down fatty-acid membranes. Divalent ions and high salt in sea water hinder membrane formation. | Adamala & Szostak 2013; Deamer, Damer & Kompanichenko 2019 |
| concentration mechanisms (evaporation, thermal gradients, pores) | Dilute feedstock does not react. Thermal gradients in mm-sized pores concentrate single nucleotides more than 10⁸-fold. | Baaske et al. 2007 |
| mineral surfaces | Catalysis, adsorption and compartments; transition-metal sulphides are discussed as early catalysts. | Martin et al. 2008 |
| chemical composition of the feedstock | The supply rate says how much, not what. Starting materials such as hydrogen cyanide and hydrogen sulphide decide which precursors can form. | Patel et al. 2015 |
| ultraviolet radiation | Both an energy source (it drives precursor synthesis) and a cause of damage. | Patel et al. 2015 |
| available time and system size | Whether a rare event happens depends on the number of "trials". | our reasoning |
| pressure | Changes the boiling point of water and reaction equilibria; deep sea and land surface differ strongly. | our reasoning |

**"Covering every known condition in the universe" cannot be achieved, and cannot be proven.** Unknown variables cannot be listed. We have one example of life. And the list already assumes water as solvent and carbon chemistry. We therefore state a scope instead of claiming completeness:

- **Scope:** water as solvent, carbon chemistry, the surface or sea floor of an Earth-like planet. Planet- and universe-scale conditions (stellar radiation, atmosphere, planet size) act through the local conditions at the reaction site. The list therefore holds local conditions only.
- **Finding gaps:** list the variables that each major scenario emphasises (deep-sea vents, land hot springs and pools, ice, mineral pores) and check each against the list. A new variable can be added later, and COMPASS's importance analysis tests whether it matters.

## 2. Independence

The list is not independent. The conditions also sit at different levels: about half are external environment, the rest are results of environment and chemistry together.

| condition | level | depends on |
|---|---|---|
| P1 energy flux | environment (split by kind: UV, redox, heat) | — |
| P2 temperature | environment | linked to heat sources |
| P3 food supply | environment | partly on energy (activating monomers) |
| P4 dilution | environment | linked to flow and mixing |
| P5 number of monomer types | environment (part of feedstock composition) | feedstock source |
| P6 catalytic probability | **derived** | monomer types, chain length, mineral surfaces, temperature |
| P7 polymerisation / hydrolysis ratio | **derived** | temperature, water activity (wet–dry cycles), pH |
| P8 copying error rate | **derived** | temperature, Mg²⁺ and other ions, monomer types, pH |
| P9 amphiphile fraction | partly environment, partly derived | feedstock, energy |
| P10 membrane permeability | **derived** | membrane composition, temperature, ions |
| P11 environmental cycling | environment | geological setting |
| P12 mixing | environment | flow, pore geometry |

Treating temperature and copying error rate as independent inputs would place physically impossible combinations on the map (hot but error-free). It would also spread the effect of temperature over several inputs, which confuses importance analysis.

**Fix: two layers.** Layer A holds primary environmental variables, the only inputs; each can be measured on site and set independently. Layer B holds derived chemical properties, computed from layer A and never set directly. Real environments correlate layer-A variables (hot springs are often hot, reducing and metal-rich). Maps are therefore drawn over the whole layer-A space first and then reported per real environment class.

**Approved layer A (13 variables):** temperature; pressure; pH; redox conditions; ultraviolet radiation; ion composition (salinity, Mg²⁺); feedstock composition and supply; dilution and flow; mixing and geometry; concentration mechanisms; wet–dry and hot–cold cycling; mineral surfaces; available time. **Layer B:** catalytic probability, polymerisation/hydrolysis ratio, copying error rate, membrane permeability, and the stage 2–3 model quantities such as membrane synthesis rate and vesicle size.

## 3. From conditions to life: a chain of four links

| link | content | basis | status in stages 1–3 |
|---|---|---|---|
| 1. definition | what counts as life: NASA's working definition, "a self-sustaining chemical system capable of Darwinian evolution"; our levels L1–L4 break it into automatically scored steps | convention; thresholds need a sensitivity analysis | defined |
| 2. environment → rates and concentrations | Arrhenius relations, mass action, polymerisation equilibria, critical vesicle concentration, permeability, thermophoretic concentration | known laws with measured parameters (where available) | **skipped: parameters were set directly** |
| 3. rates → life functions | error threshold (Eigen 1971); autocatalytic sets at about 1–2 reactions catalysed per molecule (Hordijk, Kauffman & Steel 2011); stochastic corrector (Szathmáry & Demeter 1987); chemoton coupling (Gánti 1975, 2003) | theoretical thresholds | reproduced in stages 1–3 |
| 4. functions → L1–L4 | emergence in simulation | simulation only; risk of circular design | guarded by distractors, knockouts and an independent model |

**Weak points so far:** the mechanisms were chosen from the literature, not derived; most parameters have no measured source; the models start from replicators and vesicles; nothing has been compared with experiments yet.

**Standards for stage 4:**
- Every arrow is labelled as one of: known law with measured parameters; theoretical threshold; model assumption with an uncertainty range. Unknown constants become extra conditions on the map, so their influence is visible.
- Sub-models must first reproduce known experiments before any map is drawn. Examples: high Mg²⁺ breaks fatty-acid membranes, while citrate protects them without stopping RNA copying (Adamala & Szostak 2013). Thermal gradients in pores concentrate nucleotides by more than 10⁸ (Baaske et al. 2007). Wet–dry cycles drive polymerisation and encapsulate polymers in membranes (Damer & Deamer 2020).
- The model makes testable predictions, not only maps.
- Detector thresholds are varied (three settings per level) to show how stable the maps are.

## References

Each entry was checked on 25 October 2026 against the page linked.

- Adamala, K. & Szostak, J. W. (2013). Nonenzymatic template-directed RNA synthesis inside model protocells. *Science* 342(6162), 1098. Citation checked in the [ScienceDaily reference](https://sciencedaily.com/releases/2013/11/131128141401.htm); the Mg²⁺/citrate finding in [Harvard Medical School news](https://hms.harvard.edu/news/protocell-progress).
- Baaske, P., Weinert, F. M., Duhr, S., Lemke, K. H., Russell, M. J. & Braun, D. (2007). Extreme accumulation of nucleotides in simulated hydrothermal pore systems. *PNAS* 104(22), 9346–9351. [Caltech record](https://authors.library.caltech.edu/records/095h5-wzq59).
- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7133448).
- Deamer, D., Damer, B. & Kompanichenko, V. (2019). Hydrothermal chemistry and the origin of cellular life. *Astrobiology* 19(12). [eScholarship](https://escholarship.org/content/qt1p93j20t/qt1p93j20t.pdf) (section 3.1).
- Eigen, M. (1971). Selforganization of matter and the evolution of biological macromolecules. *Die Naturwissenschaften* 58(10), 465–522. (Checked for the [stage 1 paper](../paper/COMPASS_stage1.md).)
- Gánti, T. (1975). Organization of chemical reactions into dividing and metabolizing units: the chemotons. *BioSystems* 7(1), 15–21; Gánti, T. (2003). *The Principles of Life*. Oxford University Press. (Checked for the [stage 2 paper](../paper/COMPASS_stage2.md).)
- Hordijk, W., Kauffman, S. A. & Steel, M. (2011). Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. *International Journal of Molecular Sciences* 12(5), 3085–3101. (Checked for the stage 1 paper.)
- Martin, W., Baross, J., Kelley, D. & Russell, M. J. (2008). Hydrothermal vents and the origin of life. *Nature Reviews Microbiology* 6, 805–814. [Full text](https://www.esalq.usp.br/lepse/imgs/conteudo_thumb/Hydrothermal-vents-and-the-origin-of-life.pdf).
- NASA Astrobiology, ["Forming a definition for life"](https://astrobiology.nasa.gov/news/forming-a-definition-for-life/): origin of the working definition (NASA Exobiology Discipline Working Group, about 1992; often credited to G. Joyce).
- Patel, B., Percivalle, C. et al. (J. Sutherland's group) (2015). Common origins of RNA, protein and lipid precursors in a cyanosulfidic protometabolism. *Nature Chemistry*, doi:10.1038/nchem.2202. Checked in [phys.org](https://m.phys.org/news/2015-03-chemists-riddle-life-began-earth.html) and [MRC LMB news](https://www2.mrc-lmb.cam.ac.uk/?p=7247); volume and pages not yet confirmed from the original.
- Szathmáry, E. & Demeter, L. (1987). Group selection of early replicators and the origin of life. *Journal of Theoretical Biology* 128, 463–486. (Checked for the stage 2 paper.)
