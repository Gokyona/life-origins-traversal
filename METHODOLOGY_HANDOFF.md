# Methodology handoff: Life Origins Traversal

**Purpose.** A self-contained summary of the methods, findings, limits, working rules and open questions of this project. It is written so that a researcher or AI model with no prior context can understand how the work was done and continue it. Where details matter, it points to the file in this repository that holds them.

**Author of the research:** Gokyona Wahiro (one person; pen name), independent researcher, Toronto. ORCID [0009-0008-1514-8856](https://orcid.org/0009-0008-1514-8856).
**Status at handoff:** 31 October 2026 (Toronto time). Parts I and II (Stages 1–4) are complete and published as two preprints. A Stage 4 supplementary study (robustness and measured-rate chemistry) was stopped early on the author's decision and is published. The research is paused.
**Preprints:** Paper 1 (COMPASS method) [doi:10.5281/zenodo.23825371](https://doi.org/10.5281/zenodo.23825371) · Paper 2 (Stage 4 environments) [doi:10.5281/zenodo.23825838](https://doi.org/10.5281/zenodo.23825838) · Homepage [gokyona.github.io](https://gokyona.github.io)

---

## Contents

0. How to use this document
1. Research question and framing
2. Timeline of the project
3. COMPASS, the search engine
4. Part I: validating the engine on test models
5. From conditions to environments: layers A, B and C
6. Part II (Stage 4): the bottom-up chemistry
7. The ladder of life and its detectors
8. Search protocol in Stage 4
9. Validation practice
10. Main findings
11. Limits and assumption audit
12. Stage 4 supplementary study (measured rates, robustness)
13. Side product: the Gokyona Environment platform
14. Hardware, operations and workflow rules
15. Repository and data layout
16. Open questions and suggested next steps
17. Author's standing decisions and preferences
18. Glossary
19. Key references

---

## 0. How to use this document

**Reading order for a newcomer.**
1. This document (overview, rules, open questions).
2. [paper/COMPASS_method.md](paper/COMPASS_method.md) (Paper 1) and [paper/COMPASS_stage4.md](paper/COMPASS_stage4.md) (Paper 2). The submission copies are in [paper/submission/](paper/submission/) and the PDFs in [paper/pdf/](paper/pdf/).
3. The specifications of each Stage 4 layer: [docs/stage4c1_spec.md](docs/stage4c1_spec.md) to [docs/stage4c5_spec.md](docs/stage4c5_spec.md). Each has a revisions section recording every change and when it was made.
4. The stage reports in [results/](results/), one per step, each with a "which conclusions follow from assumptions" section.
5. For the supplementary study: [docs/part3_alphabet_transfer.md](docs/part3_alphabet_transfer.md), [part3/p31/REPORT_P31.md](part3/p31/REPORT_P31.md), [stage4/supplement_robustness/SUMMARY.md](stage4/supplement_robustness/SUMMARY.md).

**If you continue the work**, follow the rules in Sections 9, 14 and 17. They are the author's standing requirements, not suggestions. In short: pre-register validation criteria, report failures as failures, test memory at worst case before every long run, make every run resumable and bit-reproducible, verify every citation against the original, state which conclusions are built in by assumptions, and keep all reports in English in this repository.

---

## 1. Research question and framing

**The question.** Most origin-of-life simulations choose one set of conditions, run a model for a long time and wait for life-like behaviour to appear. A negative run then says little. This project asks instead:

> **Under which conditions does life emerge with high probability, and under which does it never emerge?**

The goal is a *map* (a phase diagram) of condition space: the boundary between chemistry that stays chemistry and chemistry that becomes self-maintaining, hereditary, compartmentalised and evolving. Conditions are tested in an **ordered, never-repeating** sequence designed so that each test removes as much uncertainty as possible. The author asked for this ordered, non-repeating traversal from the start; it is the core idea of the project.

**Two parts.**
- **Part I (Stages 0–3): build and validate the engine.** The conditions in these models are *test conditions*. Several of them (division size, parasite advantage, number of gene types, copying error rate) are properties of the chemical system, not of the environment. Their maps validate the engine and the life detectors. They do not say which environments give rise to life. The author insisted on this distinction.
- **Part II (Stage 4): the environment search.** Inputs are environmental variables only. The chemistry is fully bottom-up: nothing living is put in. Every level of life must emerge from monomers, polymers, catalysts, templates and amphiphiles.

**The original 12 conditions (P1–P12)** from the author's project plan: P1 energy flux, P2 temperature, P3 food supply, P4 dilution, P5 number of monomer types, P6 catalytic probability, P7 polymerisation/hydrolysis ratio, P8 copying error rate, P9 amphiphile fraction, P10 membrane permeability, P11 environmental cycling (wet–dry, hot–cold), P12 mixing strength. A review before Stage 4 ([docs/condition_review.md](docs/condition_review.md)) found this list neither sufficient nor independent and replaced it by 13 primary environmental variables (Section 5).

**Scope of Stage 4.** Liquid water, carbon chemistry, the surface or sea floor of an Earth-like planet. Stage 4 asks "can life arise in this environment?", not "how long does it take?" and not "what happened on Earth?".

**What a result means.** Every map is a statement about the model's chemistry under stated assumptions, not about the history of the Earth. Each report separates conclusions that are *built in* by the model's rules from those that are *derived*.

---

## 2. Timeline of the project

All dates 2026, Toronto time.

| Date | Step |
|---|---|
| 25 Oct | Stage 1: first engine loses to random search (negative result). 1b–1d: revised engine beats random and Sobol sampling about 8× on a protocell chemostat; engine made 4× faster. 1e: autocatalytic sets (RAF), about 2× saving, threshold reproduced. |
| 25 Oct | Stage 2: stochastic corrector and coupled protocells (8× and 4× savings), multi-fidelity screening (partial negative), maps of L1–L3. |
| 25 Oct | Stage 3: knockouts, minimal sufficient conditions, Darwinian evolution (L4) with neutral-marker control, spatial cross-check, joint nested map, comparison with the straddle rule. Part I closed. |
| 25 Oct | Review of the condition list; 13 layer-A variables approved. Stage 4 plan; author chooses a fully bottom-up chemistry. |
| 26 Oct | 4a: engine at 20–30 dimensions. 4c-1 polymers (map). 4c-2 catalysis (full and focused maps; survived a real power cut). |
| 27–28 Oct | 4c-3 template copying and heredity (L2): 7 confirmed heredity environments. 4c-4 amphiphiles and vesicles: 4 search-found environments with heredity inside vesicles. |
| 29 Oct | 4c-5 coupling (L3) validation (6 of 8 tests passed) and a point study; no Darwinian evolution (L4 = 0). Author stops new experiments after step 5a. 4d-seq (sequence of levels) and 4e (real settings) done as analyses. Both papers drafted. Platform "Gokyona Environment" built and explored. |
| 30 Oct | Engine checks find that copying in the model is 10⁹–10¹¹× faster than measured chemistry. Platform v1.0 released. Part III designed (early Earth; measured rates). P3-0 literature tables; alphabet transfer analysis (no real base pair fits the heredity constants). P3-1 measured-rate copying model built and validated (9 of 9 tests). Both papers on Zenodo. |
| 31 Oct | Robustness study (P3-1b) of the 7 heredity points: heredity depends in part on the order in which the copy step fills scarce chain slots. Stopped early to save compute. Paper 2 revised. The work was reframed as a "Stage 4 supplementary study" and closed. |

---

## 3. COMPASS, the search engine

**Name.** COMPASS = Condition-space Ordered Mapping by Prioritised Active Sampling of Stochastic systems. Full description: Paper 1, Section 2.

**Problem.** A model has *d* conditions scaled to the unit cube. One simulation at condition *u* returns a binary outcome (level reached or not), so each condition has a success probability *p(u)*. Goal: classify conditions as life (*p* > 0.5) or not, with as few simulations as possible.

**Rules (final engine, default "engine4_soft --fast", plus the joint map from Stage 3e):**
- **R1, never repeat.** Every tested condition is quantised to a canonical key (4,096 levels per dimension) and archived. A key already in the archive is never simulated again.
- **R2, ordered coverage.** Candidates come from one scrambled Sobol sequence of 2¹⁹ points. The engine starts with the first 2,048; in every later round 30% of the batch continues along the sequence.
- **R3, one run per condition.** Repeated runs per condition were consistently wasteful; the map-maker averages neighbouring single runs better.
- **R4, shared map-maker.** A small neural network (two hidden layers of 256 units), binary cross-entropy, early stopping on a 10% held-out split.
- **R5, refinement by disagreement.** Five networks trained on bootstrap resamples score every untested candidate by the spread of their predictions. Spread is high where data would change the map and low where the outcome is just random.
- **R6, soft selection.** 70% of each batch is drawn with probability ∝ (spread)².
- **R7, growing batches.** Each round adds max(1,024, N/3) points, so rounds grow logarithmically with budget. The five networks train together as one batched network.
- **R8, joint map for nested levels.** When a higher level implies the lower ones, one network predicts all levels with a nesting constraint (zero contradictions, same accuracy, 3× faster training).
- **Masking (Stage 4).** Environments outside a focused region's constraints stay in the Sobol sequence but are never simulated.

```
pool <- scrambled Sobol sequence (2^19 points); archive <- {}
simulate the first 2048 pool points once each
while runs < budget:
    batch <- max(1024, runs / 3)
    30% of batch: next untested pool points in sequence order
    fit 5 bootstrap networks on all (u, outcome)
    score untested pool points by the spread of the 5 predictions
    70% of batch: draw with probability ~ score^2
    simulate the new points once each; archive their keys
final map <- early-stopped network on all data
```

**Benchmark protocol (essential for any claim of saving).**
- **Ground truth:** 4,000–8,000 held-out conditions, each simulated 64–128 times (Stage 4: 500–2,000 conditions × 16–32 seeds).
- **One map-maker for all methods.** Each method only produces an ordered list of simulated conditions. At fixed budgets the same map-maker is fitted to each list's prefix and scored on the held-out set. Methods differ only in *which* conditions they chose. (An early apparent win came from an underfitted map-maker for the baselines; see Section 4.)
- **Hard set:** conditions with 0.1 ≤ *p* ≤ 0.4 or 0.6 ≤ *p* ≤ 0.9 are reported separately.
- **Repetitions:** 3–5 per method for test models; mostly 1 in Stage 4.
- **Saving metric:** runs needed to match the baseline's accuracy; in Stage 4 also "runs to reach 95% of the achievable gain over a trivial (majority-class) classifier".
- **Equal-budget and distinct-hit checks** for rare levels (Paper 1, Section 6): truncate COMPASS to the ground truth's budget; compare distinct hits with the expected hits of a single-seed random sampler. The ground truth's 16 seeds per point make it look worse than it is.

**What was learned about the engine.**
- The saving is 1.6–8×, largest when the boundary is thin and the life region small, smallest in mostly negative spaces.
- It held at 20–30 dimensions (2.8–4.0×); engine overhead (~75 s per 64k-run map) does not grow with dimension.
- **The loop, not the scoring formula, gives the saving.** The straddle rule (Bryan et al. 2005) inside the same loop ties with COMPASS's disagreement² rule.
- COMPASS is a **map-maker, not a hit-finder.** It does not raise hits per simulation uniformly. Its unique gain is reaching levels too rare for random sampling to measure (below about 1 in 16,000 runs). For finding rare positives, a hit-oriented acquisition rule would likely be better.
- Negative results that shaped the engine: refining by *p*(1 − *p*) chases intrinsically random regions; repeated runs per condition waste budget; multi-fidelity screening helped only at small budgets (≤ 8k), because the cheap model was optimistic near the boundary.

Code: [stage1/](stage1/), [stage1b/](stage1b/) (`bench1b.py`, `surrogate.py`, imported by the Stage 4 map runners via `LOT_STAGE1B`), [stage2/](stage2/), [stage3/](stage3/), [stage4/s4a_scaling.py](stage4/s4a_scaling.py).

---

## 4. Part I: validating the engine on test models

| Stage | Model | Conditions | Outcome |
|---|---|---|---|
| 1 | protocell chemostat, Poisson tau-leaping, built-in error threshold | 6 | first engine lost to random (negative) |
| 1b–1d | same, harder; shared neural map-maker | 10 (2 distractors) | ~8× fewer simulations after fixing the map-maker fairness; engine 4× faster; soft selection default |
| 1e | Kauffman/Hordijk–Steel RAF sets, exact maxRAF on GPU | 7 | ~2×; RAF threshold reproduced within 2% (Hordijk, Kauffman & Steel 2011); map threshold within 7% |
| 2a | stochastic corrector (Szathmáry & Demeter 1987), 128 cells | 8 | ~8× |
| 2b | coupled protocell (self-made membrane triggers division; after Gánti) | 10 | ~4× |
| 2d | multi-fidelity (7.8× cheaper screen) | 10 | partial negative |
| 2e | maps of L1–L3 | 10 | COMPASS > Sobol at every level; levels nested |
| 3a | knockouts | – | compartments and selection between cells are essential |
| 3b | minimal sufficient boxes (≥ 95% life) | – | L1 17%, L2 6.9%, L3 1.9% of space |
| 3c | Darwinian evolution with 8 neutral markers as drift control | – | L4 in 11% of space, all inside L3 |
| 3d | independent spatial model (16×8 grid, leak) | 12 | conclusions hold; leak > 0.7–0.8%/step kills heredity |
| 3e | joint nested map; straddle comparison | – | joint map adopted; straddle ties |
| 4a | padded to 20/30 dimensions | 20, 30 | 3.6–4.0× (distractors), 2.8–2.9× (weak dims) |

**Methods used in Part I that carry forward:** distractor dimensions (to test that the engine ignores irrelevant inputs), knockouts (to find necessary mechanisms), minimal boxes confirmed by direct simulation, neutral markers as an internal drift control for evolution, and an independent second model as a cross-check. Reports: [results/stage1*](results/), [stage2*](results/), [stage3*](results/); papers [paper/COMPASS_stage1.md](paper/COMPASS_stage1.md) to [paper/COMPASS_stage3.md](paper/COMPASS_stage3.md).

---

## 5. From conditions to environments: layers A, B and C

The condition review ([docs/condition_review.md](docs/condition_review.md)) split the inputs into three layers. This structure is the main methodological change between Parts I and II.

- **Layer A: primary environmental variables, the only inputs.** 13 variables encoded as 19 dimensions: temperature, pressure, pH, redox (H₂ availability), UV flux, ion composition (salinity, Mg²⁺), feedstock composition and supply, dilution and flow (residence time), mixing and geometry (exchange between sites), concentration mechanisms (evaporation, thermal-gradient accumulation), wet–dry and hot–cold cycling (depth, period), mineral surfaces (amount, type), available time. Ranges and sources: [docs/stage4a_ranges.md](docs/stage4a_ranges.md); they span acid springs to alkaline vents and flowing channels to closed pools.
- **Layer B: derived chemical properties, never set directly.** Catalytic probability (P6), polymerisation/hydrolysis ratio (P7), copying error (P8), membrane permeability (P10), and rates in general are computed from layer A through known laws (Arrhenius, mass action, Li & Breaker 1999 hydrolysis for nucleotide-like chains, Radzicka & Wolfenden 1996 for peptide-like chains, critical aggregate concentration, thermophoretic concentration).
- **Layer C: unknown constants as explicit map dimensions.** Rate constants that no measurement fixes (condensation rate, catalytic probability per chain, template speed-up k_t, error rate ε, pairing enthalpy ΔH_bp, duplex melting temperature T_m8, membrane stability limits, and so on) are dimensions with literature-based ranges. Their number grows from 10 (4c-1, 29 dimensions) to 35 (4c-4, 54 dimensions); P3-1 has 57. **Rule:** a layer-C range is widened only with literature support, never to force a result.

Why: treating derived quantities as independent inputs would put impossible combinations on the map (hot but error-free) and spread the effect of one variable over several inputs, which corrupts importance analysis.

---

## 6. Part II (Stage 4): the bottom-up chemistry

Plan: [docs/stage4_plan.md](docs/stage4_plan.md); sketch [docs/stage4_sketch.md](docs/stage4_sketch.md). The simulator was built in five layers. Each layer was validated before the next was added, and **adding a layer never changes the layers below** (checked by bit-identical regression: with the new feature switched off, the new code reproduces the previous layer's stored results field by field).

**General structure.** Each run has 64 reaction sites (pores, droplets or small pools) with 512 chain slots each. All sites share one environment and differ only by chance; they exchange material according to the mixing variable. Time stepping is adaptive tau-leaping; fast cycles are cycle-averaged within a step. GPU-batched in PyTorch (Intel XPU backend), many runs in parallel. The environment is constant within a run except for square-wave wet–dry or hot–cold cycling. Molecules are represented by type and sequence, not by chemical structure: chemistry 0 is nucleotide-like, chemistry 1 amino-acid-like, plus amphiphiles.

| Layer | Adds | Key rules | Spec | Code |
|---|---|---|---|---|
| 4c-1 polymers | feedstock → monomers → activated monomers (by UV, H₂ or drying) → chains | chains start and grow, hydrolyse, wash out, adsorb on minerals, break under UV; untemplated growth capped at 4 letters per chain per step (G_CAP) | [spec](docs/stage4c1_spec.md) | `stage4/chem4c1.py` |
| 4c-2 catalysis | chains of length ≥ L_cat may catalyse joining of two other molecules by a fixed random (hash) rule | RAF detection with persistence; catalysis among "abundant species" of a site | [spec](docs/stage4c2_spec.md) | `chem4c2.py`, `cat4c2.py`, `raf4c2.py` |
| 4c-3 template copying | nucleotide-like chains of a **2-letter alphabet (0/1)** copied letter by letter into the complement with error ε | primed initiation; copy extension with processivity q = r/(r + k_off(l,T)); copies stay paired until the duplex melts (two-state ΔG = (l−1)(ΔH_bp − TΔS_bp), set by T_m8 and k₀ = 10⁶/h); one copy start per template per step; chains capped at 64 letters (int64 packing) | [spec](docs/stage4c3_spec.md) | `chem4c3.py` |
| 4c-4 amphiphiles, vesicles | amphiphiles above a critical concentration form vesicles inside a stability window | window: \|pH − pKa_app\| ≤ w_pH, Mg²⁺ and salinity below limits, hot phase below T_mem; random enclosure at formation and rehydration; osmotic growth factor β for filled vesicles (Chen, Roberts & Szostak 2004); division by shear or cycling with random sharing (Zhu & Szostak 2009); chains ≥ 4 letters cannot leave (Mansy & Szostak 2008) | [spec](docs/stage4c4_spec.md) | on the B360 (see Section 15) |
| 4c-5 coupling | no new chemistry | vesicle lineage record; L3 and L4 detection only | [spec](docs/stage4c5_spec.md) | on the B360 |

Further Stage 4 analyses: **4d-seq** pooled all map runs and read the environmental range at which each level is reached ([results/stage4dseq_sequence.md](results/stage4dseq_sequence.md)); **4e** compared the top environments with real geological settings ([docs/stage4e_environment_classes.md](docs/stage4e_environment_classes.md)). The planned **4d full joint map** and **5b focused L3 map** were not run (author's decision; listed as limits).

---

## 7. The ladder of life and its detectors

Seven levels, each scored automatically against paired controls. Exact definitions are in the layer specs; the summary below gives the final versions.

| Level | Meaning | Detector (final) | Paired control(s) |
|---|---|---|---|
| E1 | activated feedstock at working concentration | concentration threshold | – |
| E2 | persistent chains of ≥ 8 units | persistence over the window | – |
| L1 | self-sustaining catalytic network | persistent RAF set (Hordijk & Steel 2004) that raises its products above the control; *L1 strict* and *loose* recorded; L1 nested in E1 only (not in E2, because short-oligomer RAFs appear where E2 fails) | same run without catalysis |
| L2 | heredity | lineage record of abundant sequences (≥ 6 letters) and complements at 16 checkpoints; information share persists ≥ 5 molecule turnovers; *strict* = share minus the max of two random-copy controls ≥ 0.15; **L2_long** = strict and mean copy length ≥ L_H = 6 (introduced after all 77 full-map "L2" hits proved to be short-copy noise) | two runs where copying makes chains of the same length with random sequence |
| E3 | persistent vesicles that hold chains | longest run ending at the end of the simulation spanning ≥ 5 vesicle turnovers | – |
| L3 | coupled protocell | *loose*: Spearman(relative growth, osmotic load Π) ≥ 0.1, permutation p < 0.01; *strict*: content passed on by copying (ancestor 3 generations up, ≥ 90% identity match in chains made after its last checkpoint, ≥ 0.5× letters, at ≥ half of sites, and above the random-copy control by ≥ 0.25) | β = 0 (no osmotic coupling); random-sequence copying |
| L4 | Darwinian evolution | heritability of a vesicle trait plus ρ(fitness, trait) above the 99th percentile of permutation and neutral-marker nulls; a founder test (6b) checks selection directly | drift / neutral markers |

**Detector lessons.** Detector definitions were revised several times during validation, always before any map result and always recorded (for example: L2 strict margin, L2_long, the E3 window, the L3 strict criterion revised three times). Paired controls are essential: e.g. the random-copy control showed that "L3 loose" measures growth coupling by mass, while "L3 strict" measures heredity.

---

## 8. Search protocol in Stage 4

- **Full map per layer:** 1,000 random environments × 16 seeds as ground truth plus 32,000 COMPASS runs (4c-1: 2,000 × 32 ground truth and 64,000 runs each for COMPASS and a matched Sobol arm; this is the only Stage 4 layer with a matched baseline).
- **Focused maps:** regions chosen from earlier maps (for example highly catalytic short chains; heredity boxes; membrane-compatible boxes). Environments outside the region are masked. Typically 500 × 16 ground truth and 16,000 COMPASS runs. Box3 v2 drew membrane constants inside the stability window.
- **Targeted scans:** two-stage scans (`scan4c3.py`) to find a natural reference point for validation when the level was too rare for the full map.
- **Confirmation:** from 4c-3 on, every candidate for a rare level was re-run with 16 new seeds; a point counts as confirmed at ≥ 12–13 of 16. Marginal points (for example 5/16, 7/16) were reported but not counted.
- **Importance and partial dependence:** refitted maps report variable importance and 1-D/2-D partial dependence; stability across repetitions was checked (e.g. 4c-3 rep 0 vs rep 1 rank correlations).
- **Sensitivity checks on single laws:** e.g. hydrolysis floors at pH 5 and Mg²⁺ 5 mM (4c-1; the Li & Breaker law has no acid branch); an acid-artefact point (#813, 54 °C, pH 2.7) was dropped from 4c-3 for this reason.

---

## 9. Validation practice

These practices were applied throughout and are required for continuation.

1. **Pre-registration.** Every layer's validation tests and pass criteria are written into its specification *before* the tests run. Changes to criteria before any result are recorded in the spec's revisions section.
2. **Failures stay failures.** A failed test is reported as failed and its criterion is not changed after the fact. A redesigned test gets a new name (e.g. 4c-3 test 9 failed, test 9b with a criterion fixed beforehand passed; 4c-5 tests 3 and 6 failed, 6b passed).
3. **Stop rule.** After a failed pre-registered test, no maps run until the author decides. Anything run anyway is labelled a "recorded deviation" (the 4c-5 point study 5a).
4. **Bit-exact kill-and-resume.** Every program must survive a kill or power cut: atomic writes, JSONL results with fsync, checkpoints at most every 5 minutes, `--resume` with deterministic seeds, and a kill-and-resume test with bit-identical output (`stage4/kill_resume_test*.sh`, `compare_runs*.py`). 4c-2 survived a real power cut this way. GPU code must be deterministic (e.g. deterministic tie-breaking in `cat4c3.py`).
5. **Bit-identical regression across versions.** New features enter as switches that default off; with the switch off the new code must reproduce the old code's stored results in every field. This was done for every layer, the platform engine, P3-1 (85 fields vs chem4c3) and the supplementary study.
6. **Memory pre-flight before every long run.** Run the new or changed code at worst-case settings (the densest known point set, largest batch) and measure peak GPU memory (reserved and allocated), host memory and growth over at least 5 batches. Pass only if the peak is at most 75% of the caps and memory does not grow. Pick the batch size that passes (`stage4/preflight_mem.py`). Introduced after an out-of-memory crash that killed a run hours in.
7. **Controls and knockouts.** Paired controls for every detector (Section 7); null fixtures (random letters, ε = 0.5) to measure false positives; knockouts in Part I; distractor dimensions.
8. **Reproducing known theory.** RAF threshold (Hordijk, Kauffman & Steel 2011) in Stages 1e and 4c-2; Eigen's error threshold (loss near 0.5–1 errors per copy) in 4c-3 and P3-1; parasite takeover in 4c-3; exact steady-state chain-length distribution in 4c-1; nearest-neighbour duplex thermodynamics and Michaelis–Menten extension in P3-1.
9. **Ablations.** Switch one model rule at a time at fixed seeds and count which conclusions survive (supplementary study variants A0–A6, A2b–A2e; Section 12).
10. **Assumption audit in every report.** A table "which conclusions follow from assumptions" with verdicts *built in*, *partly built in*, *derived* (Paper 2, Section 6). A self-critique on 26 October found that of 20 early conclusions, 10 were built in by assumptions, 7 partly and 3 genuinely derived; since then every report carries this section.
11. **Citations verified against originals.** Each reference is checked against the publisher page, the authors' copy or a library record, with the date of the check. Unverified numbers are flagged as such. (doi.org and Crossref, and some PMC pages, were unreachable from the cloud environment; publisher or author pages were used instead.)
12. **Honest comparisons.** Fair map-maker for all methods; equal-budget and distinct-hit checks before any claim that "random sampling finds none"; numbers in papers audited against stage reports.
13. **md5 for every file.** Code and data are copied byte for byte from the run machine and checked by md5; MD5SUMS files sit next to the data.

---

## 10. Main findings

**Paper 1 (method).** COMPASS needs 1.6–8× fewer simulations than Sobol or random sampling for the same map accuracy on six models, holds at 20–30 dimensions, and its own cost does not grow with dimension. The saving comes from the loop, not the scoring rule. It found and confirmed environments for levels the ground truth never hit, but it is a map-maker, not a hit-finder.

**Paper 2 (Stage 4), within the model's assumptions:**
- The ladder narrows at every step: activated feedstock in ~45% of environments; persistent chains 8.5–16%; catalytic networks only with highly catalytic short chains (0.3% of a focused region); heredity (L2_long) in a narrow class (5 confirmed plus 2 from a repetition = 7); heredity inside vesicles at 4 search-found environments; coupled protocells at all 10 coexistence environments (298 of 435 coexistence runs); **no Darwinian evolution**.
- Catalytic networks do not need long chains; vesicles alone do not bring heredity.
- Heredity needs temperature cycling across the duplex melting point, molecules that outlive a cycle, and fewer than ~0.3 errors per copy.
- The environments at the top are cool to mild (0–25 °C), slowly flushed (residence times of years), enclosed water bodies with temperature cycles of days to weeks and near-neutral to slightly acidic water. Closest real setting: cool closed-basin lakes (Toner & Catling 2020). Hot-spring pools are too hot for copying in the model; alkaline vents do not match.
- No conflict between consecutive levels (114 comparisons), so one stable environment could carry the chemistry through all levels reached.
- The hardest step is having heredity and stable membranes in the same place (most heredity points are acidic, where no membrane forms).
- The L4 ceiling: copying families sweep once, then copying refills every vesicle to the same load, so no heritable variation remains. Prediction: sustained evolution needs a sequence-dependent trait that copying does not erase.
- Conflicts with experiment flagged: the top points are as salty as sea water (fatty-acid vesicles fail there), and three of four need membranes with pKa_app near 6 (measured mixtures form vesicles from pH ~6.5).

---

## 11. Limits and assumption audit

These limits are part of the published results. Any continuation should treat the first four as the most important.

1. **Copying is 10⁹–10¹¹× faster than measured chemistry.** The model's template-directed extension rate k_on·[M*] multiplies the uncatalysed condensation constant by the template enhancement k_t and by all layer-A concentration factors (evaporation × thermal gradient × dry phase). At coexistence point engine12538: ~1.7 × 10⁶ per molecule per hour, ~5 × 10¹⁰ letters per hour. Measured nonenzymatic copying: 30 minutes to about a day per letter (Prywes et al. 2016). With a literature-anchored extension rate, template copying stops entirely (short duplexes dissociate ~10⁶× faster than a letter is added). Templated ligation of chance oligomers and chance long primers also fail. A semi-analytic estimate says heredity from scratch would need 8-mer duplexes melting above ~50 °C, untemplated growth reduced 7–27×, and a supply of 6–8-letter oligomers at 1–10% of the monomer flux. **All heredity levels (L2, E3-joint, L3) are conditional on fast copying.**
2. **No real base pair fits the heredity constants.** The two letters were not a named base pair; pairing strength, accuracy and speed were independent layer-C constants. Translated with nearest-neighbour parameters: G·C 8-mers have T_m8 ≈ 78–92 °C, A·U ≈ 2–25 °C. The heredity points (T_m8 11.9–30.1 °C) pair like A·U but need G·C-like accuracy (ε ≤ 0.0013–0.0074) and speed. In ~30,700 runs with T_m8 70–90 °C, no L2_long run appeared (~9 expected), because hot phases hydrolyse chains. ΔH_bp of 53–56 kJ/mol (G·C) was never sampled.
3. **Template-ordering (slot-allocation) artefact.** At the heredity points, 88–99% of copy attempts find no free chain slot. The engine fills free slots in slot-index order, so the same long-lived templates are served every step (top 10% of lineages get 81–91% of copies). With random order, L2 is lost at 2 of 7 points and halves elsewhere; with fair multi-start allocation by rate, L2 survives at none (with the melt-wait floor) or at 1–3 points (without). Only rep1_1926 is robust. **Heredity at the confirmed points is largely a property of the numerical copy rule.** E3, joint and L3 were not re-tested.
4. **Biased GPU Poisson sampler.** `torch.poisson` on the Intel XPU is biased and overdispersed at small means (−1.1% at λ = 0.02, z −8.6; per-seed means ±5–8%; χ² p down to 1e-287); clean on CPU. All Stage 4 and P3-1 runs used it. An exact sampler (Cornish–Fisher start, exact inversion) leaves L2 unchanged at the 7 points; untested at coexistence points. Earlier, `torch.poisson` returned garbage above ~4 × 10⁹ on the XPU (fixed in 4c-1). **Do not trust `torch.poisson` on the XPU.**
5. **One copy start per template per step.** Starts happen with probability 1 − exp(−λ·dt); at the coexistence points 99.1–100% of template-steps have λ·dt > 3 (median ≥ 10⁶), so copy throughput is set by step size and template count, not by template properties. A copy-rate trait could not be selected; this is part of the L4 ceiling.
6. **Untemplated growth cap (G_CAP = 4 letters per chain per step)** binds in ~90% of growable chain-steps at two of three coexistence points, leaving more activated monomer for copying than mass action would. 61–82% of chains sit at the 64-letter cap, so detectors are probably only modestly affected.
7. **Cycle-averaged duplex dissociation** is an arithmetic mean over wet and dry (or hot and cold) phases, which favours melting. Phase-resolved melting kept L2 at the 3 long-period points tested.
8. **No phenotype.** Sequences affect nothing but length and mass, so L4 = 0 is built in by omission. The 4c-2 catalysis rule is also not heritable (the hash is over the exact full sequence and a copy is the complement) and is effectively inactive in the heredity regime.
9. **64-letter cap** sets inherited sequence length (~60 letters), not whether heredity happens (caps 48 and 32 keep L2).
10. **Coarse membrane laws** (stability window, exponential critical concentration); copying does not depend on Mg²⁺ (so the Mg²⁺ conflict between copying and membranes is represented only on the membrane side); vesicle counts include empty remnants (1–2%).
11. **Constant environments**, no time-varying environment, no carry-over between levels; "available time" enters only as a ratio.
12. **Sampling limits:** the layer-C ranges decide much; rare-level results rest on few points; no full joint map; the L3 region's extent is not mapped; Stage 4 maps mostly have one repetition.

Full Chinese assumption list for external review: `/mnt/project-files/part3/assumptions_review_zh.md` in the project's shared files (each assumption tagged measured / argued / layer C / simplification / built in / numerical; Section 10 lists questions for reviewers). It is not in the repository; an English version would be a useful addition.

---

## 12. Stage 4 supplementary study (measured rates, robustness)

Begun 30 October 2026 as "Part III" ([docs/part3_design.md](docs/part3_design.md)), whose larger aim was a historical question: where, when and how life began on the early Earth (4.5–3.5 Ga), via (1) an early-Earth history model of site classes and their layer-A variables over time, (2) chemistry at measured rates, and (3) a synthesis Λ(t) = Σ_s N_s(t)·λ(x_s(t)), P(T ≤ t) = 1 − exp(−∫Λ dt), tested against the geological record (liquid water by ~4.4 Ga, LUCA ~4.2 Ga, biosignatures by ~3.7–3.5 Ga) and needing a rare-event estimator for λ(x). Only the first steps were done; files keep their "part3" names.

**P3-0, measured-rate tables** ([docs/part3_p30_sources.md](docs/part3_p30_sources.md)): bridged-dinucleotide extension k_obs 0.35–84 h⁻¹ (Ding et al. 2022); activated intermediates hydrolyse with t½ ~2 h, so supply of activated blocks is a bottleneck; in-place activation works better at −3 °C (Aitken et al. 2023); ligation ~100× slower than monomer extension; duplex k_on 10⁶–10⁷ M⁻¹ s⁻¹; little land before 3.5 Ga.

**Alphabet transfer** ([docs/part3_alphabet_transfer.md](docs/part3_alphabet_transfer.md)): averaged 4-letter rates cannot be used for a 2-letter system. Each run must name the base pair its letters stand for and use that pair's measurements (nearest-neighbour thermodynamics transfer exactly as a subset; extension rates as a subset; fidelity does not transfer from 4-letter averages). Candidates: G/C (fast, accurate, needs hot separation), A/U (melts cool, slow and inaccurate), intermediate pairs A:s2U, I:s2C, D·U (Fang et al. 2025 measured A:s2U ~1.1 kcal/mol more stable than A:U per pair, s2U*s2U copied ~28× faster than U*U). Fast pairs need hot strand separation; Stage 4 did not sample **brief, rare hot pulses** in a cool setting, which could separate strands at little hydrolysis cost. The author chose **A/s2U** for P3-1, with G/C and A/U as controls.

**P3-1, measured-rate copying model** ([part3/p31/](part3/p31/), [REPORT_P31.md](part3/p31/REPORT_P31.md), spec [docs/part3_p31_spec.md](docs/part3_p31_spec.md)): nearest-neighbour duplex energies per named pair, Michaelis–Menten extension, hot pulses on top of the cycle, a reagent (activation) pool with side-product inhibition calibrated on Aitken et al. 2023 (cold-gap analysis [docs/part3_cold_gap.md](docs/part3_cold_gap.md)). All 9 tests passed (including bit-identical Stage 4 regression and kill-and-resume). No map was run. Science flags: with measured thermodynamics, 2–3-letter primers rarely extend (a G/C 4-mer's k_off ≈ 5 × 10³ h⁻¹ at 5 °C against r ≈ 30 h⁻¹), so copying needs very high primer concentrations or longer primers; in 640 random A/s2U points none made a copy ≥ 6 letters. One cold-benefit datum (24 h single addition) remains unexplained (residual ×2.2 per letter).

**P3-1b, robustness of the Stage 4 heredity points** ([stage4/supplement_robustness/](stage4/supplement_robustness/), [SUMMARY.md](stage4/supplement_robustness/SUMMARY.md), implementation notes and all review decisions in `docs/IMPLEMENTATION_P31B.md`). Seven L2 points × 16 seeds, both random-copy controls, common seed. Base reproduces the published rates (Fisher p ≥ 0.46). Variants: A0 exact Poisson sampler; A1 growth cap 16; A2 multi-start copying; A2a slot-rounding fix; A2b random slot-cut order; A2c/A2d/A2e decompositions; A4 slower hydrolysis of paired strands; A5 phase-resolved melting; A6 size-dependent washout. Robust: A0, A1, A4, A5 (3 points), A6. Not robust: slot order and multi-start (Section 11, item 3). Not run before the stop: salt correction, E3/joint/L3 groups at coexistence points, detector perturbations. Conclusion for Paper 2: heredity in Stage 4 is a property of the model's copy rule, not evidence that real chemistry can start heredity from scratch. Paper 2 and its PDF were revised on 31 October; the Zenodo record still holds the earlier version.

---

## 13. Side product: the Gokyona Environment platform

[platform/](platform/) (v1.0, release tag `platform-v1.0`, MIT). A local, installable simulation platform built on the Stage 4 engine (`gk_chem`, bit-identical to chem4c4 at one word per chain, multi-word sequences lift the 64-letter cap), with a browser UI on localhost, device auto-detection (XPU, CUDA, MPS, CPU), VRAM auto-sizing and an honesty banner. ABOUT and README are English only.

Exploratory findings from building it (29–30 October; data on the B360 in `gokyona-env/test_results_2026-10-29/`, not in the papers except where noted):
- With an imposed sequence → trait rule, copy-rate and division traits could not be selected (copy-start cap, saturated division and slots); a retention trait gave one dose-dependent shift but no sweeps.
- The copy step rationed monomers in **slot order**, biasing copying towards some templates; random-order rationing removed the bias and raised the trait plateau (the same artefact later confirmed in P3-1b).
- Slowly changing environments (100-year trend between two landscapes) gave tracking adaptation; seasonal change gave generalists; spatial gradients gave local adaptation (clines) but no regional takeovers.
- Coupling growth to the trait by "fitter chains grab monomers faster" produced positive selection but was **rejected by the author as too artificial**: coupling must come from chemistry or physics.
- A physics-only coupling exists (osmotic load correlates with division, mother–daughter correlation 0.44), but all sequences are chemically equal, so nothing sequence-dependent raises the load.
- These checks found the 10⁹–10¹¹× copying speed (Section 11, item 1).

---

## 14. Hardware, operations and workflow rules

**Hardware.** One desktop PC ("B360"): Intel Core i5-9400 (6 cores), 64 GB DDR4, 2 TB NVMe, **Intel Arc Pro B70 (32 GB)**, PyTorch on the XPU backend. Full load: GPU ~65–72 °C, 100–160 W. Device order in code: Intel XPU → CUDA → CPU.

**GPU and memory rules** (from an out-of-memory crash on 26 October: three GPU jobs at once plus XPU caching-allocator growth made the Intel driver spill ~52 GB of VRAM into system RAM as unmapped shared memory, invisible in RSS, and the PC crashed):
- **One GPU job at a time** (GPU lock in `stage4/thermal.py`; a second GPU job exits with code 3).
- VRAM cap via `LOT_XPU_FRACTION` = 0.875 (28 GB; earlier 0.5). Worst-case peak at most 75% of the cap; watch for spill (VRAM used per `xpu-smi`).
- `torch.xpu.empty_cache()` per batch (`LOT_EMPTY_CACHE_EVERY=1`), torch threads 4, each job in its own `systemd-run --user` scope with `MemoryMax=24G`.
- Thermal guard (read-only on hardware): pause at the critical temperature − 10 °C / 90 °C until below 75 °C; CPU pause 90/75 °C; log to `gpu_temp.log`. Alert at ≥ 85 °C.
- `watchdog.sh` restarts a dead or stalled queued job; `resume_all.sh` restarts the watchdog after a power cut (run by hand; no autostart without the author's approval).
- Batch size is whatever passes the pre-flight (e.g. 4c-4: 256; P3-1: 512). Do not chunk to cheat the pre-flight.

**Workflow between sessions.**
- A cloud AI session wrote specifications, reviewed results, committed to GitHub and wrote reports. A Remote Control session on the B360 implemented, ran and checked code locally and **never pushed**; it passed files back, which were committed after md5 checks.
- Running newly written code on the B360 for a long job needed the author's go-ahead; long runs started only after the pre-flight passed.
- While long jobs run, check the machine every hour and act immediately on problems (a crash once went unnoticed overnight).
- Update GitHub promptly after every step (validation results, spec changes, README, homepage), not only at the end of a stage. On the homepage, update both the research log and the project card's milestone list.
- Results must be written to permanent disk on the B360, never to a session's temporary scratch folder.
- Long or interrupted Remote Control sessions may not resume; start a fresh one with self-contained instructions.

---

## 15. Repository and data layout

```
README.md                 overview, roadmap, results table
METHODOLOGY_HANDOFF.md    this document
CITATION.cff              citation metadata (Zenodo)
paper/                    Paper 1 (COMPASS_method.md), Paper 2 (COMPASS_stage4.md), stage 1-3 reports,
  submission/, pdf/       outlines; submission copies and PDFs (built by tools/build_submission.py)
docs/                     plans and specifications: condition_review, stage2/3/4 plans, stage4a_ranges,
                          stage4c1..c5 specs, stage4e environment classes, part3_* (design, sources,
                          alphabet transfer, cold gap, P3-1 and P3-1b specs)
results/                  one report per step (stage1 ... stage4c5, 4d-seq), with JSON outputs
stage1/ stage1b/ stage1e_raf/ stage2/ stage3/   Part I code (stage1b holds bench1b.py, surrogate.py used later)
stage4/                   Stage 4 code 4a-4c3 with md5 table (README), map and validation JSONs,
                          run-safety tools (ckpt.py, thermal.py, watchdog.sh, preflight_mem.py)
stage4/supplement_robustness/   P3-1b code (chem_p31b, chem4c4b, chem4c5b), results, docs
part3/                    alphabet_tm.py; p31/ (P3-1 code, validation, logs, literature notes)
platform/                 Gokyona Environment v1.0 (package, docs, tests, screenshots)
raw/stage3/               raw results of stage 3
tools/                    paper build scripts
```

**Not yet in the repository (on the B360, `/run/media/bern/Data2T/myClaude/`):**
- Stage 1–2 raw results (`lot-run/`, `lot-raf/`, `lot-scm/`), Stage 3 run folder (`lot-s3/`).
- Most Stage 4 raw data and the code of 4c-4 and 4c-5 (`chem4c4`, `chem4c5`, run and confirm scripts) in `lot-s4/` (upload folders `upload_*` with MD5SUMS) and `lot-p3/`.
- Per-paper packages with code, data and run commands: `papers/paper1_COMPASS_method/` (256 MB) and `papers/paper2_stage4_environments/` (491 MB), each with README, MD5SUMS and `setup_workspace.sh`.
- Platform test data: `gokyona-env/test_results_2026-10-29/` (with overnight findings).
The papers promise these uploads with checksums. Uploading them is the first housekeeping task.

---

## 16. Open questions and suggested next steps

Ordered by how much they would change the conclusions.

**A. Repair the numerical copy step, then re-test what Stage 4 claims.**
1. Make copy allocation fair: several starts per template per step, shared by rate; random or rate-proportional slot and monomer rationing; uncapped untemplated growth; phase-resolved melting; an exact Poisson sampler on CPU or a validated GPU sampler. Check dt-convergence of copy length and counts.
2. Re-run the 7 heredity points and the 10 coexistence points (E3, joint, L3) under these rules (the unfinished part of P3-1b). This decides whether any Stage 4 heredity result survives as a property of the chemistry.

**B. Can heredity start at measured rates?** (The P3-1 model is built and validated but no map was run.)
3. Primer problem: test realistic ways to start copying with measured thermodynamics: longer chance oligomers, hairpin self-priming, templated ligation, oligomers made in the dry phase or on minerals, ice/eutectic concentration.
4. Environments not yet sampled: brief, rare hot pulses (convection in porous rock, geyser and hot-spring margins). Find measured pulse lengths and frequencies.
5. Run P3-1 maps for A/s2U with G/C and A/U controls; report the remaining gap in orders of magnitude if heredity never starts.
6. Missing data: fidelity of A/s2U copying, all-A/s2U stacking parameters, the Mg/Na salt correction (Owczarzy 2008 not read), a primary source for G/C-only copying fidelity.

**C. Sequence-dependent function and Darwinian evolution (L4).**
7. Add a sequence → structure → function rule that is chemically plausible for 2-letter sequences (e.g. ligase- or replicase-like motifs, amphiphile synthesis) and let coupling to reproduction *emerge* from chemistry, not from an imposed fitness rule.
8. Make catalysis heritable (motif-based, complement-symmetric) and check that it is active at the relevant copy numbers.

**D. Part III, the historical question** (design in [docs/part3_design.md](docs/part3_design.md)).
9. Rare-event estimator for λ(x) (forward flux sampling, weighted ensemble or splitting), validated on cases with known answers.
10. A parametric early-Earth model of site classes and their layer-A variables over time; then the synthesis and tests against the record.

**E. Engine (COMPASS).**
11. A hit-oriented acquisition rule for rare levels; comparison with Gaussian-process level-set methods at small budgets; matched Sobol or random arms and repetitions in more Stage 4-scale maps.

**F. Housekeeping.**
12. Upload the B360 code and raw data (Section 15) with MD5SUMS.
13. Upload the revised Paper 2 PDF as a new Zenodo version (the DOI is kept).
14. arXiv still needs an endorser (an endorsement request draft exists in the project's shared files).
15. Fix a known inconsistency: the 4c-3 and 4c-4 reports state the dimension count as "13 + 41", while the papers use 19 environmental dimensions + 35 constants.
16. Translate the assumption list for external review (Section 11) into English and add it to `docs/`.

**Questions for an external reviewer** (from the assumption list): Is A/s2U a sound 0/1 chemistry, or is there a measured 2-letter copying system? How did copying most likely start (long chance oligomers, hairpins, ligation, minerals, ice)? What explains the residual cold benefit in Aitken et al. 2023? Which early-Earth settings give brief hot pulses? What else could speed nonenzymatic copying by orders of magnitude? How to reconcile Mg²⁺ needs of copying with membranes beyond citrate? What is the simplest credible structure–function rule for 2-letter sequences? Do other numerical rules silently decide conclusions? Are the L2–L4 detectors biased? Which assumptions are missing?

---

## 17. Author's standing decisions and preferences

Respect these unless the author changes them.

- **Environmental framing.** "Conditions" means environmental conditions (originally P1–P12, now the 13 layer-A variables). Results over model-internal parameters must be labelled as such.
- **Fully bottom-up.** Nothing living is put in. Mechanisms must come from chemistry and physics; imposed rules (e.g. "fit chains grab monomers faster") are not acceptable. Prefer emergent mechanisms and scrutinise every artificial rule.
- **Two-letter alphabet (0/1).** The research question is whether the simplest alphabet can give rise to life ("2-letter life"). A 4-letter comparison is a check, never a replacement, and needs the author's approval before it is scheduled. Since the alphabet transfer analysis, the two letters must stand for a named base pair with that pair's measured properties.
- **Ladder to L4.** The aim includes Darwinian evolution, not only heredity.
- **Honesty.** Report negative results and failed tests; separate built-in from derived conclusions; do not over-claim (e.g. "random sampling finds none" needs an equal-budget note); state when a number is unverified.
- **Publication.** Publish only to this repository (GitHub `Gokyona/life-origins-traversal`), in English. The homepage is gokyona.github.io. All code, data and failed tests are public. AI use is disclosed in every paper.
- **Authorship.** Gokyona Wahiro is one person (never "they/their team"). Public texts present research interests and directions, not personal background.
- **Compute.** Use compute where it buys a result; the author stops runs early when the answer is clear enough. Times are given in Toronto time.

---

## 18. Glossary

| Term | Meaning |
|---|---|
| COMPASS | the ordered, non-repeating active-learning map engine (Section 3) |
| layer A / B / C | environmental inputs / derived properties / unknown constants as map dimensions (Section 5) |
| E1, E2, L1, L2, E3, L3, L4 | levels of the ladder (Section 7) |
| L2_long | strict heredity with mean copy length ≥ 6 letters |
| joint level | L2_long and E3 together (heredity inside vesicles) |
| coexistence points | the 10 environments where heredity and filled vesicles coexist |
| T_m8 | temperature at which an 8-pair duplex dissociates at rate k₀ = 10⁶/h; sets duplex stability with ΔH_bp |
| k_t | template speed-up relative to uncatalysed condensation (layer C) |
| ε | per-letter copying error |
| β | osmotic growth factor of filled vesicles |
| G_CAP | cap of 4 letters of untemplated growth per chain per step |
| ground truth (GT) | random environments each run with 16–32 seeds, held out for scoring |
| focused map / box | a COMPASS map restricted (masked) to a region chosen from earlier maps |
| B360, B70 | the run PC and its GPU |
| RC | Remote Control session running code on the B360 |
| P3-0, P3-1, P3-1b | the steps of the supplementary study (Section 12) |
| A/s2U | adenine / 2-thiouridine, the base pair chosen for measured-rate 2-letter copying |

---

## 19. Key references

All were checked against the original or a publisher, library or author record; check notes are in the stage reports and papers.

- Aitken, Wright, Radakovic & Szostak (2023). *JACS* 145(29), 16142–16149. In-place activation; cold benefit.
- Bryan, Nichol, Genovese, Schneider, Miller & Wasserman (2005). Active learning for identifying function threshold boundaries. *NeurIPS* 18, 163–170. (straddle rule)
- Chen, Roberts & Szostak (2004). *Science* 305, 1474–1476. Competition between model protocells.
- Ding, Zhou, Giurgiu & Szostak (2022). *Nucleic Acids Research* 50(1), 35–45. Bridged-dinucleotide extension rates.
- Eigen (1971). *Naturwissenschaften* 58, 465–522. Error threshold.
- Fang, Jia, Xing & Szostak (2025). *PNAS* 122(21), e2505720122. A/s2U, I/s2C alphabet.
- Hordijk & Steel (2004). *J. Theor. Biol.* 227, 451–461; Hordijk, Kauffman & Steel (2011). *Int. J. Mol. Sci.* 12, 3085–3101. RAF sets and thresholds.
- Li & Breaker (1999). *JACS* 121, 5364–5372. RNA hydrolysis law.
- Mansy & Szostak (2008). *PNAS* 105, 13351–13355; Zhu & Szostak (2009). *JACS* 131, 5705–5713. Protocell membranes.
- Prywes, Blain, Del Frate & Szostak (2016). *eLife* 5, e17756. Measured nonenzymatic copying rates.
- Szathmáry & Demeter (1987). *J. Theor. Biol.* 128, 463–486. Stochastic corrector.
- Toner & Catling (2020). *PNAS* 117(2), 883–888. Closed-basin carbonate lakes.
- Xia et al. (1998). *Biochemistry* 37, 14719–14735 (Turner 2004 nearest-neighbour set, NNDB).
- Full lists: Paper 1 and Paper 2 reference sections; [docs/condition_review.md](docs/condition_review.md); [docs/part3_p30_sources.md](docs/part3_p30_sources.md).

---

*Prepared on 31 October 2026 from the repository, both papers, the project's shared files and the project's working notes. The author used an AI assistant, Claude (Anthropic), to write this summary.*
