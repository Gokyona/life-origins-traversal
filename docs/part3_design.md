# Part III (design draft): where, when and how did life begin on the early Earth?

**Status:** design, 30 October 2026, Toronto time; scope decisions in Section 9. Nothing here has been run yet.
**Builds on:** the completed Stage 4 environment search ([paper](../paper/COMPASS_stage4.md)) and its main limitation (Section 7 of that paper).

## 1. The question

Stage 4 asked *which environments can give rise to life*. Part III asks a historical question:

> After the Earth formed, **where**, **when** and **by which chemical pathway** was life most likely to begin?

A simulation cannot prove what happened. It can give a probability distribution over place, time and pathway, under stated assumptions, and test that distribution against the geological record.

## 2. Why Stage 4 cannot answer it yet

Stage 4 has two gaps that matter here.

1. **Copying speed.** The model's template-directed copying is about 10⁹–10¹¹ times faster than measured nonenzymatic copying. At a literature rate, heredity never starts from scratch. Any "when" computed from the Stage 4 chemistry would therefore be wrong by many orders of magnitude.
2. **No sequence-determined function.** Sequences are copied, but what a sequence *is* does not change how well it, or its neighbours, reproduce. Without this link there is no selection for function, and Level 4 (Darwinian evolution) stays at zero.

So Part III must first repair the chemistry, then add Earth's history around it.

## 3. Structure: three layers

```
Layer 1  Early-Earth history      site classes s, their number N_s(t) and environment x_s(t)
              │                     (the 13 Stage 4 layer A variables, now functions of time)
              ▼
Layer 2  Chemistry at measured    λ(x): rate at which one site with environment x
         rates                      first reaches heredity (L2), and by which pathway
              │
              ▼
Layer 3  Synthesis                distribution of the time, place and pathway of the first origin,
                                    compared with the geological record
```

### 3.1 The synthesis in one formula

Treat the origin as a rare event that can happen at any site. With N_s(t) sites of class s, each with environment x_s(t), the total hazard is

  Λ(t) = Σ_s N_s(t) · λ(x_s(t)).

The probability that life has begun somewhere by time t is

  P(T ≤ t) = 1 − exp(−∫ Λ(t′) dt′),

starting from the last time the surface was sterilised (for example by a large impact). The share of Λ carried by each site class gives **where**; the pathway labels of the successful runs give **how**. Uncertain inputs (layer C constants, Earth-history parameters) are dimensions that COMPASS maps, so the answer is a distribution with honest error bars, not a single story.

## 4. Layer 1: early-Earth history

A parametric model, not a climate simulation. From the Moon-forming impact (about 4.5 Ga) to about 3.5 Ga, in steps of about 10 Myr:

- **Global state:** presence and depth of oceans, surface temperature, atmospheric CO₂ and H₂, surface UV, emergent land area, impact flux and sterilising impacts, hydrothermal heat flow.
- **Site classes:** deep alkaline vents, high-temperature black smokers, volcanic hot-spring fields, closed-basin lakes, impact-crater lakes and their hydrothermal systems, coastal and tidal pools, icy and freeze–thaw settings, open ocean surface, mineral pores.
- **Per class and time:** the number of sites and the distribution of the 13 layer A variables (temperature, pressure, pH, redox, UV, ions, feedstock, dilution, mixing, concentration mechanisms, cycling, mineral surfaces, time).
- **Feedstock sources by epoch:** for example a reducing atmosphere after large impacts, cyanide and sulfide chemistry in shallow water, delivery by meteorites.

Each uncertain quantity (for example how much land existed at 4.2 Ga) becomes a parameter with a range from the literature, as layer C did in Stage 4. Every range must be sourced and checked against the original.

## 5. Layer 2: chemistry at measured rates

### 5.1 Repair copying first (step P3-1)

Test, one at a time and then together, mechanisms the literature proposes to close the gap. Each enters with a measured rate range:

1. Activated monomers and bridged intermediates that speed up primer extension.
2. Templated ligation of short oligomers supplied by the environment.
3. Concentration in the liquid veins of ice (eutectic freezing).
4. Polymerisation in the dry phase of wet–dry cycles.
5. Adsorption and catalysis on mineral surfaces (for example clays).
6. Longer-lived duplexes at low temperature.

The output is a list of which combinations, at measured rates, let heredity start from scratch, and under which environments. If none does, that is itself the main result, and Part III reports how large the remaining gap is.

### 5.2 Sequence-determined function (step P3-2)

Each sequence folds into a structure by a simple pairing model. Structure decides function, starting with one function that feeds back on reproduction: ligase-like activity that joins pieces on a template. This creates selection for function and a route to Level 4.

### 5.3 Rare events

At real rates, the origin may be an event of, say, once per 10¹⁰ site-years. Direct simulation cannot see that. λ(x) must be estimated with rare-event methods (forward flux sampling, weighted ensemble or splitting), which run many short trajectories and push the promising ones forward. This is new engine work and the main use of the B70: many independent trajectory batches in parallel.

## 6. Layer 3: tests against the record

The record is used as a **test**, never as an input, so the model can fail it. Candidate checks:

- Liquid water by about 4.4 Ga (zircons).
- Possible biogenic carbon at about 4.1 Ga.
- The estimated age of the last universal common ancestor (about 4.2 Ga in recent molecular-clock work).
- Accepted biosignatures by about 3.7–3.5 Ga.

If the model puts the first origin far later than the evidence, either the chemistry is still too slow or a mechanism is missing. That is a useful result too.

## 7. What Part III can and cannot claim

- **Can:** which site classes and epochs carry most of the probability, how robust that is to each uncertain input, and which chemical pathway the successful runs use.
- **Cannot:** that life did begin there, or what the first organism was. The answer is "most likely under these assumptions".
- As in every Stage 4 report, each result will state which conclusions follow from assumptions and which are derived.

## 8. Work plan

| Step | Content | Compute |
|---|---|---|
| P3-0 | Tables of measured rates (copying, ligation, hydrolysis, duplex melting) and Earth-history ranges, each source checked | none |
| P3-1 | Copying repair: which mechanisms start heredity at measured rates | B70 |
| P3-2 | Sequence → structure → function, ligase feedback, L4 | B70 |
| P3-3 | Rare-event estimator for λ(x), validated on cases with known answers | B70 |
| P3-4 | Early-Earth history model (layer 1) | CPU, light |
| P3-5 | Maps of λ over the site classes; synthesis; sensitivity | B70 |
| P3-6 | Tests against the record; paper | none |

Existing run rules apply: memory pre-flight at worst case before every long run, checkpoints and resume, one GPU job at a time, VRAM cap 28 GB.

## 9. Decisions (owner, 30 October 2026, 12:23 Toronto time)

1. **Alphabet: 2 letters (0/1)**, as in Stage 4. Measured rates exist only for 4-letter RNA, so they are borrowed for the 2-letter chemistry. This is an explicit assumption and every report states it.
2. **"How" goes to Level 4** (Darwinian evolution), through protocells (L3), within the same Earth history.
3. **Time window: 4.5–3.5 Ga.**

## 10. Candidate references

Checked sources and values are in [part3_p30_sources.md](part3_p30_sources.md). The list below is the starting list.

To be verified before use, as for all earlier stages:

- Wilde et al. (2001), zircon evidence for oceans at 4.4 Ga, *Nature*.
- Bell et al. (2015), possible biogenic carbon in a 4.1 Ga zircon, *PNAS*.
- Moody et al. (2024), age and nature of the last universal common ancestor, *Nature Ecology & Evolution*.
- Prywes et al. (2016), nonenzymatic copying rates (already cited in the Stage 4 paper).
- Walton & Szostak (2016), imidazolium-bridged intermediates in primer extension, *JACS*.
- Attwater et al. (2010), ice as a medium for RNA replication (already checked for Stage 4e).
- Patel et al. (2015), cyanosulfidic chemistry, *Nature Chemistry*.
- Benner et al. (2020), impact-generated reducing atmospheres and RNA precursors, *ChemSystemsChem*.
- Allen, Warren & ten Wolde (2005), forward flux sampling for rare events, *Physical Review Letters*.
