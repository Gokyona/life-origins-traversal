# Stage 4c-2: where do self-sustaining catalytic networks arise? (Intel Arc Pro B70)

**Date:** 2026-10-26 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Runs:** full map 07:48–09:17, detector check 09:17–09:24, focused map 09:24–13:49 Toronto time (power cut 10:17–12:43, resumed from checkpoint) · **Specification:** [stage4c2_spec.md](../docs/stage4c2_spec.md) · **Code:** [stage4/chem4c2.py](../stage4/chem4c2.py), [cat4c2.py](../stage4/cat4c2.py), [raf4c2.py](../stage4/raf4c2.py), [run4c2_map.py](../stage4/run4c2_map.py) · **Data:** [map4c2b_result.json](../stage4/map4c2b_result.json), [map4c2_focus_result.json](../stage4/map4c2_focus_result.json), [sens4c2_result.json](../stage4/sens4c2_result.json), [validate4c2.json](../stage4/validate4c2.json), [validate4c2_t1.json](../stage4/validate4c2_t1.json)

## Bottom line

We added catalysis to the bottom-up polymer chemistry of 4c-1: any chain of at least a few units may, by a fixed random rule, catalyse the joining of two other molecules. Then we looked for **self-sustaining catalytic networks (level L1)**: sets of molecules that together catalyse their own formation from the food, that persist, and that make more of their own products than the same run without catalysis.

- **Across the whole environment space, L1 is very rare.** It never appeared in the 16,000 ground-truth runs. COMPASS, which steers towards the boundary, found 39 runs in 32,000.
- **L1 needs short chains to be good catalysts.** In a focused map where the catalytic probability p₀ is 0.01–1 and catalysts are 3–5 units long, L1 occurs in 0.33% of runs (0.3% of environments).
- **Where it occurs:** high p₀, long residence time, short catalysts, slow decay of catalytic ability with length, much feedstock, strong UV, and cold water.
- **What it looks like:** a median of 27 reactions among short, almost purely nucleotide-like oligomers, with catalysts 3–6 units long.
- **How strong it is:** catalysed joins are only about 1.5% of all joins in the network's site, yet the network's products are about 1.5 times as abundant as without catalysis (median).
- **L1 does not need long chains.** 69% of L1 runs have no persistent chains of 8 or more (E2). The assumed ladder "long polymers first, then catalytic networks" does not hold in this model.
- **The rest of the map is unchanged from 4c-1** in its drivers: residence time, temperature, feedstock and energy decide E1 and E2. Catalysis by polymers or minerals has no visible effect on E2.
- **The L1 map itself did not learn a boundary.** Its accuracy (0.997) equals that of calling every point negative. The statements about where L1 occurs come from the ground-truth counts and partial dependence, not from the classifier.

## What was simulated

The 4c-1 chemistry (feedstock → monomers → activated monomers → chains; hydrolysis, washout, UV damage, cycling) with these additions (37 dimensions in all):

- **Catalysis rule.** A chain C of length ≥ L_cat catalyses the ligation A + B → AB if a fixed hash of (C, A, B) is below p_cat(L_C) = p₀ (1 − e^(−(L_C − L_cat + 1)/λ)). The rule is random but reproducible: the same molecules have the same chemistry everywhere.
- **Rates.** A catalysed ligation runs κ times faster than uncatalysed condensation per unit catalyst concentration.
- **Temperature.** Uncatalysed condensation gets an activation energy E_a (40–120 kJ/mol); catalysed ligation uses r·E_a with r = 0.3–1.
- **Minerals.** Clay and sulphide surfaces speed up every condensation by a factor up to 1,000.
- **Where catalysis acts.** Among the abundant sequences of each site (at least 2 copies, top 16).
- **Length cap** raised from 32 to 64.

New layer C dimensions: L_cat (3–8), p₀ (10⁻⁶–1), λ (0.3–30), κ (10–10⁶), E_a (40–120 kJ/mol), r (0.3–1), mineral factors for clay and sulphide (0.1–1,000).

**L1 detector.**
1. At 16 checkpoints, the abundant species of 8 sampled sites are stored.
2. A maximal RAF set (Hordijk & Steel 2004) is computed at each checkpoint.
3. A network is *persistent* if the reactions that belong to the RAF at every checkpoint of a window still form a RAF of at least 3 reactions, and the window spans at least 5 turnover times of the members.
4. *L1 loose* = persistent RAF with E1 and liquid water. *L1 strict* additionally requires the RAF products to be more abundant at the end than in a paired run with catalysis switched off (same condition and seed).

## Validation (all passed before the maps)

| test | criterion | result |
|---|---|---|
| 1. RAF threshold | P(RAF) rises steeply at 1–2 catalysed reactions per molecule (Hordijk, Kauffman & Steel 2011; stage 1e) | crossing of P = 0.5 at f = 1.047 (n = 6) and 1.271 (n = 8); stage 1e: 1.06 and 1.25 |
| 2. catalysis off | p₀ = 0: no RAF, even with strong mineral catalysis | 0 of 8 at the reference, 0 of 64 random conditions; mineral catalysis lengthens chains (mean 29.4 against 3.4 units) |
| 3. shuffled catalysis | re-hashing catalysis every step removes persistent networks | L1 1 of 8 with fixed catalysis, 0 of 8 shuffled; mean persistence span 4.8 against 2.6 turnover times where a RAF exists |
| 4. hydrolysis link | hotter or more alkaline water shortens member lifetimes and lowers L1 | tetramer lifetime 0.33 → 0.067 → 0.013 h (+10, +20 °C) and 0.33 → 0.035 → 0.0036 h (+1, +2 pH); L1 1/8 → 0 → 0 in both |

The first validation run led to four changes to the specification before any map: persistence judged per reaction instead of per species, "self-sustaining" judged on the RAF products, L1 not nested in E2, and activation energies (see [specification, Section 2.1](../docs/stage4c2_spec.md)).

**Detector sensitivity** (sens4c2, 16 repeats at the reference): with 4 sampled sites refreshed every 8 steps, L1 strict was found in 1 of 16 runs; with 8 sites, in 6 of 16, at 112 against 114–117 ms per run. The maps use 8 sites refreshed every step. Starting catalysis at t = 0 instead of after E1 changed nothing at the reference.

## Results

**Full map** (37 dimensions; ground truth 1,000 environments × 16 runs; COMPASS 32,000 runs):

| level | share of environments | run rate |
|---|---|---|
| E1 (activated monomers at working concentration) | 0.443 | |
| E2 (chains ≥ 8 persist; E1 and E2) | 0.095 | |
| L1 strict / loose | 0 / 0 | 0 / 0 |

E2 is lower than in 4c-1 (0.158). The causes are not yet separated: the new activation energy of condensation, a different set of test environments, and longer runs when catalysis is on. The partial dependence of E2 on the activation energy is flat, so it is not the obvious cause.

Quality flags: liquid water in 95.5% of runs, 92.7% reached the planned end, 87.4% stationary, 4.7% with slot overflow.

**Focused map** (p₀ 0.01–1, κ 10⁴–10⁶, L_cat 3–5; other dimensions full range):

| level | share of environments | run rate |
|---|---|---|
| E1 | 0.445 | |
| E2 | 0.096 | |
| L1 strict | 0.003 | 0.325% (52 of 16,000) |
| L1 loose | 0.005 | 0.53% (85 of 16,000) |
| L1 strict with E2 / without E2 | 0.001 / 0.002 | |

**Where L1 occurs** (focused map; partial dependence of L1 strict from the lowest to the highest end of each range):

| input | low end | high end |
|---|---|---|
| p₀ (catalytic probability, 0.01–1) | 0.0002 | 0.0146 |
| residence time (1 h – 10⁵ h) | 0.0000 | 0.0203 |
| L_cat (3–5) | 0.0093 | 0.0006 |
| λ (0.3–30) | 0.0086 | 0.0011 |
| feedstock supply | 0.0002 | 0.0076 |
| UV flux | 0.0008 | 0.0092 |
| H₂ | 0.0019 | 0.0065 |
| temperature (0–150 °C) | 0.0058 | 0.0011 |
| κ (10⁴–10⁶) | 0.0048 | 0.0029 |

Counts of the 52 strict-positive runs by thirds of the range:

| input | low third | middle third | high third |
|---|---|---|---|
| p₀ | 1 | 18 | 33 |
| nucleotide share of feedstock | 11 | 2 | 39 |
| mineral type (none / clay / sulphide) | 17 | 35 | 0 |
| pH | 26 | 11 | 15 |
| temperature | 16 | 15 | 21 |

The temperature counts are nearly even, while the partial dependence favours cold water; the counts are few, so the trend is weak. No L1 run had sulphide minerals.

**What the networks are made of** (52 strict-positive runs; 10th / 25th / 50th / 75th / 90th percentiles):

| quantity | values |
|---|---|
| RAF size (reactions) | 5 / 12 / 27 / 35 / 56 |
| persistence span (turnover times) | 6.8 / 7.9 / 8.3 / 18.7 / 19.6 |
| catalyst length | 3 / 4 / 5 / 6 / 6 |
| nucleotide share of RAF products | 0.52 / 0.98 / 1.0 / 1.0 / 1.0 |
| mean chain length in the run | 3.4 / 4.9 / 6.1 / 7.7 / 34.5 |
| catalysed share of joins in the RAF site | 0.5% / 0.6% / 1.4% / 2.9% / 5.8% |
| product abundance, catalysis / control | 1.06 / 1.29 / 1.54 / 1.99 / 2.62 |

Across all runs of the focused map, the median catalysed share of joins is 0 (90th percentile 0.09%); in E2-positive runs, 0.09%.

**L1 and E2.** 31% of L1 strict runs also have E2. Where both occur, the delay of L1 after E2 is bimodal: nearly simultaneous in some runs, thousands of hours later in others.

**Accuracy along the budget** (share of test environments classified correctly; overall / hard set):

| level | 2k | 4k | 8k | 16k | 32k |
|---|---|---|---|---|---|
| E1, full map | .885 / .60 | .918 / .64 | .931 / .64 | .953 / .72 | .957 / .78 |
| E2, full map | .927 / .73 | .928 / .73 | .952 / .73 | .950 / .67 | .955 / .60 |
| L1, focused map | .994 / .33 | .997 / .33 | .997 / .33 | .997 / .33 | .997 / .33 |

The E2 hard set has only 15–16 environments and the L1 hard set 3, too few to use. There was no Sobol arm in 4c-2, so no saving is claimed here.

**What matters for E1 and E2** (full map, drop in map accuracy when one input is scrambled): for E2, temperature 0.076, residence time 0.067, feedstock 0.055, UV 0.033, UV yield 0.029, H₂ 0.022, pH 0.016; every catalysis and mineral dimension is at or below 0.004. Partial dependence of E2 on p₀, κ, mineral availability and the mineral factors varies by at most 0.01. Polymer catalysis therefore does not help long chains form in this model, and mineral catalysis hardly does.

**Composition of E2 chains** (1,529 positive runs, full map) is as in 4c-1: the nucleotide share follows the feedstock (medians 0.17, 0.51 and 0.90 for the three thirds of feedstock nucleotide share) and falls in hot water (0.49, 0.51 and 0.006 for 0–50, 50–100 and 100–150 °C). 14.7% of chains are longer than 32 units, so the old cap of 32 would have bound.

**Power cut, survived.** The power failed at 10:17 during the focused map. Every result file ended with a complete line and the checkpoints loaded. Only the batch running at the moment was lost (≤ 512 runs, about a minute of compute). The run resumed at 12:43. The 17 finished batches were byte-for-byte unchanged and the selections of the first 6 COMPASS rounds identical. GPU temperature stayed at 64–74 °C with no thermal pause.

## Which conclusions follow directly from the assumptions

| conclusion | how it arises | verdict |
|---|---|---|
| L1 needs high p₀ and short catalysts | p₀ and L_cat are inputs that set how much catalysis exists; with little catalysis no RAF can form (test 1) | **built in**: the result is the threshold, not the direction |
| L1 networks are almost purely nucleotide-like | catalysis acts only among sequences with repeated copies; with alphabet 2 that means short, abundant oligomers; nucleotide-like chains are short-lived in the conditions where L1 forms, so they turn over and repeat | **partly built in**: may reflect the abundance filter and persistence measured in turnover units, not chemistry |
| L1 does not need E2 | short oligomers suffice for a RAF in this model; long chains reduce repeated copies | **partly built in**: follows from the abundance filter |
| residence time, feedstock and energy favour L1 | they raise the concentration of short oligomers, which the detector needs | **partly built in** |
| cold water favours L1 | partial dependence says so, but the counts by temperature third (16 / 15 / 21) do not | **not established** |
| catalysis is weak (1.5% of joins) yet raises products 1.5× | measured on the emergent network against a paired control | **derived** |
| polymer and mineral catalysis do not help E2 | measured over the full space | **derived**, within the model's rate laws |
| no L1 with sulphide minerals | 0 of 52, few counts | **not established**: needs more runs |

## Limits

- **Abundance filter.** Catalysis is evaluated among sequences with at least 2 copies (top 16 per site). With an alphabet of 2 this means short oligomers. Long, unique sequences can neither catalyse nor be catalytically joined.
- **κ saturates** above about 10⁴: catalysis is limited by the copies available.
- **No uncatalysed joining of two chains**, as in 4c-1.
- **Persistence in turnover units** favours fast-turnover molecules. A version in absolute time is planned.
- **Screening fidelity:** 64 sites × 512 slots.
- **One repetition, COMPASS only.** The L1 boundary is not resolved by the map.
- **Catalysis is assigned at random**, not derived from chemistry. The result is the catalytic ability that polymers would need, for comparison with experiments on random sequences.

## References

- Hordijk, W., Kauffman, S. A. & Steel, M. (2011). Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. *International Journal of Molecular Sciences* 12(5), 3085–3101.
- Hordijk, W. & Steel, M. (2004). Detecting autocatalytic, self-sustaining sets in chemical reaction systems. *Journal of Theoretical Biology* 227, 451–461.
- Full references and checks: [stage4c2_spec.md](../docs/stage4c2_spec.md), [stage4c1_spec.md](../docs/stage4c1_spec.md).
