# Stage 4c-1: where do polymers form and persist? (Intel Arc Pro B70)

**Date:** 2026-10-26 · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Map run:** 02:29–05:01 (ground-truth rerun to 07:00) Toronto time · **Specification:** [stage4c1_spec.md](../docs/stage4c1_spec.md) · **Code:** [stage4/chem4c1.py](../stage4/chem4c1.py), [stage4/run4c1_map.py](../stage4/run4c1_map.py) · **Data:** [stage4/map4c1_result.json](../stage4/map4c1_result.json) (published version), [stage4/sens4c1_result.json](../stage4/sens4c1_result.json) (sensitivity)

## Bottom line

We simulated polymer chemistry from feedstock upward, with nothing living put in, across the full stage 4 environment space.
- **Activated monomers (E1)** reach a working concentration in 44% of that space.
- **Persistent chains of 8 or more units (E2)** form in 16%.
- **What decides whether chains persist:**
  - the residence time of the water;
  - temperature (cold is better);
  - feedstock supply;
  - an energy source (UV, or H₂ from below).
- Mineral surfaces and salinity hardly matter at this layer. Wet–dry cycling helps at every temperature but does not rescue hot water.
- **The chains' composition mostly follows the feedstock.** Hot environments shift it towards peptide-like chains, because peptide bonds resist heat better.
- **COMPASS needs about 2× fewer simulations than Sobol sampling** for the E1 map and about 1.6× for E2. That is less than in the test models, because here a large part of the space is plainly negative.

## What was simulated

This is the first layer of the bottom-up simulator. Nothing living goes in. Each run is one environment: a point in 29 dimensions.
- 19 dimensions encode the 13 environmental variables of layer A, with ranges from [stage4a_ranges.md](../docs/stage4a_ranges.md).
- 10 are constants that no measurement fixes yet (layer C): the condensation rate, the dry-phase boost, the water activity when dry, the UV and H₂ yields, the deactivation energy, mineral adsorption and protection, the salinity exponent of hydrolysis, and UV bond damage.

One run has 64 reaction sites with 512 chain slots each. This is screening fidelity, at 45 ms per run; full fidelity has 256 × 2,048 and costs 559 ms per run.

In each site:
- feedstock becomes monomers;
- monomers are activated, using UV or H₂ energy, and drying;
- activated monomers start and lengthen chains;
- chains hydrolyse, wash out, adsorb to minerals and break under UV.

Hydrolysis follows Li & Breaker (1999) for nucleotide-like chains and Radzicka & Wolfenden (1996) for peptide-like chains. Wet–dry and hot–cold cycling, liquid-water limits and the available time are modelled as in the specification.

Two levels are scored:
- **E1:** activated monomers reach a working concentration (0.01 mM effective) at some time during the run.
- **E2:** at least half the sites hold a chain of at least 8 units, averaged over the last quarter of the run. Lengths 4 and 16 are also recorded.

The map uses E2 together with E1, so the two levels nest. Without that rule, E2 occurs without E1 in 0.2% of runs: these are conditions where slow hydrolysis lets chains accumulate from a supply below the threshold.

## Validation (all passed before the map)

| test | criterion | result |
|---|---|---|
| (a) steady state | chain-length distribution matches the exact steady state of the same kinetics | total-variation distance 0.0009 / 0.0028 / 0.0026 at growth/break ratios 0.5 / 1 / 1.5; mean length 2.251 / 2.498 / 2.755 against 2.250 / 2.500 / 2.750 |
| (b) hydrolysis law | slope 10^0.07 per °C; reference rate 1.4 × 10⁻⁷ /min | slope 0.07014 (0.2% off); 1.390 × 10⁻⁷ /min at 37 °C; peptide half-life 401 years at 25 °C, 5.54 weeks at 100 °C |
| (c) cycling | wet–dry cycling lengthens chains at the same mean water activity | at 80 °C mean length 3.12 against 2.12, and sites with chains ≥ 8: 0.78 against 0; at 90 °C chains only with cycling |
| (d) washout, no feed | τ = 1 h kills E2; no feedstock kills E1 | 4 of 4 each |
| (e) power cut | kill −9 and resume gives the same result | simulations, selections and accuracy curves bit-identical |
| inactive variable | amphiphile fraction has no effect in this layer | output bit-identical |

*Correction to the specification:* test (a) first named the Flory distribution. That distribution belongs to reversible step-growth, where any two chains can join. Here chains grow by adding monomers only, so the test compares with the exact steady state of this kinetics. Such a steady state exists only when growth is less than twice the per-bond break rate.

## Results

**Life fractions** (ground truth, 2,000 environments × 32 simulations; an environment counts as positive when more than half its runs are):

| level | share of environments |
|---|---|
| E1 (activated monomers at working concentration) | 0.4385 |
| E2 (chains ≥ 8 persist; E1 and E2) | 0.1575 |

**Quality flags:**
- Liquid water in 95.9% of runs.
- 99.4% reached the planned end of the run.
- 87.6% were stationary.
- 4.3% had at least one fragment lost to a full slot table.

**Accuracy of the map along the budget** (share of test environments classified correctly; COMPASS / Sobol, one repetition each):

| level | 2k | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| E1, all | .892 / .888 | .927 / .932 | .955 / .938 | .963 / .947 | .965 / .964 | **.968** / .962 |
| E1, hard (140) | .643 / .657 | .721 / .743 | .779 / .729 | .807 / .714 | .814 / .807 | **.814** / .800 |
| E2, all | .934 / .935 | .956 / .942 | .969 / .968 | .976 / .968 | .977 / .972 | **.980** / .978 |

**Savings:**
- COMPASS reaches Sobol's 64k accuracy at about 15k simulations for E1 and 40k for E2.
- For E1, Sobol's own curve is not monotone: it is already at that level by about 30k. The honest saving is therefore about 2× for E1 and 1.6× for E2.
- The E2 hard set has only 11 environments, too few to use.

Both methods spent about 42–44 minutes each. Simulation dominates, so the engine's overhead is small.

**What matters for E2** (drop in map accuracy when one input is scrambled):

| input | importance | | input | importance |
|---|---|---|---|---|
| residence time | 0.146 | | H₂ | 0.039 |
| temperature | 0.115 | | pH | 0.032 |
| feedstock supply | 0.080 | | H₂ yield (C) | 0.021 |
| UV flux | 0.059 | | cycling depth | 0.016 |
| UV yield (C) | 0.053 | | dry-phase boost (C) | 0.011 |

All other inputs are at or below 0.011. Mineral availability and type, adsorption, mineral protection, salinity and dry-phase water activity have no visible effect: their partial-dependence range is at most 0.002. The amphiphile fraction acts only in a later layer; its importance is 0.000 for E2, as required.

**Partial dependence of E2** (probability of E2 when the input is set to the given position in its range, averaged over everything else):

| input (range) | 0 | ¼ | ½ | ¾ | 1 |
|---|---|---|---|---|---|
| residence time (1 h – 10⁵ h) | .000 | .016 | .125 | .304 | .465 |
| temperature (0–150 °C) | .354 | .266 | .152 | .057 | .013 |
| feedstock supply | .018 | .082 | .157 | .247 | .335 |
| UV flux | .085 | .108 | .151 | .208 | .268 |
| H₂ | .118 | .122 | .140 | .191 | .279 |
| pH (2–12) | .193 | .189 | .176 | .140 | .072 |
| cycling depth | .130 | .146 | .161 | .174 | .186 |
| Mg²⁺ | .175 | .170 | .162 | .154 | .143 |

**Two slices for E2:**
- **Temperature × cycling.** Without cycling: .33, .22, .11, .04 and .01 from cold to hot. With full cycling: .39, .31, .19, .08 and .02. Cycling helps everywhere, but not enough to make hot water work.
- **pH × Mg²⁺.** At low Mg²⁺: .21, .21, .20, .18 and .12 from acid to alkaline. At high Mg²⁺: .20, .19, .18, .13 and .05. Alkaline, Mg-rich water is worst, as the hydrolysis law predicts.

**Composition of the chains where E2 holds** (10,059 positive runs):
- **Chemistry.** The nucleotide-like share of chain mass has a median of 0.55 (10th–90th percentile 0.05–0.96). 18% of runs are mostly nucleotide-like (> 0.9) and 14.5% mostly peptide-like (< 0.1).
- **Feedstock.** The share follows the feedstock: medians 0.11, 0.52 and 0.90 for the three thirds of the nucleotide fraction in the feed.
- **Temperature.** The nucleotide-like share falls in hot water: 0.57, 0.56 and 0.37 for 0–50, 50–100 and 100–150 °C. Mean chain length falls with it: 15.2, 9.1 and 6.9.
- **Minerals.** Sulphide minerals lengthen chains slightly: mean length 13.5, against 10.8 without minerals and 10.9 with clay.
- **Sequence diversity.** Chains of 8–15 units are nearly all different sequences: median 673 distinct sequences, 9.3 bits. That is expected with no selection at this layer.
- **Length cap.** 17.6% of chains reached the length cap of 32, so the cap binds. It is raised to 64 from layer 4c-2.

**Power-cut test, applied for real.** At 03:01 one batch stalled. On this GPU, Poisson draws with rates above about 4 × 10⁹ return a garbage value, and very concentrated environments reach such rates. The fix (exact draws up to rate 10⁴, a normal approximation above, rates capped at 10⁷) changes nothing else. The run resumed from its checkpoint. The first 1,056 ground-truth environments had used the old code, so the whole ground truth was rerun with the fixed code:
- E1 changed at 1 of those 1,056 environments and E2 at none;
- the other 944 environments were bit-identical.

The published numbers are from that rerun. The GPU stayed at 68–78 °C, with no thermal pause.

## Which conclusions follow directly from the assumptions

"Bottom-up" means nothing living is put in, but the chemical rules are. Some conclusions are therefore fixed by those rules, and the simulation only gives them a size. We sort them into three groups:
- **built in:** the direction was fixed when the rules were written;
- **partly built in:** the mechanism is derived, but its rank or strength depends on the chosen ranges;
- **derived:** no rule states it; it appears only when several processes interact.

| conclusion | assumption it mainly comes from | verdict |
|---|---|---|
| cold water is better | hydrolysis rises steeply with temperature; condensation has no temperature dependence in this layer | built in |
| alkaline water is worse, acid water harmless | the Li & Breaker law covers base catalysis only | built in (the acid side is untested) |
| composition follows the feedstock | two independent chemistries with the same condensation rate | built in |
| minerals hardly matter | minerals act only through adsorption and protection here, and adsorbed fractions are mostly small in the chosen ranges | built in for this layer |
| salinity has no effect | salinity enters only through a hydrolysis exponent spread symmetrically over −1 to 1 | built in |
| available time hardly matters | time compression: pools start in their state after the available time | built in |
| an energy source is required | monomer synthesis and activation are proportional to the energy terms | built in |
| the shares 44% (E1) and 16% (E2) | ranges, log-uniform sampling and detector thresholds | built in as absolute numbers |
| residence time ranks first | washout competing with growth is derived; the rank depends on its five-decade range | partly built in |
| hot water favours peptide-like chains | the two measured hydrolysis laws | partly built in |
| cycling helps but does not rescue hot water | lower water activity and faster growth when dry are built in; the size of the effect is derived | partly built in |
| saturation of thermophoretic concentration at about 100× | interaction with other rate limits | derived |
| where the E2 boundary lies in residence time × temperature × feedstock × energy | competition of growth, hydrolysis and washout | derived |

The shares of the space should be read only within the model. The layer's main result is where the boundary lies: no single rule shows that. Layer 4c-2 adds an activation energy for condensation and mineral catalysis, to test the first and fourth rows.

## Limits

- **Hydrolysis outside the measured range.** The Li & Breaker law covers base catalysis only, and Mg²⁺ from 5 to 50 mM. In acid water a different, acid-catalysed path, which the law leaves out, can dominate; below 5 mM Mg²⁺ the law's rate falls towards zero. The map may therefore make nucleotide chains look too stable in acid or Mg-poor water. A sensitivity run put floors on both terms: pH at 5 and the Mg²⁺ term at 5 mM. It used the same 2,000 environments with 8 paired seeds each.
- E2 changed at 1.15% of environments, and the E2 share went from 0.157 to 0.152.
- The E2 rate at pH 2–5 did not change (0.178 against 0.181).
- In the lowest third of Mg²⁺ it fell from 0.160 to 0.148.
The map is robust to both extrapolations.
- **Screening fidelity.** 64 sites instead of 256; points near the boundary are to be rerun at full fidelity.
- **Time.** Available time enters only through slow processes; see the specification, Section 3.
- **One repetition** of each method.
- **No temperature dependence of condensation** in this layer: temperature acts through hydrolysis, deactivation and water's boiling point only. Layer 4c-2 adds an activation energy for condensation.
- **Slot overflow** in 4.3% of runs: fragments that found no free slot were returned to the monomer pool.

## References

- Li, Y. & Breaker, R. R. (1999). *JACS* 121(23), 5364–5372.
- Radzicka, A. & Wolfenden, R. (1996). *JACS* 118, 6105–6109.
- Full references and checks: [stage4c1_spec.md](../docs/stage4c1_spec.md).
