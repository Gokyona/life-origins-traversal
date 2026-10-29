# Stage 4c-5: coupled protocells, and why Darwinian selection stops (Intel Arc Pro B70)

**Dates:** 28 October 2026, 22:40, to 29 October 2026, 05:29 Toronto time · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Specification:** [stage4c5_spec.md](../docs/stage4c5_spec.md), with its revisions, test records and implementation choices in §6 · **Validation:** [stage4c5_validation_b70.md](stage4c5_validation_b70.md) · **Code:** chem4c5.py (md5 c747dc68), validate4c5.py, validate4c5_6b.py, validate4c5_6c.py, run4c5_5a.py; kept with md5 sums on the benchmark machine and to be added to [stage4/](../stage4/README.md) · **Data:** map4c5_5a result.json (md5 c5013592) and runs.jsonl (ea2c8c4b)

**Status:**
- **Stopped after the point study, by the owner's decision.**
- **Two pre-registered validation tests failed (3 and 6).** The point study ran as a recorded deviation from the stop rule (see Validation).
- The focused map of the plan (5b) was not run. The extent of the coupled-protocell region and its drivers are left to the joint map of stage 4d.

**Alphabet:** 2 letters (0/1), as before.

## Bottom line

4c-5 adds no chemistry. It records what the vesicles of 4c-4 do and looks for two levels:
- **L3, coupled protocell:** vesicles whose contents help them grow, and which pass those contents on by copying;
- **L4, Darwinian evolution:** a heritable vesicle property rises beyond drift.

1. **Coupled protocells arise wherever heredity and vesicles coexist.** At all 10 known coexistence environments, the point study found L3 strict: content-dependent growth together with inheritance by copying over many vesicle generations. That is 298 of 435 runs in which heredity and vesicles coexisted (69%), with 16 of 16 at engine12538, engine15233 and engine15539. The paired random-copy control never shows inheritance (0 of 608).
2. **The coupling is physical, not informational.** Filled vesicles grow faster because of osmotic pressure. Without it (β = 0), coupling vanishes in all 608 control runs. With random-sequence copying, growth coupling stays (515 of 608), because the load counts letters whatever their sequence, but nothing is inherited. Heredity is what turns a filled vesicle into a coupled protocell.
3. **Darwinian selection among vesicles does not continue in this model** (L4 = 0 in all 435 runs).
   - Selection does happen, as one fast sweep. In a test with copying and non-copying founder vesicles, the copying families leave 100% of the descendants (test 6b, 16 of 16 against a drift null).
   - After the sweep, copying refills every surviving vesicle to the same load within a few hundred hours. No heritable variation in content remains for selection to act on (test 6c calibration).
   - At the real points, vesicles with more lineage per volume even leave slightly fewer descendants (ρ(w, z) from −0.07 to −0.32).
   - Sustained evolution would need a heritable trait that copying does not erase. Chains in this model have none: sequence does not change copy speed or function, and every chain grows to the 64-letter cap.
4. **The bottleneck of the ladder is therefore not L3 but the step before it**: the rare coexistence of heredity and stable membranes found in 4c-4. Beyond L3, the limit is the model itself, which has no sequence-dependent phenotype.

## What was detected

- **L3 loose:** the relative growth of each vesicle between checkpoints correlates with its osmotic load (Spearman ρ ≥ 0.1). The test uses a permutation null that shuffles contents among the vesicles of the same site and checkpoint (p < 0.01).
- **L3 strict:** L3 loose, plus inheritance by copying.
  - A descendant must hold at least 0.5 times its ancestor's dominant lineage letters in chains copied after the ancestor's last checkpoint, with at least 90% sequence identity.
  - This must hold in at least half of the sites.
  - The share of sites passing must exceed that of a paired random-copy run by at least 0.25.
  - The ancestor used was on average about 12 generations up (spec §6, choice 13), which is stricter than the specified 3.
- **L4:** L3 strict, plus heritability (a positive parent–daughter regression of z) and selection beyond drift. Selection is measured as the correlation of an ancestor's descendant count with its lineage load z, against a permutation (neutral-marker) null.
- **Paired controls at every point:** β = 0 (no osmotic coupling) and random-sequence copying. Both are scored without the coexistence gate.

## Validation

[Details](stage4c5_validation_b70.md).
- **Before any GPU result,** CPU tests and diagnostics revealed problems that led to three revisions of the inheritance criterion and one fix of the test fixtures:
  - re-enclosure of released chains on rehydration;
  - a closed fixture with full chain slots;
  - ancestors so far back that copy errors had changed every sequence.
- **Passed:** memory pre-flight, bit-identity with 4c-4, record consistency, inheritance (16 of 16 with copying against 0 of 16 without), false positives (0 of 80), and kill and resume.
- **Test 3 failed** on its pre-registered effect size. Coupling was detected in 16 of 16 runs, but the median ρ was 0.27 against a guessed threshold of 0.3.
- **Test 6 failed** (the L4 positive control, 3 of 16). A bit-identical replay showed that the observation window opened only after selection had finished. Test 6b, written before it ran, observed from t = 0 and passed (16 of 16, 1 of 16 false positives).
- **Test 6c** was meant to test the power of the map's L4 statistic. Its CPU calibration showed that differences in content load do not persist (above), so under its pre-registered rule it was not run.
- **By the rules of §3, the two failures meant no maps.** The point study ran on the coordinator's decision as a recorded deviation, after test 6b had passed. The deviation is written into its output files.

## Point study (5a)

16 seeds per point and membrane set, screen fidelity, 608 runs. The 9a points use the first 8 window-compatible membrane sets each (choice 17). "Gated" means heredity and filled vesicles coexist in the run (L2_long and E3 strict).

| Point | Gated | L3 loose | **L3 strict** | L4 | Site share: main / random copy | ρ(g, Π) | ρ(w, z) | Generation time |
|---|---|---|---|---|---|---|---|---|
| engine12538 | 16/16 | 16 | **16** | 0 | 0.63 / 0.03 | 0.51 | −0.32 | about 400–470 h |
| #1558 (8 sets) | 111/128 | 75 | **74** | 0 | 0.72 / 0.03 | 0.43 | −0.22 | 1.2–2.6 kh |
| #3709 (8 sets) | 39/128 | 39 | **38** | 0 | 0.84 / 0.00 | 0.43 | −0.14 | 209 h |
| C (8 sets) | 78/128 | 78 | **44** | 0 | 0.59 / 0.03 | 0.41 | −0.25 | 1.2–2.6 kh |
| rep1_1926 (8 sets) | 128/128 | 80 | **75** | 0 | 0.70 / 0.02 | 0.42 | −0.25 | about 400–470 h |
| engine4527 | 10/16 | 10 | **10** | 0 | 0.83 / 0.02 | 0.49 | −0.22 | 1.2–2.6 kh |
| engine8083 | 8/16 | 6 | **6** | 0 | 0.61 / 0.02 | 0.32 | −0.27 | 1.2–2.6 kh |
| engine11230 | 13/16 | 13 | **3** | 0 | 0.39 / 0.02 | 0.47 | −0.07 | 1.2–2.6 kh |
| engine15233 | 16/16 | 16 | **16** | 0 | 0.69 / 0.02 | 0.49 | −0.32 | 1.2–2.6 kh |
| engine15539 | 16/16 | 16 | **16** | 0 | 0.67 / 0.00 | 0.44 | −0.23 | 1.2–2.6 kh |
| **Total** | **435/608** | **349** | **298** | **0** | | | | |

| Control (all 608 runs, scored without the gate) | L3 loose | L3 strict | L4 |
|---|---|---|---|
| β = 0 | 0 | 0 | 0 |
| Random-sequence copying | 515 | 0 | 0 |

- **Typical coupled protocell:** 3–12 enclosed chains of at least L_H letters, all at the 64-letter cap. It divides at 1.8–4.1 times the membrane of a newborn vesicle and runs through about 290–300 vesicle generations in a run.
- At #3709 and C, the coexistence gate depends strongly on the membrane set (#3709: 0–15 of 16 per set), as found in 4c-4.
- engine11230 holds heredity inside vesicles but passes the inheritance test in only 3 of 13 runs (site share 0.39, below 0.5). Its lineage keys change faster, consistent with its higher Mg²⁺ and lower β; this is not tested.
- The vesicle-id table overflowed in 27 of 608 runs (flagged per run). This can only make L3 strict harder to reach.

## Which conclusions follow directly from the assumptions

| Conclusion | How it arises | Verdict |
|---|---|---|
| Content helps growth (L3 loose) | Osmotic growth with β > 0 is put in (4c-4) | **built in**: the mechanism is a rule. That it acts at every coexistence point and vanishes with β = 0 is a check, not a discovery |
| Contents are inherited by copying | Measured on chains copied after the ancestor, against a random-copy control | **derived**: the model has no rule that vesicle families keep lineages; it follows from copying inside compartments and random sharing at division |
| L3 follows coexistence almost always | Counted in the point study | **derived**, at 10 points only; the extent of the region is not mapped |
| No sustained Darwinian selection | Content load is refilled by copying; chains have no other trait | **built in by omission**: the model gives sequences no phenotype, so no heritable variation can persist. A finding about the model, not about environments |
| Selection by sweep of copying families | Measured in test 6b | **partly built in**: it follows from β and copying, both put in; its speed is derived |

## What this step cannot show

- **How large the L3 region is, and which variables drive it.** 5b was not run; this goes to the 4d joint map.
- **Evolution beyond the first sweep.** It needs heritable variation that copying does not erase, for example sequence-dependent copying, catalysis by chains, or membrane-binding chains. Each would be a new chemistry layer with its own specification and validation.
- **The layer C caveats of 4c-4 remain.** Three of the four search-found coexistence points need membranes with pKa_app near 6, at the edge of the modelled range and below measured fatty-acid mixtures.

## References

As in [4c-4](stage4c4_vesicles_b70.md): Chen, Roberts & Szostak 2004 for osmotic competition, and Zhu & Szostak 2009 for growth and division. Stages 3a and 3c are in this repository.
