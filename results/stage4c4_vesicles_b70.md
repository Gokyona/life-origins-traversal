# Stage 4c-4: can heredity live inside vesicles? (Intel Arc Pro B70)

**Dates:** 27 October 2026, 21:25, to 28 October 2026, 23:19 Toronto time · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Specification:** [stage4c4_spec.md](../docs/stage4c4_spec.md), with the implementation choices in its §8 and the focused searches in its §9 · **Validation and run log:** [stage4c4_validation_b70.md](stage4c4_validation_b70.md) · **Code:** chem4c4.py (md5 b79360f6), validate4c4.py, run4c4_full.py, scan4c4mem.py, run4c4_box3.py, run4c4_box3v2.py, confirm4c4.py and confirm4c4v2.py; kept with md5 sums on the benchmark machine and to be added to [stage4/](../stage4/README.md) · **Data:** result JSONs with md5 sums on the benchmark machine, to be added to the repository

**Alphabet:** 2 letters (0/1), as in 4c-3.

## Bottom line

We added membranes to the heredity chemistry of 4c-3:
- amphiphiles from the feedstock;
- a stability window for bilayers (pH near the apparent pKa, limits on Mg²⁺, salt and heat);
- vesicles that enclose nearby chains at random, grow (faster when filled, by osmotic pressure), divide under shear or cycling, and break down outside the window.

Nothing says "be a cell". Then we asked where **persistent vesicles that hold chains (E3)** arise, and whether they can arise **together with heredity (L2)**.

1. **Vesicles and heredity can live together, and then heredity lives inside the vesicles.** At every confirmed point where both occur, essentially all of the inherited lineage is enclosed (86–100%). Copies are still about 60 letters long.
2. **The search found four such environments by itself,** each confirmed with 16 new runs:
   - engine12538, engine11230 and engine15233 in 16 of 16 runs;
   - engine15539 in 13 of 16.

   Random sampling found none: 0 in 1,000 conditions over the whole space, 0 in 500 conditions in each of the two focused regions.
3. **The four environments are cool or mild pools with long temperature cycles** (0–25 °C, cycles of 2.6 to 41 days, slow washout). They are **near-neutral to slightly acidic** (pH 4.4–7.6).
   - Only engine12538 (pH 7.6) lies inside the pH range measured for fatty-acid and alcohol vesicles, which is from about pH 6.5 (Jordan et al. 2019).
   - The other three (pH 4.4–5.9) need a membrane with an apparent pKa near 6, at the lower edge of the model's range, and a pH window up to ±1.75 wide. Whether real amphiphile mixtures allow that is the main open question for these points.
   - The strongly acidic heredity environments of 4c-3 (pH below 4) cannot hold any membrane in the model. The warm, slightly acidic #3709 (39 °C, pH 4.6) leaves only 0.4% of the membrane-constant space compatible, and even then heredity survives enclosure in only 28% of runs.
4. **Enclosure costs little at three of the four membrane-compatible 4c-3 heredity points.** With membranes, the joint rate is 98% at rep1_1926 (7.8 °C, pH 6.7) against 16 of 16 without, 86% at C against 14 of 16, and 64% at #1558 against 9 of 16 at the same seed.
5. **Vesicles alone are rare and do not bring heredity.** In the full 54-dimensional space only 0.2% of conditions give persistent vesicles. There, copying is active but copies are only about 3 letters long.
6. **Vesicles do not change the chemistry below them.** A paired rerun without amphiphiles shows they do not raise polymer formation (E2). They slightly lower it, because they use part of the feedstock.

## What was simulated

The 4c-3 chemistry plus a membrane layer ([specification §1–§3](../docs/stage4c4_spec.md)):
- **Amphiphiles:** made from the amphiphile share of the feedstock. Above a critical concentration (CVC) they assemble into membrane.
- **Stability window:**
  - |pH − pKa_app| ≤ w_pH;
  - Mg²⁺ ≤ Mg_crit;
  - salinity ≤ sal_crit;
  - hottest phase ≤ T_mem.
  Outside the window, vesicles break down.
- **Vesicles:** up to 8 per site. At formation, and at each rehydration after a dry phase, a vesicle encloses each nearby chain with probability equal to its volume share. A vesicle grows by area, and faster when filled, through the osmotic load β (Chen, Roberts & Szostak 2004). It divides under shear or cycling, sharing its contents at random (Zhu & Szostak 2009). Activated monomers cross the membrane; chains of 4 or more letters do not (Mansy & Szostak 2008). Copying happens only between chains in the same compartment.
- **New layer C constants:** 12 (CVC₀, γ_T, pKa_app, w_pH, Mg_crit, sal_crit, T_mem, φ₀, f_rehyd, β, s_shear, s_cycle). The model has 54 dimensions: 13 environmental (layer A) and 41 unknown constants (layer C).
- **Detector:**
  - E3 loose: vesicles persist in at least half of the sites over at least 5 vesicle turnovers.
  - E3 strict: also, they hold chains of length ≥ L_H.
  - The joint level is E3 strict together with L2_long (heredity with copies of at least 6 letters) in the same run.

## Validation

All nine pre-registered tests passed before any map ([details](stage4c4_validation_b70.md)). Without amphiphiles, the model is bit-identical to 4c-3. Encapsulation and division follow the binomial statistics they should (χ² tests). Osmotic competition appears with β > 0 and not with β = 0. Kill and resume are bit-identical.

Three rules were changed after the memory pre-flight and before any validation result, because they made E3 impossible at exactly the heredity points (spec §8, rows 8, 10, 14):
- rehydration is averaged over the cycles in a time step;
- the E3 window is built like the L2 window;
- the run extension was dropped, so the horizon equals that of 4c-3.

## Results

### Full map (54 dimensions)

Ground truth 1000 × 16 runs, COMPASS 32,000 runs, 10.4 h.

| Level | Ground-truth fraction | COMPASS accuracy | 4c-3 |
|---|---|---|---|
| E1 | 0.450 | 0.975 | 0.466 |
| E2 | 0.107 | 0.957 | 0.085 |
| L1 | 0 | – | 0.001 |
| L2_long | 0 | – | 0 |
| E3 strict | 0.002 | – | – |
| E3 strict ∩ L2_long | 0 (0 of 32,000 COMPASS runs) | – | – |

E2 differs from 4c-3 because the 54-dimensional ground truth has different points, not because of vesicles. A paired rerun without amphiphiles gives 0.113 against 0.106 with them.

### Membrane scan at the heredity points

For each 4c-3 heredity point, 64 sets of membrane constants were drawn inside the window, with 4 seeds each.

| Point | Membrane space compatible | Both (E3 strict ∩ L2_long) | Without membranes (4c-3 model, same seed) |
|---|---|---|---|
| #3709 (38.7 °C, pH 4.57) | 0.4% | 0.28 | 16/16 |
| #1558 (7.6 °C, pH 7.96) | 18% | 0.64 | 9/16 |
| C (26.7 °C, pH 4.73) | 1.4% | 0.86 | 14/16 |
| rep1_1926 (7.8 °C, pH 6.73) | 11% | 0.98 | 16/16 |
| A, B, rep1_9823 (pH 2.6–3.8) | 0 | impossible | – |

### Focused maps

| Map | Region | Ground truth: joint | COMPASS: joint candidates | Confirmed (16 runs) |
|---|---|---|---|---|
| box3 | 4c-3 heredity region inside the membrane window (97% masked) | 0 of 8,000 runs | 2 | engine12538 16/16 |
| box3 v2 | the same, with membrane constants drawn inside the window | 0 of 8,000 runs | 5 | engine11230 16/16, engine15233 16/16, engine15539 13/16; engine8083 7/16 (marginal) |

**The confirmed environments:**

| Point | Environment | Both / 16 | Margin | Enclosed share of L2 lineage |
|---|---|---|---|---|
| engine12538 | 4.5 °C (T_m8 13.0 °C), cycles of 63 h, pH 7.58, Mg²⁺ 1.0 mM, 0.74 M salt, τ 11,178 h; amphiphile fraction 0.40, β 9.7 | 16 | 0.57 | 1.0 |
| engine11230 | 7.1 °C (cycling 0–17.1 °C, T_m8 20.1 °C), cycles of 484 h, pH 5.85, Mg²⁺ 38.7 mM, 0.76 M salt, τ 63,326 h; amphiphile fraction 0.26, pKa_app 6.03 ± 0.34, β 1.7 | 16 | 0.36 | 1.0 |
| engine15233 | 25.2 °C (19.2–31.3 °C, T_m8 30.1 °C), cycles of 242 h, pH 4.39, Mg²⁺ 0.35 mM, 0.83 M salt, τ 29,686 h; amphiphile fraction 0.21, pKa_app 6.06 ± 1.75, β 8.6 | 16 | 0.57 | 1.0 |
| engine15539 | 0.0 °C (0–10.0 °C, T_m8 11.9 °C), cycles of 992 h, pH 5.84, Mg²⁺ 0.55 mM, 0.39 M salt, τ 20,978 h; amphiphile fraction 0.33, pKa_app 6.21 ± 0.48, β 7.0 | 13 | 0.43 | 0.86 |

At every confirmed point, E3 strict holds in all 16 runs. The limit is heredity, not the membrane.

## Which conclusions follow directly from the assumptions

| Conclusion | How it arises | Verdict |
|---|---|---|
| Acidic heredity points (pH < 4) cannot hold membranes | The window needs pH ≥ pKa_app − w_pH, and the ranges give pKa_app ≥ 6, w_pH ≤ 2 | **built in** by the layer C ranges. These are more permissive than measured fatty-acid and alcohol mixtures (from about pH 6.5; Jordan et al. 2019), so the conclusion holds for fatty-acid-like membranes. Phospholipid-like membranes are not modelled |
| Heredity ends up inside vesicles | Chains of 4 or more letters cannot leave; copying needs the same compartment; vesicles enclose chains at formation and rehydration | **partly built in**: the enclosure rules are put in. That nearly all lineage mass ends up enclosed, and that heredity survives being split over 8 compartments, is derived |
| Joint points are cool or mild, near-neutral to slightly acidic pools with long cycles | Found by the search and confirmed | **derived**, within the layer C ranges. It combines the membrane window (built in) with the heredity conditions of 4c-3 (derived). Three of the four rely on pKa_app near 6, at the edge of the range and below measured fatty-acid mixtures |
| Enclosure costs heredity at #3709 | Measured against the same point without membranes | **derived**: the main factor is a large amphiphile fraction (logistic fit) |
| Vesicles alone do not bring long copies | Measured in the full map | **derived** |
| Vesicles do not raise E2 | Paired rerun without amphiphiles | **derived** |
| Copies of about 60 letters | They fill the 64-letter cap | **built in**, as in 4c-3 |

## What this layer cannot show

- **The membrane laws are coarse.** A stability window and an exponential CVC stand in for real phase diagrams. Only 0.4–18% of the membrane-constant space suits each heredity point, and three of the four search-found points need pKa_app near 6. engine11230 also has 39 mM Mg²⁺, far above the 4 mM tolerated by the best-studied fatty-acid mixture (Mansy & Szostak 2008) and allowed only through the Mg_crit range, which stands for chelation. So which membrane chemistries are realistic decides much of the answer.
- **Copying does not depend on Mg²⁺, while membranes do.** The known conflict (Adamala & Szostak 2013) is only half modelled.
- **No content-dependent advantage is tested yet.** Whether filled vesicles outgrow empty ones and pass their contents on (L3), and whether that leads to selection (L4), is the subject of [4c-5](../docs/stage4c5_spec.md).

## References

- Adamala, K. & Szostak, J. W. (2013). Nonenzymatic template-directed RNA synthesis inside model protocells. *Science* 342(6162), 1098–1100.
- Chen, I. A., Roberts, R. W. & Szostak, J. W. (2004). The emergence of competition between model protocells. *Science* 305(5689), 1474–1476.
- Jordan, S. F., Rammu, H., Zheludev, I. N., Hartley, A. M., Maréchal, A. & Lane, N. (2019). Promotion of protocell self-assembly from mixed amphiphiles at the origin of life. *Nature Ecology & Evolution* 3, 1705–1714. doi:10.1038/s41559-019-1015-y. (Checked on 28 October 2026 on the author's site and in the UCL accepted manuscript: mixed C10–C15 fatty acid and 1-alkanol vesicles form from about pH 6.5 to above 12.)
- Mansy, S. S. & Szostak, J. W. (2008). Thermostability of model protocell membranes. *PNAS* 105(36), 13351–13355.
- Zhu, T. F. & Szostak, J. W. (2009). Coupled growth and division of model protocell membranes. *JACS* 131(15), 5705–5713.

Checks of the first, second, fourth and fifth references are noted in the [specification](../docs/stage4c4_spec.md).
