# Part III, P3-1: why the model misses the cold benefit seen by Aitken et al. (2023)

**Status:** analysis, 30 October 2026, Toronto time (CPU only, on the B360; report `lot-p3/p31/cold_gap_report.md`, script `cold_gap.py`, numbers `cold_gap.json`, to be uploaded with the P3-1 data) · **Spec:** [part3_p31_spec.md](part3_p31_spec.md), validation test 6

The owner asked for this gap to be understood before anything else.

## 1. The experiment, exactly

Source: Aitken, Wright, Radakovic & Szostak (2023), *JACS* 145, 16142–16149, doi:10.1021/jacs.3c04635, main text and SI (University of Chicago repository; MD5 sums match the repository's).

| Metric | Room temperature | −3 °C | Ratio | Protocol |
|---|---|---|---|---|
| Primers extended by ≥1 letter at 24 h | 65% | 87% | 1.34 | CMP and GMP 10 mM each, 30 mM Mg²⁺, pH 8, activating reagents (200 mM) added once |
| Primers extended by ≥7 letters after 5 days | 7% | 60% | 8.6 | as above plus six helper trimers; reagents re-added every 24 h |
| Pre-activated control, ≥7 letters at 24 h | 42% | 8% (+3 °C) | 0.19 | no in-place activation |

- The −3 °C solutions were **liquid** ("no freezing or eutectic"), so ice concentration plays no part.
- 1.34⁷ ≈ 7.8, close to 8.6: the large 5-day effect is a per-step gain of about 1.3 compounded over seven letters.
- With pre-activated building blocks, **cold slows copying**. The cold benefit therefore comes from the activation chemistry, not from the copying step.
- The bridged A\*A intermediate peaks at about 2.9 mM after 80 min at room temperature and falls to 1.5 mM by 400 min; at −3 °C it rises slowly to 3.4 mM and stays there.
- The authors' explanation: side reactions that waste the activating reagents (AppA and Passerini products) are suppressed in the cold, and A\*A hydrolyses more slowly.

## 2. Model against experiment, same protocol and metrics

Measured A\*A curves were used per reagent addition; one rate scale was fitted at room temperature; −3 °C was predicted.

| Model variant | ≥+1 at 24 h (observed 0.87) | ≥+7 at 5 days (observed 0.60) |
|---|---|---|
| Extension temperature term only (E_ext = 18.4 kJ/mol, from the pre-activated control) | 0.41 | 0.002 |
| Cold A\*A formation only | 0.74 | 0.17 |
| Cold A\*A persistence only | 0.72 | 0.17 |
| Cold A\*A + E_ext 18.4 | 0.54–0.56 | 0.02–0.04 |
| Cold A\*A + E_ext 0 | 0.78–0.80 | 0.31–0.41 |
| Cold A\*A + E_ext −10 | 0.89–0.91 | 0.69–0.78 |
| Binding that tightens in the cold (ΔH_bind −50 to −100 kJ/mol), E_ext 18.4 | 0.63–0.65 | 0.06–0.07 |

## 3. Diagnosis

1. **The model has no activating-reagent chemistry.** In P3-1 the activating agent is a constant supply a. In the experiment the reagent is consumed, mostly by side reactions, in about 2 h at room temperature but lasts more than 7 h at −3 °C. This is the main structural gap.
2. **Activated-block lifetime:** its temperature dependence explains about half of the 24 h gain (apparent activation energy of A\*A loss ≥ 21 kJ/mol, probably near 50).
3. **What is left:** about a 2.8-fold gain per step is not covered by any existing model term. It matches only if the in-place extension had an apparent activation energy near −9 kJ/mol. The pre-activated control shows the real copying step has about +18. So the missing factor is not in copying; it must come from the reagent chemistry (less waste, fewer inhibiting side products), which item 1 adds.
4. The model's 0 °C floor (a Stage 4 rule) also prevented any test at −3 °C.

No direct measurement of the binding enthalpy of bridged dinucleotides or of the activation energy of their hydrolysis was found in two searches.

## 4. Proposed change (not yet applied)

- **(a) Activating-reagent pool:** supply rate or periodic additions, temperature-dependent loss k_r(T), and a temperature-dependent share y(T) of reagent spent on side products. Calibrated on the Aitken time courses (room temperature and −3 °C), labelled as a calibration on one experiment.
- **(b)** Activation energy of activated-block hydrolysis: range 20–100 kJ/mol (was 50–100).
- **(c)** Extension activation energy kept at the measured band 0–20 kJ/mol. No negative values: the cold gain must come from (a), not from a fudge factor in copying. If (a) cannot close the gap, the remaining factor is reported as unexplained.
- **(d)** Temperature floor lowered to −10 °C; liquid only above the freezing point of the solution (about −3.7 °C per mole of dissolved salt per litre).

**Acceptance test (replaces test 6):** with (a)–(d), the model must reproduce both metrics of Table 1 within their uncertainty, and the pre-activated control (cold slows copying), using one parameter set.
