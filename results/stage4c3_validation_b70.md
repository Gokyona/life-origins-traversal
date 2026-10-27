# Stage 4c-3 validation: template copying and heredity (L2)

**Status:** validation finished; the full map has not run yet. This page will become the 4c-3 report once the map is done.
**Dates:** 26 October 2026, 16:00 Toronto time, to 27 October 2026, 05:12 · **Hardware:** Intel Arc Pro B70 · **Spec:** [stage4c3_spec.md](../docs/stage4c3_spec.md), with the revisions in its §6
**Alphabet:** 2 letters (0/1), on purpose. The question is whether the simplest alphabet can give rise to heredity, and where. It is not a simulation of modern 4-letter RNA.

Code and raw JSON will be added with the map results, once the last code change (a memory fix, below) is final.

## Bottom line so far

1. **The model can show heredity, and it obeys the error threshold.** Under deliberately ideal copying conditions, strict heredity appears in 15 of 16 runs at an error rate of 1.7% per letter. It disappears completely at 5%. This is an Eigen-type error threshold: about 0.2 → 0.6 errors per copy of a 12-letter master.
2. **Short parasites take over, as known from experiments and theory.** In a long-versus-short contest, short templates end with 99.4% of the lineage mass. Two equal lengths stay near 50:50.
3. **In realistic parameter ranges, heredity was not found.** Two random samples of 4,096 environments each found none, and a targeted scan found none that reproduces. These searches were made before the start-rule correction (spec §6, row 2). A targeted scan under the final rules is still to come.
4. **Under the final rules, the bottleneck is environmental, not chemical.** At the best reference environment of the scan, copying is fast and processive: about 99.999% per letter in the cool phase, and nothing changes when the copy rate is raised 100× or 1000×. But a template makes at most one full copy per temperature cycle. That cycle took 53 hours, while chains lived a few hours at pH 11.4. A sequence was therefore copied in full only about 0.1 times per lifetime, where at least 1 is needed. Heredity in this model needs temperature cycles that are fast compared with molecule lifetime. That is a condition one can compare with real settings.
5. **Several rules as first written suppressed heredity by construction.** Finding and fixing them was most of the work (spec §6). The biggest was the copy start rule: copies began at a random position, so a full copy had probability only 1/L.

## Validation tests

Tests at fixed conditions. They run at validation fidelity, on the final code unless marked.

| # | Test | Result | Numbers |
|---|---|---|---|
| 1 | Copy fidelity: fixed template, no decay | Pass | Measured error 0.0494 vs ε = 0.05 (z = −1.13); 0 primer mismatches; copy-length distribution vs processivity formula χ² p = 0.98; 15k primed and 14k unprimed starts |
| 2 | No copying (k_t = 0) *(earlier code)* | Pass | 0 strict in 16 runs at the box reference; run before the start-rule and statistic revisions, to be repeated with the map's code check |
| 3 | Random-sequence copying as the main run *(earlier code)* | Pass | 0 strict in 16; superseded by test 8, which checks the same null with 128 runs |
| 4 | Error threshold, positive control | Pass, with one deviation | See below |
| 5 | Strand separation | Pass | (a) Below T_m8, no cycling: 99.9% of templates paired, copying stalls (1.4% of tries in the last quarter succeed). (b) Cycling across T_m8: 0.01% paired, 40% succeed. (c) Far above T_m8: copies average 3.0 letters vs 7.8 in (a) |
| 6 | Power cut: kill and resume | Pass | Bit-identical records, selections, life fractions, composition and accuracy curves (after the determinism fix, spec §6 row 8) |
| 7 | Parasite takeover: 24-mer vs 8-mer, resolved cycling, decay on | Pass | Short share 0.994 ± 0.006 (43 sites); equal-length control 0.44 ± 0.13 (16 sites, within 0.5 ± 3 SE); with 2% copy errors the short share is 1.00 |
| 8 | Null false positives: copied letters random | Pass | Strict 1/64 and 0/64 at two points at margin 0.15 |

**Test 4 in detail.** This test needs an environment where heredity exists. No realistic one had been found, so it uses an explicit condition: a seeded 12-letter master, period 0.1 h, hydrolysis 0.01/h, washout τ = 50 h, 100 h run, copy speed-up k_t = 10⁹, and uncatalysed growth slowed to k_c = 10⁻⁶ so that growth does not overwrite lineages. **k_t and the period are far outside the map's ranges.** The test shows only that the detector can see heredity when it exists. It does not show that heredity is plausible.

| ε · L_H | Strict / 16 | Mean margin | Master share | S |
|---|---|---|---|---|
| 0.1 | **15** | 0.21 | 0.367 | 0.59 |
| 0.3 | 0 | 0.04 | 0.180 | 0.41 |
| 1 | 0 | −0.02 | 0.094 | 0.34 |
| 3 | 0 | −0.01 | 0.083 | 0.36 |

- The pre-set criteria were strict ≥ 9/16 at 0.1, 0 at 3, and master share at 3 below 0.2 × the share at 0.1.
- The third criterion missed: 0.083 against a bound of 0.073. The master share counts any 6–12-letter segment of the master or its complement, and with 2 letters many random 6-mers match by chance, so this measure has a floor.
- We count the test as passed on the strict-L2 criteria and report the miss here.

Two earlier designs of this test failed for reasons of test design, and both are recorded:

1. **Normal growth rate:** growth appended about 15 random letters per chain lifetime and erased the master.
2. **A 30 h run:** shorter than the 5 turnover times the detector requires.

## Searches for a natural reference (before the start-rule correction)

| Search | Environments | Copies ≥ 6 letters | Strict L2 | Reproduced in 16 reruns |
|---|---|---|---|---|
| Random over all 42 dimensions | 4,096 | ≈ 0% | 4 | not tested (noise level) |
| Targeted: T just below T_m8, resolved cycling across T_m8, k_t ≥ 10⁴ | 4,096 | 21% | 6 | 0 of 6 |

In the targeted box, copy length and strand separation were mostly solved. The 8-site statistic then in use also had a 15–30% false-positive rate (spec §6, row 6), so these strict counts are at noise level.

## Sensitivity at the box reference (final start rule)

| Setting | Loose / 16 | Strict / 16 | Mean margin |
|---|---|---|---|
| base (k_t = 1.7 × 10⁴) | 16 | 0 | +0.017 |
| k_t = 10⁶ | 16 | 0 | +0.017 |
| k_t = 10⁷ | 16 | 0 | +0.017 |
| primers 2–6 letters | 16 | 2 | −0.025 |
| k_t = 10⁷ and primers 2–6 | 16 | 1 | +0.003 |

- No bound on copying chemistry is binding here.
- The reference has T = 30.7 °C, T_m8 = 35.6 °C, a hot phase of 36.1 °C, a period of 53 h and pH 11.4.
- That gives about 0.1 full copies per template lifetime.

## Still to do

- A targeted scan under the final rules, restricted to environments where chain lifetime ≥ cycle period. The first attempt ran out of GPU memory: the copy step's tensors grow with the number of templates. The fix processes templates in chunks, gives identical results, and is being checked.
- Then the full map: 42 dimensions, ground truth 1000 × 16, COMPASS 32k, about 9 hours.
- Then the report with the assumption table.
