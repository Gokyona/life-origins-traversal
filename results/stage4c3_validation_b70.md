# Stage 4c-3 validation: template copying and heredity (L2)

**Status:** validation finished; first natural heredity found; full-space map done (27 October, 14:19 Toronto); focused map of the heredity region done; 5 confirmed heredity environments (27 October, evening). This page will become the 4c-3 report.
**Dates:** 26 October 2026, 16:00 Toronto time, to 27 October 2026, 05:12 · **Hardware:** Intel Arc Pro B70 · **Spec:** [stage4c3_spec.md](../docs/stage4c3_spec.md), with the revisions in its §6
**Alphabet:** 2 letters (0/1), on purpose. The question is whether the simplest alphabet can give rise to heredity, and where. It is not a simulation of modern 4-letter RNA.

Code: [stage4/](../stage4/README.md) (md5-checked copies of the versions that ran). Raw JSON will be added later.

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

## Update, 27 October 2026: natural heredity under the final rules

**Targeted scan (box2).** 4,096 environments. Each has T just below T_m8, a hot phase above it, chain lifetime ≥ cycle period, and k_t ≥ 10⁴. 204 were loose and 10 strict at screening fidelity. Eight candidates were rerun 16 times at validation fidelity. Three gave strict heredity in 16 of 16 runs, with copies of about 60 letters and margins of 0.4–0.6:

| Point | pH | T / T_m8 (°C) | Period (h) | ε (per letter) | Long copies per template per lifetime (R) | Status |
|---|---|---|---|---|---|---|
| #1558 | 7.96 | 7.6 / 10.2 | 118 | 0.0036 | 0.98 | **stands**: a cold, near-neutral pool |
| #3709 | 4.57 | 38.7 / 46.3 | 33 | 0.0063 | 1.05 | **stands**: unchanged with a pH floor of 5 for hydrolysis (16/16, margin 0.619) |
| #813 | 2.68 | 54.2 / 69 | 378 | – | 0.42 | **artifact**: 0/16 with a pH floor of 5 |

- The hydrolysis law (Li & Breaker) has no acid branch, so lifetimes at pH < 5 may be too long. #813 owes its heredity to that gap and is dropped.
- The reproduction number R separates the classes: 0.4–1.05 at the reproducible points, ≤ 0.05 at points without heredity.

**Copy-length cap.** At #3709, 99.97% of chains sit at the model's length cap of 64 letters, which is a hard limit of the sequence encoding. At #813 the share is 6%. A check with caps of 48 and 32 is queued: it tells whether the cap sets only lineage length or also whether heredity happens.

**Test 9: error threshold at a natural reference.** S is the heredity statistic.

| Point | Reference | ε·L_H = 0.1 | 0.3 | 1 | 3 | Pre-registered verdict |
|---|---|---|---|---|---|---|
| #1558 | 16/16, S 0.42 | 0/16, S 0.069 | S 0.002 | S 0 | S 0 | **fail** |
| #3709 | 15/16, S 0.57 | 2/16, S 0.14 | S 0.004 | S 0 | S 0 | pass |

- #1558 fails the criterion "rate at 0.1 > rate at 3" (0 vs 0) because heredity is already gone at the first grid point.
- *Post-hoc reading, not a verdict:* the relevant length is the copy length (about 62 letters), not L_H = 6. ε·L_copy is 0.22 and 0.39 at the references, where heredity survives, and about 1 at ε·L_H = 0.1, where it is gone. That is Eigen's threshold ε·L ≈ 1 at both points.
- A follow-up with a finer grid (ε·L_copy = 0.15–3) and a criterion fixed before running is queued (test 9b).

**Full-space map (map4c3).** 42 dimensions; ground truth 1000 × 16 runs; COMPASS 32,000 runs; 5.6 h.

| Level | Ground-truth fraction | COMPASS accuracy | Top variables |
|---|---|---|---|
| E1 | 0.466 | 0.961 | residence time τ, T, Y_UV (C), UV, thermal cycling, H₂, feedstock |
| E2 | 0.085 | 0.959 | T, τ, feedstock, Y_UV (C), UV, pH |
| L1 | 0.001 | – | – |
| L2 (strict) | 0.000 | – | – |

- COMPASS found 77 strict-L2 runs, but **all are short-copy noise**: median copy length 2.4 letters, none ≥ 6. In the full 42-dimensional space, real heredity is too rare for 32,000 samples to hit; it lives in the narrow box2 region. The L1 and L2 maps are therefore uninformative ("no" everywhere).
- **Detector consequence.** A new level, L2_long = strict L2 and mean copy length ≥ L_H, keeps noise-level copying from counting as heredity. A focused map of the box2 region with L2_long (ground truth 500 × 16, COMPASS 16,000, infeasible environments masked) is running.

## Update, 27 October 2026 (evening): focused map, confirmations, cap check, test 9b

**Focused map of the box2 region (map4c3_box2).** Level L2_long = strict L2 and mean copy length ≥ L_H. Environments with chain lifetime < cycle period are masked as infeasible and not simulated. The run had ground truth of 500 × 16 runs, COMPASS of 16,000 runs, and passed a kill-and-resume test and a memory pre-flight first.
- Ground truth: E1 0.452, E2 0.184, L1 0.002, L2_long 0. 46% of the region is infeasible.
- COMPASS accuracy at 16,000: E1 0.968, E2 0.95. The L1 and L2 maps stay uninformative because positives are too rare. Importance is led by residence time τ, pH and cycle period.
- COMPASS found 3 L2_long runs. Ground truth found none, consistent with about 3 long-copy points per 4,096 in the box.

**Confirmation of the 3 new points.** 16 reruns each at validation fidelity, with both copy_random controls, plus a hydrolysis pH floor of 5:

| Point | T / T_m8 (°C) | Period | pH | k_t, ε, L_t | Strict / 16 (as is → pH floor 5) | Margin | Copy length | R |
|---|---|---|---|---|---|---|---|---|
| A | 14.2 / 26.2 | 1.14 h | 3.83 | 5.2e4, 0.0017, 10 | 16 → 16 | 0.54 → 0.55 | 62.5 | 0.93 |
| B | 21.4 / 25.2 | 3.8 h | 2.64 | 1.1e4, 0.0074, 12 | 12 → 12 | 0.33 → 0.30 | 62.4 | 0.98 |
| C | 26.7 / 37.2 | 633 h | 4.73 | 2.5e4, 0.0013, 5 | 16 → 16 | 0.61 → 0.56 | 62.4 | 0.96 |

- **Why the pH floor barely matters here.** It raises hydrolysis per bond by ×14, ×209 and ×1.8 at A, B and C. Even so, hydrolysis of a 64-mer is only 2e-4 to 8e-4 of the washout rate, because the water is cool or mild. Chain lifetime is set by washout, and the conclusion holds even if hydrolysis were about 100× faster.
- At #813, by contrast (54 °C), hydrolysis at pH 5 is about 30% of all loss, which explains its collapse.
- Caveat: pH 2.6 (B) is still an extrapolation of the hydrolysis law.

**Heredity points so far: 5.** #1558, #3709, A, B and C. All have:
- cycling across the strand melting temperature;
- chain lifetime set by slow washout (τ = 2×10⁴ – 8×10⁴ h);
- copy errors below about 0.3 per full copy;
- R ≈ 1.

Cycle periods range from 1 h to 26 days, and temperatures from 8 to 39 °C.

**Copy-length cap.** The model's cap of 64 letters is a hard limit of the sequence encoding. With caps of 48 and 32, heredity at #1558 and #3709 remains (15–16 of 16 strict), copy length follows the cap, and the margin rises (#3709: 0.61 → 0.74 → 0.86). The cap sets lineage length, not whether heredity happens. Chains sit at the cap in every arm.

**Test 9b: error threshold, criterion fixed before running — PASS at both points.** x = ε·L_copy, with L_copy = 62.

| x | #1558 strict / 16 | #1558 S | #3709 strict / 16 | #3709 S |
|---|---|---|---|---|
| 0.15 | 16 | 0.542 | 16 | 0.873 |
| 0.3 | 16 | 0.381 | 16 | 0.697 |
| 0.5 | 5 | 0.243 | 16 | 0.497 |
| 0.7 | 1 | 0.156 | 16 | 0.329 |
| 1.0 | 0 | 0.096 | 1 | 0.157 |
| 1.5 | 0 | 0.030 | 0 | 0.050 |
| 3 | 0 | 0.003 | 0 | 0.005 |

- The criteria: S non-increasing, strict ≥ 12/16 at x ≤ 0.3, and strict ≤ 2/16 at x ≥ 1.5.
- Heredity is lost between about 0.5 and 1 errors per copy, as Eigen's error threshold predicts.
- The original test 9 (grid in ε·L_H) stays recorded as a fail at #1558.

**Caveat: saturated regime.** At the heredity points nearly all copy slots hold chains at the length cap, and R ≈ 1 may be limited by the slot count. Heredity here is maintenance of a saturated population, not growth.

## Still to do

- The 4c-3 report with the assumption table, written as this page's final form.
