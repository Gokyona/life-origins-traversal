# Stage 4c-3: where does heredity arise? (Intel Arc Pro B70)

**Dates:** 26 October 2026, 16:00, to 27 October 2026, 20:50 Toronto time · **Device:** Intel Arc Pro B70 32 GB (PyTorch XPU) · **Specification:** [stage4c3_spec.md](../docs/stage4c3_spec.md), with the revisions in its §6 · **Validation and run log:** [stage4c3_validation_b70.md](stage4c3_validation_b70.md) · **Code:** [stage4/](../stage4/README.md) (chem4c3.py, cat4c3.py, run4c3_map.py, validate4c3.py, scan4c3.py and others, md5-checked) · **Data:** the result JSONs are kept with md5 sums on the benchmark machine and will be added to the repository.

**Alphabet:** 2 letters (0/1), on purpose. The question is whether the simplest alphabet can give rise to heredity, and where. It is not a simulation of modern 4-letter RNA.

## Bottom line

We added template copying to the chemistry of 4c-2. A free nucleotide-like chain can be copied letter by letter into its complement, with errors. The copy stays paired to its template until the duplex melts. Then we looked for **heredity (level L2)**: sequences that outlive the molecules carrying them because they are copied. We compare against a paired control in which copying makes chains of the same length but with random sequence.

1. **Heredity exists in the model, but only in a narrow set of environments.** Five environments give heredity reproducibly (12–16 of 16 independent runs), with copies of about 60 letters.
2. **What the five have in common:**
   - The temperature cycles across the melting point of short duplexes: the cool phase pairs and copies, the hot phase separates copy from template, so both can serve again.
   - Chains live long compared with a cycle. In the three points examined in detail, washout sets the lifetime, not hydrolysis.
   - Copying is fast and accurate enough: fewer than about 0.3 errors per full copy.
3. **Random sampling does not find it.**
   - In the full 42-dimensional space, 16,000 ground-truth runs and 32,000 COMPASS runs found no real heredity. The 77 runs that passed the statistical test copied only about 2 letters per copy and were noise.
   - Even inside the focused region, where temperature, melting point, lifetime and copy speed are in a favourable range, ground truth found none in 8,000 runs. COMPASS found 3 in 16,000.
4. **The error threshold holds.** At two of the points, heredity disappears when the errors per full copy rise from 0.3 to 1–1.5. This is Eigen's error threshold (Eigen 1971), here emerging from the chemistry rather than being put in.
5. **Short parasites take over** when a long and a short template compete, as known from experiments and theory: the short one ends with 99.4% of lineage mass.
6. **One of the first points was an artifact.** It depended on hydrolysis being too slow in strong acid, where the rate law has no acid branch, and was dropped. The other five are robust to this.

## What was simulated

The 4c-2 chemistry (feedstock, activation, condensation, hydrolysis, washout, UV, cycling, catalysis) plus template copying ([specification §1](../docs/stage4c3_spec.md)):

- **Templates and copies.** Only nucleotide-like chains of at least L_t letters act as templates. A copy starts at the template's first letter (primed or unprimed), adds one complementary letter per step, and is wrong with probability ε at each step.
- **Copy length.** Set by processivity: at each step the copy continues with probability r / (r + k_off(l, T)), where k_off is the dissociation rate of the partial duplex.
- **Duplexes.** A finished copy stays paired and blocks its template until the duplex melts at rate k_off(l, T). This uses a two-state duplex energy fixed by ΔH_bp and T_m8, the melting temperature of an 8-pair duplex.
- **New layer C constants:** L_t, k_t (copy speed-up), ε, ΔH_bp, T_m8. The model has 42 dimensions: 13 environmental (layer A) and 29 unknown constants (layer C).

**Detector.**
- L2 strict: a lineage persists for at least 5 molecule lifetimes, and its information share S beats both of two random-copy controls by at least 0.15 (null false-positive rate 1/64).
- L2_long: L2 strict plus a mean copy length of at least 6 letters. It was added after the full map showed that very short copies can pass the strict test at noise level.

## Validation

All tests passed before the results below were accepted, with two recorded exceptions. Full details are in the [validation page](stage4c3_validation_b70.md).

| Test | Result |
|---|---|
| Copy fidelity | measured error 0.0494 against ε = 0.05; length distribution matches processivity (χ² p = 0.98) |
| No copying / random copying | no strict L2 *(earlier code)* |
| Error threshold, positive control | strict 15/16 at ε·L_H = 0.1, 0 at 0.3–3; one secondary criterion missed (0.083 against 0.073), reported |
| Strand separation | copying stalls without cycling below T_m8 and works with cycling across it |
| Kill and resume | bit-identical runs, selections and curves |
| Parasite takeover | short template's share 0.994 ± 0.006 |
| Null false positives | strict 1/64 and 0/64 |
| Error threshold at natural references (test 9) | as pre-registered: pass at #3709, **fail at #1558**, because heredity was already gone at the first grid point |
| Test 9b, finer grid, criterion fixed before running | **pass at both points** (table below) |

## Results

### Full-space map (42 dimensions)

Ground truth 1000 × 16 runs, COMPASS 32,000 runs, 5.6 h.

| Level | Ground-truth fraction | COMPASS accuracy | Most important variables |
|---|---|---|---|
| E1 (activated feedstock) | 0.466 | 0.961 | residence time τ, T, Y_UV (C), UV, thermal cycling, H₂, feedstock |
| E2 (long chains) | 0.085 | 0.959 | T, τ, feedstock, Y_UV (C), UV, pH |
| L1 (catalytic network) | 0.001 | – | too rare to map |
| L2 (heredity) | 0 | – | too rare to map |

### Focused map of the heredity region (box2)

The region: T up to 15 °C below T_m8, a hot phase above T_m8, chain lifetime ≥ cycle period (otherwise infeasible and not simulated), high k_t, strong thermal forcing. Ground truth 500 × 16, COMPASS 16,000.
- Ground truth: E1 0.452, E2 0.184, L1 0.002, L2_long 0. 46% of the region is infeasible.
- COMPASS accuracy for E1 and E2: 0.968 and 0.95. Importance is led by τ, pH and cycle period.
- COMPASS found 3 L2_long runs (A, B, C below). Each was confirmed with 16 new runs.

### The five heredity environments

| Point | T / T_m8 (°C) | Cycle period | pH | ε (per letter) | Strict / 16 | Margin | R |
|---|---|---|---|---|---|---|---|
| #1558 | 7.6 / 10.2 | 118 h | 7.96 | 0.0036 | 16 | 0.42 | 0.98 |
| #3709 | 38.7 / 46.3 | 33 h | 4.57 | 0.0063 | 16 (15 in test 9) | 0.61 | 1.05 |
| A | 14.2 / 26.2 | 1.1 h | 3.83 | 0.0017 | 16 | 0.54 | 0.93 |
| B | 21.4 / 25.2 | 3.8 h | 2.64 | 0.0074 | 12 | 0.33 | 0.98 |
| C | 26.7 / 37.2 | 26 days | 4.73 | 0.0013 | 16 | 0.61 | 0.96 |

- R is the number of long copies per template per molecule lifetime. It separates the classes: 0.4–1.05 where heredity was found, ≤ 0.05 elsewhere.
- Copies are about 62 letters long in all five.
- **Acid check.** The hydrolysis law has no acid branch, so lifetimes at pH < 5 might be too long. We reran the acidic points with hydrolysis computed at pH 5:
  - #3709, A, B and C are unchanged (#3709 16/16, margin 0.619).
  - At A, B and C, hydrolysis of a 64-mer is only 2×10⁻⁴ – 8×10⁻⁴ of the washout rate even at pH 5, so it does not matter.
  - The dropped point #813 (54 °C, pH 2.7) went from 16/16 to 0/16, because at that temperature hydrolysis is about 30% of all loss.

### Copy-length cap

The model can hold at most 64 letters per chain, a hard limit of the sequence encoding, and at the heredity points nearly all chains sit at it. With the cap at 48 and at 32, heredity at #1558 and #3709 stays (15–16 of 16), copy length follows the cap, and the margin rises (#3709: 0.61, 0.74, 0.86). **The cap sets how long the inherited sequence is, not whether heredity happens.**

### Error threshold (test 9b)

x = errors per full copy = ε × 62.

| x | #1558 strict / 16 | #1558 S | #3709 strict / 16 | #3709 S |
|---|---|---|---|---|
| 0.15 | 16 | 0.542 | 16 | 0.873 |
| 0.3 | 16 | 0.381 | 16 | 0.697 |
| 0.5 | 5 | 0.243 | 16 | 0.497 |
| 0.7 | 1 | 0.156 | 16 | 0.329 |
| 1.0 | 0 | 0.096 | 1 | 0.157 |
| 1.5 | 0 | 0.030 | 0 | 0.050 |
| 3 | 0 | 0.003 | 0 | 0.005 |

The criteria were written into the output file before the first run: S non-increasing, strict ≥ 12/16 at x ≤ 0.3, strict ≤ 2/16 at x ≥ 1.5. All hold at both points.

## Which conclusions follow directly from the assumptions

| Conclusion | How it arises | Verdict |
|---|---|---|
| Heredity needs temperature cycling across T_m8 | Duplexes block templates by rule; only melting frees them (test 5) | **partly built in**: the need for strand separation is put in. Which cycle periods and depths work is derived. |
| Chains must live longer than a cycle | One full copy per template per cycle at most, so lifetime / period bounds R | **partly built in**: follows from the copying rule. The focused region was chosen with this bound. |
| Error threshold near 0.5–1 errors per copy | Measured on emergent lineages against random-copy controls | **derived**: the threshold is not coded, ε is only a per-letter error rate |
| Short parasites win | Shorter templates finish copies sooner | **derived** from the copy kinetics; agrees with known results |
| Heredity is very rare in environment space | Counted over ground truth and COMPASS | **derived**, within the layer C ranges. Wider ranges for k_t or ε would change the frequency. |
| Five quite different environments work (1 h to 26 d cycles, 8–39 °C, pH 2.6–8) | Confirmed with 16 runs each | **derived**; the common factors are the two above |
| Acidic points are real | Hydrolysis is negligible next to washout there | **derived** for A, B, C and #3709; depends on the washout and hydrolysis laws |
| Inherited sequences are about 60 letters | They fill the 64-letter cap | **built in**: the cap sets the length (cap test) |

## Limits

- **Copying ability is set by layer C constants** (k_t, ε, L_t, ΔH_bp, T_m8), not derived from chemistry. The result is how fast and accurate copying would have to be. The five points need ε of 0.0013–0.0074 per letter and k_t of 10⁴ or more. Whether prebiotic nucleotides could copy that well is an experimental question.
- **Saturated regime.** At the heredity points nearly all copy slots hold chains at the length cap, and R ≈ 1 may be limited by the slot count. Heredity here is the maintenance of a full population, not growth into empty space.
- **No acid branch in the hydrolysis law**, and pH 2.6 (B) is an extrapolation. The acid check above bounds the effect for the five points.
- **Paired strands are not protected from hydrolysis**, and there is no explicit Mg²⁺ or salt effect on duplex stability.
- **Lineages are keyed on whole sequences**, which undercounts heredity when a copied chain later grows.
- **Two-letter alphabet**, so random sequence matches are frequent; the detector's controls account for this.
- **Few points.** The L2 maps do not resolve a boundary. Statements about where heredity occurs come from the five confirmed points and the targeted region, not from a classifier. A second COMPASS repetition of the focused map is running to find more points and test importance stability.

## References

- Eigen, M. (1971). Selforganization of matter and the evolution of biological macromolecules. *Die Naturwissenschaften* 58(10), 465–522.
- Li, Y. & Breaker, R. R. (1999). Kinetics of RNA degradation by specific base catalysis of transesterification involving the 2′-hydroxyl group. *Journal of the American Chemical Society* 121(23), 5364–5372.
- Further references and checks: [stage4c3_spec.md](../docs/stage4c3_spec.md), [stage4c2_spec.md](../docs/stage4c2_spec.md), [stage4c1_spec.md](../docs/stage4c1_spec.md).
