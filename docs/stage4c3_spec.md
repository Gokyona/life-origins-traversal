# Stage 4c-3 specification: template copying and heredity (level L2)

**Status:** specification, drafted 26 October 2026 (15:35 Toronto time), revised during validation on 26–27 October (see §6) · **Plan:** [stage4_plan.md](stage4_plan.md), step 4c-3 · **Builds on:** [4c-1](stage4c1_spec.md) (polymer chemistry), [4c-2](stage4c2_spec.md) (catalysis, L1)

The third layer adds template-directed copying to the chemistry of 4c-2. Catalysis stays on. The question is whether, and in which environments, sequence information is copied well enough that sequences outlive the molecules that carry them. This is level L2, "heredity". No rule says "replicate": copying is a chemical reaction on a template, with errors, and it needs the two strands to come apart again.

## 1. Copying rule

- **Who copies.** Only nucleotide-like chains (chemistry 0) act as templates or are made by copying. Amino-acid-like chains do not pair. This follows known chemistry: base pairing is a property of nucleotides.
- **Templates.** A free (single-stranded) nucleotide-like chain of length ≥ L_t (layer C, 4–12) can serve as a template.
- **Copying step.** A copy starts at a random position of the template and is extended with activated monomers of chemistry 0, one position at a time, each the complement of the template letter (alphabet 2: 0 ↔ 1). Each position is wrong with probability ε (layer C, 10⁻³–0.3, logarithmic).
- **Rate.** Copies start at rate k_c · f_min · k_t · Arr_cat(T) · c_eff · (M*₀/V) per free template (k_t: templated speed-up relative to uncatalysed condensation, layer C, 10–10⁵, logarithmic; same temperature dependence as catalysed ligation). The monomers come out of the M*₀ pool, as in growth.
- **Copy length.** Extension continues position by position with probability q = r / (r + k_off(l, T)) per position, where l is the length copied so far, r = k_t-scaled extension rate and k_off the dissociation rate of the l-long partial duplex (below). The copy stops at the template end or when it falls off. The copy is released as a new chain once it has at least 2 letters.
- **Duplexes (product inhibition).** A finished copy of length l stays paired to its template. While paired, the template cannot be copied again and neither strand can be ligated or act as a catalyst. The duplex comes apart at rate k_off(l, T) = k₀ · exp(ΔG(l, T) / RT), with a two-state duplex energy ΔG(l, T) = (l − 1)(ΔH_bp − T ΔS_bp) (ΔG < 0 for a stable duplex). Two layer C constants fix it: ΔH_bp (25–45 kJ/mol per pair, magnitude) and T_m8, the melting temperature of an 8-pair duplex at reference conditions (10–90 °C), from which ΔS_bp follows. k₀ = 10⁶ /h. In cycling environments the hot phase (T + 10 d) melts duplexes that the cool phase made. This is how the map can show whether copying needs temperature cycling, which is a known problem for non-enzymatic copying.
- **Hydrolysis and washout.** Paired strands are hydrolysed like single strands (no protection: a stated simplification) and washed out together.

New layer C constants (5): L_t, k_t, ε, ΔH_bp, T_m8. The map then has 42 dimensions.

## 2. Detector for L2

L2 means that sequence information persists longer than the molecules carrying it, because it is copied.

1. **Lineage record.** At the 16 checkpoints of each sampled site (as for L1), store the abundant nucleotide-like sequences of length ≥ L_h = 6 (≥ 2 copies, top 128), together with their complements.
2. **Persistence of information.** For the window from checkpoint k to the end, compute the share of chain mass at the end whose sequence (or its complement) was already present at checkpoint k. Express the window length in molecule turnover times (mean lifetime of a length-L_h chain against hydrolysis and washout, as for L1).
3. **Score L2.** L2 holds if E1 holds and:
   - **Lasting:** the longest window with information share ≥ 0.2 spans ≥ 5 turnover times;
   - **Copied, not just stable:** the information share is higher than in a paired control run (same condition and seed) where copying makes chains of the same length and composition but with random sequence (copying on, information transfer off). Two versions are recorded: *L2 strict* (lasting, and share above control by ≥ 0.1) and *L2 loose* (lasting only). The map targets L2 strict.
4. L2 is linked to E1 only. Its overlap with E2 and L1 is reported as a result, as for L1.

Composition record fields added by this layer: copies made and their mean length, measured error rate, share of templates paired at the end, longest persistent lineage (in turnover times), the most abundant lineage's sequence length, the order of first L1 and first L2.

## 3. Validation before any map

1. **Copy fidelity.** A fixed template, no hydrolysis or washout: the measured error rate of copies equals ε within sampling error; copies are complements of the template segment; the copy-length distribution matches the processivity formula.
2. **No copying.** With k_t → 0, L2 strict never occurs.
3. **Random-copy control.** With random-sequence copying as the main run, L2 strict never occurs.
4. **Error threshold.** At a reference condition where L2 occurs, raising ε so that ε · L_h goes from 0.1 to 3 lowers the information share and removes L2 (an Eigen-type error threshold: information is lost when more than about one error is made per copy).
5. **Strand separation.** At a temperature where 8-pair duplexes are stable (T well below T_m8) and without cycling, templates end up paired and copying stalls; with temperature cycling that crosses T_m8, copying continues. Far above T_m8, copies fall off before they are long enough, and copying again fails.
6. **Power cut.** Kill and resume gives bit-identical results, as in 4c-1.

Ranges to check against original sources before the map: the duplex parameters (nearest-neighbour thermodynamics of nucleic acids), and published rates and error rates of non-enzymatic template copying. Until checked, the spec uses wide layer C ranges only.

## 4. Map

A joint map of E1, E2, L1 and L2 over the 42 dimensions, COMPASS only, ground truth 1000 × 16, budget 32k, as in 4c-2. The random-copy control runs only where the main run has a lasting lineage. If L2 turns out as rare as L1, a focused map follows, as for 4c-2. One GPU job at a time; every job in its own memory-limited scope.

## 5. What this layer cannot show

Copying ability is set by layer C constants, not derived from chemistry. The map therefore says how fast and how accurate copying would have to be, in each environment, for heredity to arise. Whether real prebiotic nucleotides could copy that well is an experimental question; the map's threshold for k_t and ε is the number to compare with experiments. Paired strands are not protected from hydrolysis, and there is no explicit Mg²⁺ or salt effect on duplex stability; both are limits of this layer.

## 6. Revisions made during validation (26–27 October 2026)

Validation showed that several rules in §1–§2, as first written, suppressed heredity by construction or made the detector unreliable. Each change below was made before any map was run. Results are in [stage4c3_validation_b70.md](../results/stage4c3_validation_b70.md).

| # | Rule as drafted | Problem found | Rule now |
|---|---|---|---|
| 1 | Copies start at the template without a primer | The first pair is unstable (k_off(1) = k₀), so q(1) ≈ 0 and almost no copy starts | Primed initiation: a 2–3-letter chain from the pool binds and is extended; start rate ∝ n_primer · 2^(−l_p); the primer is consumed. Unprimed starts remain possible |
| 2 | Copies start at a random position of the template | A full-length copy needs a start at position 0, probability 1/L; no lineage can re-make itself faster than it decays, so heredity is impossible whatever the copy chemistry | Copies start at the template's first letter and run to its end; a primer must be the complement of the template's first l_p letters (the textbook picture of a primer at one end) |
| 3 | Copying finishes within one time step | Copy time does not depend on length, so the classic advantage of short parasites (a shorter time per copy) is missing | Finite copy time: a template stays busy for l / r hours |
| 4 | Copies and new chain starts drew on free slots independently | Copying could starve chain initiation, or the reverse | Free slots are shared in proportion to the two rates |
| 5 | Control: a separate run with random-sequence copying | Not paired with the main run; biased (9/16 false positives with copying off) | The control is a rerun of the same candidate with identical random streams; only the copied letters differ |
| 6 | Score: share of chain mass at the end whose sequence was present at checkpoint k, best of 8 sampled sites, margin ≥ 0.1 over one control | With copied letters fully random (no information), 15–30% of runs still passed: the best-of-8 statistic is noisy | S is computed over all sites (64 at screening). Strict requires S(main) − max(S(control 1), S(control 2)) ≥ 0.15, with two independent random-copy controls. The margin was calibrated on 2 × 64 null runs (false positives 1/64 and 0/64; smallest margin with ≤ 2/64 is 0.142) |
| 7 | Turnover time from hydrolysis and washout only | Adsorbed chains are protected from hydrolysis | Turnover time includes adsorption protection (changes the 4c-2 L1 result almost nothing) |
| 8 | Results reproducible from seeds | Some GPU reductions are not ordered, so runs with heavy copying differed by 0.1–1% between processes | Deterministic GPU algorithms and a fixed tie-break in the choice of abundant sequences; kill and resume is now bit-identical |

**Lineage identity.** Lineages are still keyed on whole sequences. Uncatalysed growth and ligation append letters at the end of a chain, so a faithfully copied chain that later grows counts as a new sequence. This lowers the L2 rate, and it is a stated limitation. Keying on the first k letters was tested and rejected: with a 2-letter alphabet, k ≤ 8 gives too many chance matches (no margin keeps the null false-positive rate ≤ 2/64 at k = 8), and prefixes of k = 10 are never abundant enough to detect anything.

**Cost.** 0.31 s per run for the main simulation and 0.65 s with both controls at screening fidelity on the B70, above the 0.3 s target; a full map (48k runs) takes about 9 hours.
