# Part III, step P3-1b: corrections after external review, ablations and detector robustness

**Status:** specification, 30 October 2026, Toronto time · **Builds on:** [part3_p31_spec.md](part3_p31_spec.md) and the validated code `part3/p31/code/chem_p31.py` (md5 cfc1e3db…) · **Why:** an external review of the assumption list (owner's request, 30 October) and our own flags in [Paper 2](../paper/COMPASS_stage4.md), Sections 7–8

No map is run in this step. P3-1b changes the model in three places and then measures how much the Stage 4 results depend on each built-in assumption. Item 4 of the review plan, alternative ways to start copying, gets its own specification after this step.

## 1. Remove the cap of one copy start per template and step

**Problem.** A template can start at most one copy per time step, with probability 1 − exp(−λ·dt). At the Stage 4 coexistence points, 99.1–100% of template-steps have λ·dt > 3. Every template therefore copies once per step whatever its rate. A trait that changes the copy rate cannot be selected (Paper 2, Section 7). This is a property of the numerical scheme, not of chemistry.

**New rule** (switch `MULTI_START`, default on in P3-1b, off reproduces chem_p31 bit for bit):

- A copy cycle of a template has three parts: waiting for a start (1/λ), copying (t_copy = l/r, as now) and waiting for the duplex to melt (t_melt).
  - t_melt is 1/k_off of the full-length duplex. In cycle-averaged steps, t_melt is the expected waiting time to the next phase that melts it, and at least half a period. In steps with pulses, it uses the pulse hazard of P3-1 §3.3.
- The number of copies a template makes in a step is n ~ Poisson(dt_free / (1/λ + t_copy + t_melt)), where dt_free is the part of the step in which the template is not already busy.
- **Limits.** n is capped by the template's occupancy, n ≤ ceil(dt_free / (t_copy + t_melt)).
  - Across a site, all copies together draw on free slots, activated building blocks and (for primed starts) matching primers.
  - Where demand exceeds supply, every template's n is scaled down in proportion, as the growth step already does for monomers.
- **Copies.** Each copy is drawn independently: start type (unprimed, or primer of length l_p), length (the existing processivity rule) and errors.
  - The copy made last in the step stays paired to its template at the end of the step if the duplex has not melted; the others are released.
  - The template is busy into the next step for the remaining copy time, as now.
- **Records** (no change to dynamics): per run, the share of template-steps that hit the occupancy cap and the share limited by slots, building blocks or primers. This replaces the old saturation counter.

**Validation (criteria fixed now):**

1. **Off switch.** With `MULTI_START` off, the code is bit-identical to chem_p31 (cfc1e3db) on the 16 stored runs of test 4 and on 8 P3-1 runs (AS2U).
2. **Low-rate limit.** Where λ·dt ≤ 0.05 for all templates, mean copies per template-step agree with the old rule within 3% (two-sided z-test, |z| < 3).
3. **Rate sensitivity.** In a fixture with two template sequences in one site, one with twice the extension rate of the other, at λ·dt ≥ 10:
   - the ratio of copies made matches the renewal formula within 5%;
   - with `MULTI_START` off, the ratio is 1.00 ± 0.02. The old rule's blindness to rate is reproduced.
4. **Resource limits.** With slots or building blocks clamped low, total copies never exceed supply, and the cut falls on templates in proportion to their n.
5. **Kill and resume:** bit-identical.
6. **Memory pre-flight** at the worst case (template-dense Stage 4 coexistence points, batch 512, `LOT_EMPTY_CACHE_EVERY=1`): peak ≤ 80% of the 28 GB VRAM cap, no growth over 5 batches. Multi-start can raise the number of copies per step by orders of magnitude, so this check is required before anything else runs.

## 2. Copying accuracy tied to the named pair, and salt in duplex stability

### 2.1 Accuracy

In P3-1, ε₀ and ε₁ are free layer C constants (10⁻³–10⁻¹). In a 2-letter alphabet, an error is a mismatch: A on A or s2U on s2U for AS2U (G on G or C on C, A on A or U on U for the controls).

- **First task: find measured data.** Extract every misincorporation or mismatch-extension measurement in Fang et al. (2025, main text and SI, already on the B360) and in the sources of [part3_p30_sources.md](part3_p30_sources.md). Look for the relative rate of a mismatched bridged dinucleotide or monomer, or for a measured fidelity.
- **If measured values exist** for the pair, ε per letter becomes ε = k_mis/(k_mis + k_match) at the same activated-block concentrations, with the measured ratio. Its uncertainty becomes a narrow layer C range around the measured value.
- **If none exist,** ε stays a layer C constant and the report says so. For each pair, the P3-1b report then also gives the largest ε at which heredity survives (from test 7 of P3-1: ε·L_copy ≈ 0.85), as the number for experiments to beat.
- Each source is read in the original before use, as for all earlier steps.

### 2.2 Salt

P3-1 uses nearest-neighbour parameters for 1 M NaCl with no salt correction. Duplex stability depends strongly on ionic strength, and salinity (0–1.2 M) and Mg²⁺ (0.1–100 mM) are layer A variables.

- **Correction.** Apply a published salt correction for RNA or nucleic-acid duplexes to ΔG°(T). It uses the run's salinity and Mg²⁺ and enters through ΔS, so the correction grows with duplex length. Mg²⁺ enters by a published monovalent-equivalent rule.
- **Which formula** is chosen after reading the originals. Candidates to read: SantaLucia (1998) for the monovalent entropy term; Owczarzy et al. (2008) for Mg²⁺; any RNA-specific correction that NNDB cites. The chosen formula, its stated validity range and the extrapolation outside it are recorded in the implementation table.
- **Validation.** For the Fang et al. 9-pair duplexes, the code reproduces the formula's ΔG° shift at 3 salt conditions within 0.1 kcal/mol. With salt correction off, the code is bit-identical to chem_p31.

## 3. Ablations and detector robustness at the Stage 4 points (`PAIR=None`)

This answers the reviewer's main request: how much do the Stage 4 levels depend on each built-in rule? It runs on the Stage 4 chemistry (`PAIR=None`, the 4c-3/4c-4/4c-5 path inside chem_p31), so the comparison is with the published results.

**Points:**
- the 7 confirmed L2 points (#1558, #3709, A, B, C, rep1_1926, rep1_9823);
- the 10 coexistence point groups of the 4c-5 point study, for E3, joint and L3.

**Runs:** 16 seeds per point and variant, with the paired controls of each level (random-copy ×2; β = 0 where L3 is scored).

| Variant | Change | What it tests |
|---|---|---|
| base | none (reproduces the published rates) | reference; must match the stored results statistically |
| A1 growth cap | untemplated growth cap 4 → 16 letters per chain and step | whether L2 relies on the cap leaving more activated monomer for copying |
| A2 multi-start | §1 on | whether the one-start cap shapes L2, L3 and L4 |
| A3 no extra enclosure | f_rehyd = 1 | whether vesicle content (E3 strict, joint, L3) depends on rehydration enclosure |
| A4 duplex protection | paired strands hydrolyse 10× more slowly | the stated simplification of 4c-3 |
| A5 phase-resolved melting | time step limited to period/4 at cycling points (no cycle averaging of k_off) | the arithmetic-mean melting flag |
| A6 size-dependent washout | chains of length ≥ 4 wash out at 1/τ × 0.5 | the uniform-washout simplification |

A5 is costly at long periods. It runs at the 3 L2 points with the shortest periods first, and at the others only if the pre-flight cost estimate is under 4 hours in total.

**Detector perturbation** (no new runs: all scores are recomputed from the stored per-checkpoint records of the base runs; where a record is missing, the base runs store it):

| Level | Thresholds varied |
|---|---|
| L2 strict | margin 0.10 / 0.15 / 0.20; persistence 4 / 5 / 6 turnovers; L_H 5 / 6 / 7 |
| L3 strict | identity 85 / 90 / 95%; letter fraction 0.4 / 0.5 / 0.6; ancestor 2 / 3 / 4 generations up |
| L3 loose | ρ threshold 0.05 / 0.1 / 0.2 |

For every perturbation, report:
- the share of points that keep the level;
- the false-positive rate in the paired controls at the same threshold.

**Reading, fixed now.** A level counts as **robust** to a variant if its rate at a point changes by less than half its base value, or by less than 3 of 16 runs. Every level and variant whose result changes beyond that is listed in the report as "depends on X". The Paper 2 assumptions table is then updated to match, with the date of the revision.

## 4. Order of work and stop rule

1. §1 implementation and its validation tests 1–6.
2. §2 literature extraction, then implementation and validation.
3. §3 base runs, then the ablations, then detector perturbation.

- If any validation test fails, the run stops and the failure is reported. Nothing is redesigned without review.
- Implementation choices are recorded in a table, as in earlier steps, before any result is read.
- No P3-1 map runs in this step.

## Sources to verify before use

- Fang, Z., Jia, X., Xing, Y. & Szostak, J. W. (2025). *PNAS* 122(21), e2505720122 (full text and SI on the B360).
- SantaLucia, J. (1998). *PNAS* 95(4), 1460–1465. Unified nearest-neighbour parameters, including the salt correction (to be read).
- Owczarzy, R. et al. (2008). *Biochemistry* 47(19), 5336–5353. Mg²⁺ correction for DNA duplexes (to be read; applicability to RNA to be checked).
