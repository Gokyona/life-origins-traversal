# Part III, step P3-1: heredity at measured rates with a named 2-letter chemistry and hot pulses

**Status:** specification, 30 October 2026, Toronto time · **Design:** [part3_design.md](part3_design.md) · **Why:** [part3_alphabet_transfer.md](part3_alphabet_transfer.md) · **Rates:** [part3_p30_sources.md](part3_p30_sources.md)

## 1. Question

Can heredity (L2) start from scratch when every copying-related rate is a measured value for a named 2-letter chemistry, and the environment can deliver short, rare hot pulses?

The main chemistry is **A/s2U** (adenosine and 2-thiouridine; letters 0 and 1). **G/C** and **A/U** are run as controls at the two ends of pairing strength.

## 2. Base code

`stage4/chem4c3.py` (md5 cab66933…), the Stage 4c-3 chemistry with copying, duplexes, hydrolysis (Li & Breaker 1999 law), washout, primed initiation and the L2 detectors with paired random-copy controls. New file `part3/chem_p31.py`; Stage 4 files are not changed. Everything not listed below stays as in 4c-3, including the 4c-1 and 4c-2 layers.

## 3. Changes

### 3.1 Named pair chemistry (replaces layer C constants ΔH_bp, T_m8, k_t, and K0)

A run sets `pair ∈ {AS2U, GC, AU}`. The pair fixes:

- **Duplex stability, nearest-neighbour.** ΔH and ΔS of a duplex are computed from the actual sequence, with Turner 2004 stacks for GC and AU (NNDB table values, as in `part3/alphabet_tm.py`), the intermolecular initiation term (+3.61 kcal/mol ΔH, +4.09 ΔG°37) and the terminal AU penalty for A/U ends. For AS2U: AU stacks plus a per-pair offset δ (layer C, 0.6–1.6 kcal/mol stabilisation per pair at 25 °C, centred on 1.1 from Fang et al. 2025, Table 1), applied to ΔH with ΔS unchanged. The duplex dissociation rate is k_off = k_on · exp(ΔG°(T)/RT) with ΔG° for a 1 M standard state.
- **Association.** k_on layer C, 10⁶–10⁷ M⁻¹ s⁻¹ (log), Ouldridge et al. 2013. It replaces K0. It also sets the primer-binding rate.
- **Extension rate per letter.** r = k_max · c / (K_m + c) · f_Mg · f_T, where c is the concentration of activated building blocks of the letter required (mM). k_max and K_m are the measured values for that letter on that template letter (bridged homo-dinucleotides, room temperature, 100 mM Mg²⁺, pH 8):

  | pair | substrate on template | K_m (mM) | k_max (h⁻¹) | source |
  |---|---|---|---|---|
  | AS2U | s2U\*s2U on A | 28 | 26 | Fang et al. 2025, Fig. 3 |
  | AS2U | A\*A on s2U | 0.071 | 16 | same |
  | GC | C\*C on G | 0.39 | 39 | same |
  | GC | G\*G on C | 0.0060 | 24 | same |
  | AU | U\*U on A | 34 | 0.92 | same |
  | AU | A\*A on U | 0.38 | 10 | same |

  f_Mg = ([Mg²⁺]/100 mM)^m_Mg, with m_Mg layer C 0–1.5 (no measured exponent found yet; stated assumption). f_T = exp(−E_ext/R (1/T − 1/296 K)), E_ext layer C 20–80 kJ/mol (not yet measured for these substrates; stated assumption).
- **Errors.** ε per letter, layer C 10⁻³–10⁻¹ (log), separately for each letter. No 2-letter fidelity measurement exists for AS2U; the G/C-only first-position fidelity (99.7–99.9%) is quoted second-hand. The map gives the ε threshold, as in Stage 4.

### 3.2 Activated building blocks (replaces the constant activated-monomer supply of the copying step)

Activated building blocks are made and destroyed explicitly:

- **Activation:** monomers become activated at rate k_act · a, where a (mM) is the activating-agent supply (new layer A variable "activation chemistry", 0–500 mM, log, with 0 allowed) and k_act layer C. Calibration point: about 2 mM bridged A\*A in 2 h from 20 mM AMP with 200 mM isocyanide and aldehyde at room temperature (Aitken et al. 2023). The range of k_act is set so that this point lies inside it, with one decade either side.
- **Loss:** activated building blocks hydrolyse with half-life 2–4 h at room temperature (Ding et al. 2022: about 2 h for A\*U and A\*C; 245 min for C\*C in another buffer), Arrhenius with E_act,h layer C 50–100 kJ/mol; washed out like other species.
- **Cold:** no special rule. The cold benefit seen by Aitken et al. (60% of primers +7 letters at −3 °C against 7% at room temperature) must come out of the temperature terms, and this is checked in validation test 6.

### 3.3 Hot pulses (new environmental variables)

The Stage 4 square-wave temperature cycle stays available. A second cycle shape is added:

- base temperature T (layer A as before),
- pulse temperature T_p (T to 100 °C),
- pulse length τ_p (1 s to 1 h, log),
- pulse fraction f_p (10⁻³ to 0.5, log), so the pulse period is τ_p / f_p.

Cases with f_p = 0.5 and τ_p equal to half the Stage 4 period reproduce the Stage 4 square wave.

**Integration.** Pulses can be much shorter than a model step. Within a step, the effect of n pulses is applied in closed form: each duplex melts with probability 1 − exp(−n · k_off(T_p) · τ_p); each bond breaks with probability 1 − exp(−(n · τ_p · k_h(T_p) + (Δt − n τ_p) · k_h(T))), and the same for activated-block hydrolysis. Copying and activation run at the base temperature only (stated simplification: extension during a pulse is neglected).

Physical settings for pulses (thermal convection in porous rock, geyser and hot-spring margins) are to be sourced in P3-0; until then the pulse ranges are argued, not measured, and are labelled as such.

### 3.4 Time scale and rare events

At measured rates, heredity may need far longer than one run. Each run reports, besides the levels, a **time-to-L2 estimate**: the copy reproduction number R (long copies per template lifetime) and the lineage growth rate. Runs with R > 1 but no L2 in the run time are flagged "growing". The rare-event estimator (P3-3) takes over where this is not enough.

## 4. Layer C and layer A summary

New or changed dimensions: pair (categorical: one map per pair), δ (AS2U only), k_on, m_Mg, E_ext, ε_0, ε_1, k_act, E_act,h, activation supply a (layer A), T_p, τ_p, f_p (layer A). Removed: ΔH_bp, T_m8, k_t, ε (single), K0.

## 5. Validation tests (before any map)

1. **Duplex thermodynamics:** the code's ΔG°(T) for the Fang et al. 9-pair duplexes (A:U and G:C centre, all-canonical flanks) matches the NNDB calculation within 0.3 kcal/mol; for AS2U the offset reproduces the A:s2U − A:U difference of Table 1 (−1.1 kcal/mol at 25 °C) by construction.
2. **Extension kinetics:** a fixed primer/template with clamped activated-block concentration gives the per-letter rate of the table (Michaelis–Menten curve over c = 0.1–10 K_m) within sampling error.
3. **Pulse integration:** with one simulated duplex population and one bond population, the closed-form pulse step equals an explicit fine-step simulation (pulse resolved with 20 sub-steps) within 2% for τ_p = 1 s, 1 min and 1 h.
4. **Stage 4 regression:** with pair chemistry switched off (old constants, square wave), the code reproduces chem4c3 bit-identically on 16 stored runs.
5. **Activation calibration:** at the Aitken et al. conditions, the activated-block concentration after 2 h lies within the calibrated k_act range.
6. **Cold benefit (check, not a fit):** at −3 °C against 22 °C, copying yield with in-place activation rises; report the factor and compare with Aitken et al. (no pass/fail threshold, because two terms carry stated assumptions).
7. **Error threshold:** heredity is lost near ε · L_copy ≈ 1, as in Stage 4 test 9b.
8. **Null false positives:** with copy_random, L2_strict ≤ 2 of 64.
9. **Kill and resume:** bit-identical, as in Stage 4.

## 6. Pre-flight and runs

- Memory pre-flight at the worst case (dense templates, longest runs, largest batch) before any long run; peak ≤ 75% of the 28 GB VRAM cap and of the 24 GB host scope; no growth over 5 batches.
- One GPU job at a time, thermal guard and checkpoints as in Stage 4.
- **Stop after validation.** Maps (one per pair, COMPASS, focused on pulse settings) start only after the owner sees the validation results.

## 7. Which conclusions will follow from assumptions

- Extension during pulses is neglected: this can only make copying slower, so it does not create heredity.
- The AS2U offset is applied uniformly per pair; stacking data for all-A/s2U duplexes are not available.
- m_Mg and E_ext are not measured; the maps report how results depend on them.
- Fidelity is a layer C constant; the result is a threshold to compare with future experiments.
- Pulse ranges are argued until physical sources are added.

## Sources

As in [part3_p30_sources.md](part3_p30_sources.md) and [part3_alphabet_transfer.md](part3_alphabet_transfer.md): Fang, Jia, Xing & Szostak (2025), *PNAS* 122, e2505720122; Ding, Zhou, Giurgiu & Szostak (2022), *Nucleic Acids Research* 50(1), 35–45; Aitken, Wright, Radakovic & Szostak (2023), *JACS* 145, 16142–16149; Li & Breaker (1999), *JACS* 121(23); Ouldridge et al. (2013), *Nucleic Acids Research* 41(19), 8886–8895; Xia et al. (1998), *Biochemistry* 37, 14719–14735 (via NNDB).
