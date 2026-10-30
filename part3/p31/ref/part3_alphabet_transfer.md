# Part III: can measured 4-letter RNA rates be used for the 0/1 alphabet?

**Status:** first analysis, 30 October 2026, Toronto time · **Design:** [part3_design.md](part3_design.md) · **Script:** [part3/alphabet_tm.py](../part3/alphabet_tm.py) (CPU, seconds)

The owner asked for this to be settled before any further work, because every later inference depends on it.

## 1. The short answer

**Rates averaged over four letters cannot be used for a 2-letter system.** Real 2-letter chemistry is always one specific base pair (or a modified pair), and its properties differ by orders of magnitude depending on which pair it is. What *can* be used is the part of the 4-letter data that belongs to that pair. Nearest-neighbour thermodynamics, dinucleotide extension rates and per-base fidelities are all measured pair by pair, so a G/C or an A/U alphabet has its own measured values, taken as a subset of the 4-letter tables.

Doing this exposes a problem in the Stage 4 results: **the heredity environments found in Stage 4 need pairing as weak as A·U but copying as accurate as G·C.** No single real base pair has both.

## 2. Which quantities transfer, and how

| Quantity | Transfers? | How |
|---|---|---|
| Duplex stability (ΔH, ΔS, melting) | **Yes, exactly** | The nearest-neighbour model is a sum over stacks. A G/C alphabet uses only the GG/CC, GC/CG and CG/GC stacks; an A/U alphabet uses only AA/UU, AU/UA and UA/AU (Turner 2004 parameters, Xia et al. 1998) |
| Duplex association rate | Yes, roughly | About 10⁶–10⁷ M⁻¹ s⁻¹ for short strands at high salt, weakly sequence-dependent (Ouldridge et al. 2013) |
| Extension rate per letter | **Yes, as a subset** | Ding et al. (2022) measured all 16 bridged dinucleotides. G/C: C\*G 84 h⁻¹ down to slower G/C combinations; A/U: down to U\*U 0.35 h⁻¹ (at 5 mM) |
| Copying fidelity | **No, not from 4-letter averages** | In a 2-letter system the only possible error is the other letter, while 4-letter error rates are dominated by G·U wobble and A:G mismatches. Fidelity must come from 2-letter experiments or from the matching mismatch types |
| Hydrolysis | Yes | The backbone is the same; base effects are smaller than pH and temperature effects (Li & Breaker 1999) |
| Catalysis (ribozymes) | Partly | A 2-letter ligase ribozyme exists (2,6-diaminopurine and uracil; about 36,000-fold rate enhancement; Reader & Joyce 2002), so 2-letter catalysis is possible, but with weaker structures |

## 3. Result: what the Stage 4 constants correspond to

Stage 4 described pairing by two constants: ΔH_bp (enthalpy per pair, 25–45 kJ/mol) and T_m8. T_m8 is the temperature at which an 8-pair duplex comes apart at rate k₀ = 10⁶ per hour. For real RNA, the same quantity is the temperature at which k_on · exp(ΔG/RT) = k₀. We computed it for all 256 possible 8-mers of each alphabet, using the nearest-neighbour parameters from the NNDB and k_on = 10⁶–10⁷ M⁻¹ s⁻¹ (5th–95th percentile across sequences):

| Alphabet | ΔH per pair (kJ/mol) | Equivalent T_m8 (°C) | Tm of 8-mer at 1 µM (°C) |
|---|---|---|---|
| G/C | 53–56 | **78–92** | 65–72 |
| A/U | 30–35 | **2–25** | −15 to 1 |
| 4-letter (reference) | 38–52 | 34–80 | 19–59 |
| **Stage 4 confirmed heredity points** | layer C range 25–45 | **11.9–30.1** | — |

The Stage 4 points match the **A/U** alphabet. G/C pairing lies outside the Stage 4 ΔH range altogether, and its T_m8 is far above every confirmed point.

Now set this against what copying needs:

| Requirement at the Stage 4 points | A/U alphabet | G/C alphabet |
|---|---|---|
| Error per letter ε ≤ 0.0013–0.0074 | **fails**: in 4-letter copying, A and U are copied with only about 54–88% fidelity at the first position, and A/U-only data are missing | **plausible**: a G/C-only system reached 99.7–99.9% (ε about 0.001–0.003) at the first position (as reported in the NERPE-Seq preprint; primary source still to be checked) |
| Fast copying | **fails**: U\*U is the slowest bridged dinucleotide (0.35 h⁻¹ at 5 mM); without helper strands A and U are at least 1,000 times slower than G and C (Prywes et al. 2016) | **best case**: C\*G 84 h⁻¹ |
| Duplexes that come apart in cool cycles (T_m8 12–30 °C) | **fits** | **fails**: strands separate only near 80–90 °C, where RNA backbone cleavage is about 10⁴ times faster than at 25 °C (0.07 log units per °C; Li & Breaker 1999) |

So the Stage 4 heredity result rests on a combination of constants that no measured base pair provides. This is a second, independent reason, beside the copying speed, to treat the Stage 4 heredity levels as conditional.

## 4. What this means for Part III

1. **The 0/1 alphabet must be a named chemistry.** Each run states which pair the letters stand for, and every rate comes from that pair's measurements. "Borrowed 4-letter rates" is dropped as a formulation.
2. **Three candidates, to be tested side by side:**
   - **G/C:** fast and accurate, but separation needs hot phases; this moves heredity towards hot, cycling settings (hot-spring fields), where hydrolysis is the limit.
   - **A/U:** separates in cool cycles, but slow and inaccurate.
   - **An intermediate pair**, such as 2,6-diaminopurine·U (the pair of the only known 2-letter ribozyme, Reader & Joyce 2002) or 2-thio-U·A, which pair more strongly than A·U. Whether their stability, speed and fidelity fall in the window Stage 4 needs is the key open question.
3. **The question becomes testable:** is there *any* 2-letter chemistry with real measured properties for which heredity can start in some early-Earth environment? If only the intermediate pairs work, that is a concrete prediction for experiments.
   *First numbers for D·U (added the same day):* replacing A by 2,6-diaminopurine stabilises a duplex by 0.4–1.2 kcal/mol per substitution at 37 °C (2′-O-methyl RNA/RNA, 0.1 M NaCl; Pasternak et al. 2007, *Nucleic Acids Research* 35, 4055–4063, doi:10.1093/nar/gkm421). Applying this crudely to every pair of an A/U 8-mer gives an equivalent T_m8 of about 28–44 °C at 0.4 kcal/mol, 55–70 °C at 0.9 and 72–86 °C at 1.2. Only the weakest stabilisation reaches the Stage 4 window (12–30 °C), so D·U is not an obvious fit either. Its copying speed and fidelity are not yet known to us. A 2-thiouridine study (Ding et al. 2024, *JACS*, doi:10.1021/jacs.3c11158) reports that s2U pairs with A and also with itself, the latter about as strongly as a native A·U pair, with little cost in copying fidelity; numbers not yet read. Fang, Jia, Xing & Szostak (2025), *PNAS* 122, e2505720122, doi:10.1073/pnas.2505720122 (abstract read; full text not yet): in an alphabet of s2U, s2C, inosine and A, the I:s2C pair is isomorphic and isoenergetic with A:s2U; A:s2U is stronger than A:U and I:s2C weaker than G:C, so all pairs have similar energies, between A·U and G·C. Either of these pairs is therefore a candidate "intermediate" 2-letter chemistry, with 2-thiopyrimidines proposed as plausible precursors of the canonical pyrimidines. Their numbers (ΔG, rates, fidelity) are in the full text, which is the next source to obtain.


## 5. What the Stage 4 maps already say about G/C-like pairing

A CPU-only re-analysis of the existing Stage 4 results on the B360 (30 October 2026, 13:05–13:3x Toronto time; no new simulations; summary, tables, script and MD5 sums in `lot-p3/alphabet_check/` on the B360, to be uploaded with the Part III data):

- In about 30,700 simulated runs with T_m8 between 70 and 90 °C, **no run reached lasting heredity (L2_long)**. At the hit rate seen at T_m8 10–30 °C, about 9 would have been expected. The region was searched, so the absence is a result, not a gap in sampling (Poisson chance of 0 when 9 are expected: about 10⁻⁴).
- The cause is hydrolysis: duplexes that melt only near 70–90 °C need hot phases at which chains break down faster than they are copied.
- ΔH_bp of 53–56 kJ/mol, the G/C value, was never sampled (the range stopped at 45), so this covers a G/C-like melting temperature, not the full G/C chemistry.

**Reading.** Within the Stage 4 model, a G/C-like alphabet does not give heredity, and an A/U-like one gives it only with copying accuracy and speed that A/U does not have. The intermediate pairs (A:s2U, I:s2C, D·U) are the remaining candidates for a 2-letter chemistry, and they become the main line of the alphabet question.

## 5b. Measured values for the intermediate pairs (Fang et al. 2025)

Full text and supporting information read on 30 October 2026 (PDFs from the University of Chicago repository; MD5 sums match the repository's; copies on the B360 in `lot-p3/literature/`).

**Duplex stability** (Table 1): a 9-pair RNA duplex 5′-CUGA X GUAG-3′ with one varied central pair; 1 M NaCl, pH 8, 5 µM total strands; ΔG at 25 °C.

| Central pair | Tm (°C) | ΔG°25 (kcal/mol) | Difference from A:U |
|---|---|---|---|
| G:C | 55.6 | −15.2 | −1.6 |
| A:s2U | 53.1 | −14.7 | −1.1 |
| I:s2C | 52.0 | −14.6 | −1.0 |
| G:s2C | 50.1 | −13.7 | −0.1 |
| A:U | 46.8 | −13.6 | 0 |
| I:C | 46.2 | −13.0 | +0.6 |

**Copying** (Fig. 3, values read from the figure): bridged homo-dinucleotides filling a 2-letter gap; room temperature, 100 mM MgCl₂, pH 8.

| Substrate on template | K_m (mM) | k_obs,max (h⁻¹) |
|---|---|---|
| U\*U on AA | 34 | 0.92 |
| s2U\*s2U on AA | 28 | 26 |
| A\*A on UU | 0.38 | 10 |
| A\*A on s2U s2U | 0.071 | 16 |
| C\*C on GG | 0.39 | 39 |
| G\*G on CC | 0.0060 | 24 |
| s2C\*s2C on II | about 6.9 | 74 |
| I\*I on s2C s2C | 0.14 | 24 |

No fidelity data, no data on 2- or 3-letter alphabets, no deamination rates.

**What this means.**

1. **Replacing U by 2-thio-U turns the slowest pair into a fast one.** s2U\*s2U is copied about 28 times faster than U\*U, and A\*A binds 5 times more tightly on an s2U template. A 2-letter A/s2U alphabet would copy at 16–26 h⁻¹, close to the G/C pair.
2. **The intermediate pairs are much closer to G:C than to A:U.** Per substituted pair, A:s2U gains about 1.1 kcal/mol over A:U, about 70% of the G:C gain. Applied crudely to all eight pairs of an 8-mer (as for D·U above), this puts an A/s2U or I/s2C alphabet at an equivalent T_m8 of roughly 65–90 °C. That is the G/C-like regime in which the Stage 4 maps found no heredity.
3. **So the fast pairs need hot strand separation, and in Stage 4 hot phases destroyed the chains.** The Stage 4 hydrolysis law is the measured one (Li & Breaker 1999), so the failure is not an artefact of the rate law. It may be an artefact of the cycle shape. Stage 4 used square-wave cycles with periods of 1–1,000 h and equal hot and cool halves, so chains spent at least half an hour, and half of all time, at the hot temperature. A G/C-like 8-mer duplex at T_m8 comes apart in milliseconds (k₀ = 10⁶ h⁻¹). At 85 °C and pH 7, the measured law gives roughly 0.02 breaks per bond per hour. Short, rare hot pulses (seconds to minutes, a few % of the time), as in thermal convection in rock pores, would separate strands at almost no hydrolysis cost. **Stage 4 did not sample this kind of cycle.**

## 5c. Conclusion so far

- A 2-letter chemistry that is fast and accurate enough exists on paper: **A/s2U** (or I/s2C, or G/C). Its duplexes need hot strand separation.
- Whether heredity can start with it depends on an environment Stage 4 did not test: **brief, rare hot pulses** in an otherwise cool setting.
- That is the first concrete simulation of Part III: the Stage 4c-3 copying chemistry with A/s2U-like pairing, the measured rates above, and asymmetric hot pulses (duty fraction and pulse length as new environmental variables), with G/C and A/U as controls.

## 6. Next checks

- Stacking (nearest-neighbour) data for all-A/s2U duplexes; the 65–90 °C estimate assumes single-pair gains add up.
- Fidelity of A/s2U copying (the s2U:s2U self-pair is outcompeted by A:s2U; Ding et al. 2024).
- Pulse-heating environments on the early Earth (thermal gradients in porous rock, geyser and hot-spring margins): pulse length and frequency, from the literature.
- Primary source for the G/C-only copying fidelity, and any A/U-only fidelity data.
- Salt and Mg²⁺ correction of the duplex values (the NNDB parameters are for standard conditions; the salt condition is not stated on the pages read).

## Sources

- Xia, T. et al. (1998). *Biochemistry* 37, 14719–14735. Parameters as tabulated in the NNDB, Turner 2004 set (rna.urmc.rochester.edu/NNDB), read 30 October 2026; the original paper was not read.
- Ding, D., Zhou, L., Giurgiu, C. & Szostak, J. W. (2022). *Nucleic Acids Research* 50(1), 35–45. doi:10.1093/nar/gkab1202.
- Prywes, N., Blain, J. C., Del Frate, F. & Szostak, J. W. (2016). *eLife* 5, e17756. doi:10.7554/eLife.17756.
- Leu, K. et al. (2013). Cascade of reduced speed and accuracy after errors in enzyme-free copying of nucleic acid sequences. *JACS* 135, 354–366. doi:10.1021/ja3095558. RNA: correct C:G extension 4 h⁻¹ against 0.1–0.33 h⁻¹ for mismatches; a first mismatch slows extension 10–100-fold.
- NERPE-Seq preprint, bioRxiv 10.1101/2021.01.02.425068 (authors not on the pages read): 11% mismatches overall; first-position fidelity rG 95.1%, rC 94.2%, rU 88.4%, rA 54%; the G/C-only figures (99.9%, 99.7%) are quoted there from earlier work.
- Ouldridge, T. E. et al. (2013). *Nucleic Acids Research* 41(19), 8886–8895. doi:10.1093/nar/gkt687.
- Li, Y. & Breaker, R. R. (1999). *JACS* 121(23).
- Reader, J. S. & Joyce, G. F. (2002). A ribozyme composed of only two different nucleotides. *Nature* 420, 841–844. doi:10.1038/nature01185 (abstract).
