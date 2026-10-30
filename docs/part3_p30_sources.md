# Part III, step P3-0: measured rates and early-Earth ranges (first batch)

**Status:** batches 1–2, 30 October 2026, Toronto time · **Design:** [part3_design.md](part3_design.md)

Every value below was read in the original or in an author or repository copy on the date given; where only part of a source could be read, this is said. Values not yet found are listed under "Still missing". PubMed Central pages could not be read from the working environment (browser check), so other hosts were used.

## 1. Copying chemistry

| Quantity | Value | Conditions | Source |
|---|---|---|---|
| Primer extension by imidazolium-bridged dinucleotides, k_obs | 0.35 h⁻¹ (U\*U) to 84 h⁻¹ (C\*G) at 5 mM intermediate; maximum at saturation 5.7–85 h⁻¹ | 100 mM MgCl₂, 200 mM Tris pH 8.0, 12-nt primer, 26-nt template, 12-nt downstream helper; temperature not found in the text read | Ding, Zhou, Giurgiu & Szostak (2022), *Nucleic Acids Research* 50(1), 35–45, doi:10.1093/nar/gkab1202 |
| Apparent binding constant K_m of the bridged dinucleotide | 0.0071 mM (G\*G) to 54 mM (U\*U) | as above; K_m taken as an approximation of K_d | same |
| Effect of base pairs on rate | 15–240-fold range across the 16 dinucleotides; G·C-rich pairs fastest | as above | same |
| Lifetime of the activated intermediate | half-life about 120–135 min (A\*U, A\*C); C\*C 0.17 h⁻¹ (half-life about 245 min) in another buffer | as above / 90 mM Tris, 90 mM MgCl₂ | same |
| Threo-nucleotide addition | 10-fold slower than a ribonucleotide | with an activated downstream oligonucleotide | Nucleic Acids Research 49(2), 646–656 (2021; authors not legible in the record read) |
| Spontaneous RNA backbone cleavage | about 1.3 × 10⁻⁹ min⁻¹ | 23 °C, pH 6.0, 3.16 M K⁺ | Li & Breaker (1999), *JACS* 121(23) |
| pH dependence (pH 9–13) | log k = 0.983·pH − 14.8 | 23 °C, 3.16 M K⁺ | same |
| Temperature dependence | +0.07 log units per °C; Ea 29 kcal mol⁻¹ | — | same |
| Mg²⁺ dependence | log k rises with slope 0.80 in log[Mg²⁺] (0.005–0.05 M); about 3-fold at 5 mM, pH 9.5 | — | same |

**What this changes.** The Stage 4 check used about one letter per hour (Prywes et al. 2016). With pre-formed bridged dinucleotides, measured extension reaches about 1 letter per minute for G·C-rich pairs. That is up to about 100 times faster than the rate used in the Stage 4 check, but still far slower than the Stage 4 model. The size of the remaining gap will be recomputed in P3-1 with the same units as the Stage 4 paper. The measured rates need conditions that are hard to justify on the early Earth: 100 mM Mg²⁺ (sea water about 53 mM), millimolar activated intermediates that hydrolyse within hours, and a downstream helper strand. **So the bottleneck for P3-1 shifts from the copying step to the supply of activated building blocks.** The model must make and destroy these intermediates explicitly, with rates tied to the environment.

**For the 0/1 alphabet:** we map the two letters to the G·C pair (fastest, tightest binding) as the optimistic case and to the A·U pair as the pessimistic case, and report both.

### 1b. Making activated building blocks in place, ligation and duplexes (batch 2)

| Quantity | Value | Conditions | Source |
|---|---|---|---|
| In-place activation of AMP | bridged A\*A up to 2.0 mM after 2 h (2.7–3.3 mM with other catalysts); at −3 °C A\*A stays above 3 mM for 10 h | 20 mM AMP, 30 mM MgCl₂, pH 8, 200 mM methyl isocyanide and 2-methylbutyraldehyde, 50 mM catalyst, 10 mM 2-aminoimidazole, room temperature | Aitken, Wright, Radakovic & Szostak (2023), *JACS* 145, 16142–16149, doi:10.1021/jacs.3c04635 |
| Copying with in-place activation | ≥1 letter added to 58–68% of primers in 24 h; 87% at −3 °C; at pH 7.0 only 13% | as above, 10 mM each of CMP and GMP | same |
| Longer copying | up to 60% of primers extended by ≥7 letters after 5 days at −3 °C with periodic reagent addition; 7% at room temperature | as above | same |
| Activation of an oligonucleotide | 25% after 6 h (plateau); >90% at 250 mM 2-aminoimidazole and 250–500 mM isocyanide | 30 mM Mg²⁺, pH 8 | bioRxiv 10.1101/2025.07.11.664453 (preprint; authors not on the pages read) |
| Templated ligation without ribozyme | about 1.5–2% in 3 h, about 26% in 24 h | as above | same |
| Ligation vs. monomer extension | extension "almost 100 times more rapidly" than ligation | 100–200 mM Mg²⁺, pH 8–9 | Prywes, Blain, Del Frate & Szostak (2016), *eLife* 5, e17756, doi:10.7554/eLife.17756 |
| Copying A and U | without helper oligomers at least 1,000 times slower than G and C; about 10 times slower with helpers; average fidelity 98% | as above | same |
| Duplex association | about 10⁶–10⁷ M⁻¹ s⁻¹ for short DNA and RNA near room temperature at high salt; falls with temperature, while dissociation rises exponentially | [Na⁺] ≥1 M or [Mg²⁺] ≥0.01 M | Ouldridge et al. (2013), *Nucleic Acids Research* 41(19), 8886–8895, doi:10.1093/nar/gkt687 (review of measurements, introduction) |

**What this adds.** Activation chemistry that runs in place exists, but it needs very high concentrations of reagents (hundreds of mM of isocyanide), and it works much better in the cold: at −3 °C, more than half the primers gained 7 or more letters in 5 days, against 7% at room temperature. This agrees in direction with the Stage 4 preference for cold water. Templated ligation is about 100 times slower than monomer extension. A and U letters copy far more slowly than G and C, which supports running the 0/1 alphabet as G·C (optimistic) and A·U (pessimistic) cases.

## 2. Early-Earth timeline

| Constraint | Value | Source |
|---|---|---|
| Oldest zircon; earliest sign of liquid water | 4,404 ± 8 Myr; oxygen isotopes indicate interaction with a liquid hydrosphere ("earliest evidence for continental crust and oceans") | Wilde, Valley, Peck & Graham (2001), *Nature* 409, 175–178, doi:10.1038/35051550 (abstract) |
| Possible biogenic carbon | graphite in a 4.10 ± 0.01 Ga zircon, δ¹³C = −24 ± 5‰; the authors say "not all light carbon isotopic signals are biogenic" and stress the limits of a single grain | Bell, Boehnke, Harrison & Mao (2015), *PNAS* 112(47), 14518–14521, doi:10.1073/pnas.1517557112 |
| Age of LUCA | about 4.2 Ga (4.09–4.33 Ga) | Moody et al. (2024), *Nature Ecology & Evolution* 8, 1654–1666, doi:10.1038/s41559-024-02461-1 (abstract) |
| Window for RNA formation in one impact scenario | about 4.36 ± 0.1 Ga, about 120 ± 100 Myr after a ~10²³ kg impactor that temporarily reduced the atmosphere | Benner et al. (2020), *ChemSystemsChem* 2, e1900035, doi:10.1002/syst.201900035 (arXiv:1908.11327; preprint read) |
| Emergent land before 3.5 Ga | probably no exposed continental crust; mainly hotspot volcanic islands, about 8 × 10⁵ km² in total if there were about 50 as today (authors' estimate) | Bada & Korenaga (2018), *Life* 8(4), 55, doi:10.3390/life8040055 |
| Impact melting during the late heavy bombardment (about 3.9 Ga) | about 1.5–2.5 vol.% of the upper 20 km of crust melted, only 0.3–1.5% molten at any time; about 5–10% of the surface covered by melt sheets deeper than 1 km | Abramov, Kring & Mojzsis (2013), *Chemie der Erde* 73(3), 227–248, doi:10.1016/j.chemer.2013.08.004 (abstract) |
| Global sterilisation | not reached even in the most extreme bombardment simulated (10 times the base case) | Abramov & Mojzsis (2009), *Nature* 459, 419–422, doi:10.1038/nature08015; result read only in the University of Colorado press release of 20 May 2009, not in the paper |

**What this changes.** (1) Land was scarce. Pond and lake settings, which the Stage 4 results favour, existed mainly on volcanic islands, so N_s(t) for these classes is small and must come from this kind of estimate. (2) The record gives a test with a narrow window: water by about 4.4 Ga and LUCA by about 4.2 Ga leave about 200 Myr or less for the whole climb to L4, if LUCA's age is right.

## 3. Still missing

- Duplex dissociation rates as a function of length and temperature for 4–12-mers (measured, not simulated).
- Temperature of the Ding et al. (2022) and Prywes et al. (2016) experiments.
- Primary sources for isocyanide or other activating agents on the early Earth, and their likely concentrations.
- Abramov & Mojzsis (2009) read in the paper itself.
