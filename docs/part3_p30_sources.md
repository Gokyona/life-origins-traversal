# Part III, step P3-0: measured rates and early-Earth ranges (first batch)

**Status:** first batch, 30 October 2026, Toronto time · **Design:** [part3_design.md](part3_design.md)

Every value below was read in the original or in an author or repository copy on the date given; where only part of a source could be read, this is said. Values not yet found are listed under "Still missing". PubMed Central pages could not be read from the working environment (browser check), so other hosts were used.

## 1. Copying chemistry

| Quantity | Value | Conditions | Source |
|---|---|---|---|
| Primer extension by imidazolium-bridged dinucleotides, k_obs | 0.35 h⁻¹ (U\*U) to 84 h⁻¹ (C\*G) at 5 mM intermediate; maximum at saturation 5.7–85 h⁻¹ | 100 mM MgCl₂, 200 mM Tris pH 8.0, 12-nt primer, 26-nt template, 12-nt downstream helper; temperature not found in the text read | Ding et al., bioRxiv 10.1101/2021.10.06.463393 |
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

## 2. Early-Earth timeline

| Constraint | Value | Source |
|---|---|---|
| Oldest zircon; earliest sign of liquid water | 4,404 ± 8 Myr; oxygen isotopes indicate interaction with a liquid hydrosphere ("earliest evidence for continental crust and oceans") | Wilde, Valley, Peck & Graham (2001), *Nature* 409, 175–178, doi:10.1038/35051550 (abstract) |
| Possible biogenic carbon | graphite in a 4.10 ± 0.01 Ga zircon, δ¹³C = −24 ± 5‰; the authors say "not all light carbon isotopic signals are biogenic" and stress the limits of a single grain | Bell, Boehnke, Harrison & Mao (2015), *PNAS* 112(47), 14518–14521, doi:10.1073/pnas.1517557112 |
| Age of LUCA | about 4.2 Ga (4.09–4.33 Ga) | Moody et al. (2024), *Nature Ecology & Evolution* 8, 1654–1666, doi:10.1038/s41559-024-02461-1 (abstract) |
| Window for RNA formation in one impact scenario | about 4.36 ± 0.1 Ga, about 120 ± 100 Myr after a ~10²³ kg impactor that temporarily reduced the atmosphere | Benner et al. (2020), *ChemSystemsChem* 2, e1900035, doi:10.1002/syst.201900035 (arXiv:1908.11327; preprint read) |
| Emergent land before 3.5 Ga | probably no exposed continental crust; mainly hotspot volcanic islands, about 8 × 10⁵ km² in total if there were about 50 as today (authors' estimate) | Bada & Korenaga (2018), *Life* 8(4), 55, doi:10.3390/life8040055 |
| Impact sterilisation | Abramov & Mojzsis (2009), *Nature* 459, 419–422, doi:10.1038/nature08015 — citation checked; the numbers have not yet been read | — |

**What this changes.** (1) Land was scarce. Pond and lake settings, which the Stage 4 results favour, existed mainly on volcanic islands, so N_s(t) for these classes is small and must come from this kind of estimate. (2) The record gives a test with a narrow window: water by about 4.4 Ga and LUCA by about 4.2 Ga leave about 200 Myr or less for the whole climb to L4, if LUCA's age is right.

## 3. Still missing

- Duplex association and dissociation rates for short oligomers (needed for the melting step).
- Templated ligation rates of activated oligomers.
- Rates for making activated nucleotides or bridged intermediates under prebiotic conditions.
- Abramov & Mojzsis (2009): fraction of crust melted, global sterilisation.
- Temperature of the Ding et al. experiments, and its final journal version.
- Prywes et al. (2016) values, re-checked for this table.
