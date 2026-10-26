# Stage 4a: ranges of the 13 environmental variables, with sources

**Status:** first version, 26 October 2026 (02:30 Toronto time) · **Plan:** [stage4_plan.md](stage4_plan.md), Section 3 · **Scope:** liquid water, carbon chemistry, surface or sea floor of an Earth-like planet

Each range must cover the real settings that origin-of-life scenarios discuss, without being so wide that most of the map lies in impossible conditions. Every anchor below was checked against the source named, on the date given. Where no source fixes a bound, we say so and give our reasoning ("argued").

| variable | range | anchors (checked) | basis |
|---|---|---|---|
| temperature | 0–150 °C | Lost City alkaline vent fluids about 40–90 °C (Kelley et al. 2005). Yellowstone spring waters 85–94 °C, where water boils at about 92 °C (USGS Yellowstone Volcano Observatory). The upper limit of known life is 122 °C, for *Methanopyrus kandleri* strain 116 grown under high pressure (Takai et al. 2008). | The upper bound of 150 °C is argued: it lies above the known limit of life, so the map can show where chemistry fails even with liquid water. Above 100 °C, water is liquid only under pressure (next row). |
| pressure | 1–1,100 bar | Hydrostatic pressure rises by about 1 bar per 10 m of water; the deepest ocean, about 11 km, gives about 1,100 bar. | Physics. Pressure decides whether water at a given temperature is liquid: saturation vapour pressure from a standard fit. |
| pH | 2–12 | Acid springs: pH 1.0–2.5 at Monument Geyser Basin, 2.6–6.1 at Norris, near neutral at Beryl Spring. Alkaline pools: pH 9–10 at Imperial Geyser (USGS YVO). Lost City vent fluids: pH 9–11 (Kelley et al. 2005). | The range covers acid springs to alkaline vents. The lowest extreme (pH ~1) is left out and noted. |
| redox (H₂ availability) | log scale, 6 decades below the vent maximum | Lost City fluids carry dissolved H₂ from below 1 to 15 mmol/kg (Kelley et al. 2005). | The upper end is anchored at vent values; the lower end (oxidised settings) is argued. |
| ultraviolet radiation | none to full unshielded surface flux (log, 4 decades) | On a planet with a CO₂ atmosphere, CO₂ blocks wavelengths below about 204 nm and water below about 168 nm, so the surface receives broadband UV from about 200 nm upward (Ranjan & Sasselov 2016). | The deep sea floor receives none. The absolute flux scale enters through the UV yield and damage constants (layer C). |
| ion composition | salinity 0–1.2 M; Mg²⁺ 10⁻⁴–0.1 M (log) | Sea water contains about 0.47 mol/kg Na⁺, 53 mmol/kg Mg²⁺ and 10 mmol/kg Ca²⁺ (lecture table, University of Hawaiʻi; the primary reference composition is Millero et al. 2008, whose individual values we have not yet opened). Fatty-acid vesicles formed in Yellowstone hot-spring water but in none of the tests in sea water (Milshteyn et al. 2018). | The range covers fresh hot-spring water through sea water to brines. |
| feedstock composition and supply | supply rate on a log scale (4 decades, relative); composition as fractions of nucleotide-like, amino-acid-like and amphiphile precursors | — | Argued. No source fixes absolute prebiotic supply rates; the relative scale sets the ratio of supply to loss, which is what the chemistry depends on. |
| dilution and flow | residence time 1 h – 10⁵ h (log) | — | Argued: from flowing vent channels (hours) to closed pools (years). |
| mixing and geometry | exchange rate between reaction sites (log), from isolated pores to well-mixed | — | Argued. |
| concentration mechanisms | evaporation factor 1–10³ (log); thermal-gradient factor 1–10⁸ (log) | Thermal gradients in mm-sized pores concentrate single nucleotides more than 10⁸-fold (Baaske et al. 2007). | The thermal-gradient range is anchored at the experiment; the evaporation range is argued. |
| cycling | depth 0–1 (fraction of water lost in the dry phase); period 1–10³ h (log) | Wet–dry cycling in hot-spring pools drives polymerisation and encapsulation (Damer & Deamer 2020). | The period range is argued: from geyser splash cycles to day–night and seasonal cycles. |
| mineral surfaces | availability (log, 3 decades) and type (none, clay, metal sulphide) | Transition-metal sulphides are discussed as early catalysts at vents (Martin et al. 2008). | Argued. |
| available time | 1–10⁶ years (log) | — | Argued. It enters as the ratio of available time to the slowest time scale; see the [4c-1 specification](stage4c1_spec.md). |

## Open items

- Open the primary sea-water reference composition (Millero et al. 2008) to replace the lecture-table values.
- Find measured ranges of prebiotic feedstock supply where they exist, for example cyanide production in the cyanosulfidic scenario (Patel et al. 2015, volume and pages still to be confirmed).
- Find measured residence times for vent channels and hot-spring pools.

## References

- Baaske, P. et al. (2007). Extreme accumulation of nucleotides in simulated hydrothermal pore systems. *PNAS* 104(22), 9346–9351. (Checked for the [condition review](condition_review.md).)
- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452. (Checked for the condition review.)
- Kelley, D. S. et al. (2005). A serpentinite-hosted ecosystem: the Lost City hydrothermal field. *Science* 307(5714), 1428–1434. Checked 26 October 2026 against a copy of the paper: fluid temperature, pH and dissolved H₂.
- Martin, W., Baross, J., Kelley, D. & Russell, M. J. (2008). Hydrothermal vents and the origin of life. *Nature Reviews Microbiology* 6, 805–814. (Checked for the condition review.)
- Millero, F. J., Feistel, R., Wright, D. G. & McDougall, T. J. (2008). The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale. *Deep-Sea Research I* 55(1), 50–72, doi:10.1016/j.dsr.2007.10.001. Citation checked 26 October 2026 in the IAPWS seawater release and the TEOS-10 publication list; ion values not yet read from the paper (the TEOS-10 Manual table of the composition, Appendix D, could not be retrieved).
- Milshteyn, D., Damer, B., Havig, J. & Deamer, D. (2018). Amphiphilic compounds assemble into membranous vesicles in hydrothermal hot spring water but not in seawater. *Life* 8(2), 11. Checked 26 October 2026 against the University of Minnesota research record.
- Ranjan, S. & Sasselov, D. (2016). Influence of the UV environment on the synthesis of prebiotic molecules. *Astrobiology* 16(1). Checked 26 October 2026 against the arXiv record (1511.00698) and the author's 2017 Harvard dissertation.
- Takai, K. et al. (2008). Cell proliferation at 122 °C and isotopically heavy CH₄ production by a hyperthermophilic methanogen under high-pressure cultivation. *PNAS* 105(31), 10949–10954. Citation checked 26 October 2026 in the Japanese National BioResource Project reference record; the 122 °C limit, the strain and its origin (Kairei field, Central Indian Ridge) in the JAMSTEC press release of 29 July 2008.
- University of Hawaiʻi, OCN 623 lecture "Salinity" (2017), table "Major ions in seawater of salinity 35". Checked 26 October 2026.
- USGS Yellowstone Volcano Observatory, "Complicated plumbing: hot springs and steam vents in Yellowstone National Park". Checked 26 October 2026.
