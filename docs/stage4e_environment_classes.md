# Stage 4e: which real settings match the environments at the top of the ladder?

**Status:** first version, 29 October 2026, Toronto time · **Plan:** [stage4_plan.md](stage4_plan.md), step 4e · **Inputs:** the confirmed environments of [4c-3](../results/stage4c3_heredity_b70.md), [4c-4](../results/stage4c4_vesicles_b70.md) and [4c-5](../results/stage4c5_coupling_b70.md)

This step adds no simulation. It compares the environments where the model reaches heredity, vesicles and coupled protocells with settings described in the literature. Every source was checked on the date given; values not read in the original are marked.

## 1. The environments at the top

The 4 environments found by the search where heredity and filled vesicles coexist, and where coupled protocells (L3) appear, are engine12538, engine11230, engine15233 and engine15539.

| Variable | Range at the 4 points | Notes |
|---|---|---|
| Temperature | 0–25 °C, cycling by about 10–12 °C across the duplex melting point | three of four at 0–7 °C |
| Cycle period | 63–992 h (2.6–41 days) | from multi-day weather to monthly or seasonal |
| Residence time τ | 11,000–63,000 h (1.3–7 years) | closed or slowly flushed water |
| pH | 4.4–7.6 | only engine12538 (7.6) is inside the measured fatty-acid vesicle range |
| Salinity | 0.39–0.83 M | similar to sea water (about 0.47 mol/kg Na⁺) |
| Mg²⁺ | 0.35–1.0 mM at three points; 39 mM at engine11230 | sea water about 53 mM |

The two cold 4c-3 heredity points that hold membranes, #1558 (7.6 °C, pH 8.0) and rep1_1926 (7.8 °C, pH 6.7), fit the same pattern.

## 2. Comparison with real settings

| Setting | What the literature gives | Match with the top environments |
|---|---|---|
| **Volcanic hot-spring pools** (wet–dry cycling at the edges) | Damer & Deamer (2020) propose hot-spring pools with wet–dry cycles. Milshteyn et al. (2018) tested water from Yellowstone: Bison Pool (pH 7.9, 68 °C, Na⁺ 12.5 mM, Mg²⁺ 0.08 µM) and a Midway Geyser Basin pool (pH 3.3, 53 °C). Fatty-acid vesicles formed in both | **Partial.** The low salt and Mg²⁺ suit membranes, and cycling is present. But the water is far hotter than any confirmed point (0–25 °C), and the hot pools are much fresher than the points found (0.4–0.8 M). In the model, hot water above the duplex melting point blocks copying. A cooler pool fed by such water, at the outer edge of a hydrothermal field, would come closer |
| **Alkaline hydrothermal vents** | Lost City fluids: 40–90 °C, pH 9–11, dissolved H₂ up to 15 mmol/kg (Kelley et al. 2005); vents as early reactors (Martin et al. 2008) | **No match.** No confirmed heredity point lies above pH 8, and vent flow implies short residence times, while the top points need years |
| **Closed-basin, carbonate-rich lakes** | Toner & Catling (2020): closed-basin lakes accumulate phosphate up to about 0.1 molal at steady state, more on evaporation. Brines are buffered to pH 6.5–9 under a CO₂-rich early atmosphere; modelled at 0–50 °C; seasonal wet–dry cycles at lake margins | **Good on several axes.** Closed basins give long residence times, temperatures can be cool, pH reaches the lower end of our range (6.5–7.6), and salinity is high. Their brines can be far saltier than our range (up to 1.2 M) |
| **Cold and icy settings** (freeze–thaw) | Attwater et al. (2010): ice promotes the activity of an RNA polymerase ribozyme, protects it from hydrolysis and enables "exceptionally long replication products", through concentration in the eutectic phase | **Partial, and outside the model.** Three of the four points sit at 0–7 °C, so the model favours cold water. But the model has no ice and no eutectic concentration; its coldest setting is liquid water at 0 °C. The match is in direction, not in mechanism |
| **Open sea water** | Milshteyn et al. (2018): no tested amphiphile formed stable vesicles in sea water. They attribute this to 10 mM Ca²⁺ and 54 mM Mg²⁺ forming insoluble soaps at pH 8.1, and to the osmotic effect of about 0.5 M NaCl | **Conflict.** The top points have sea-water-like salinity (0.39–0.83 M). They allow membranes only because the model's salt tolerance (sal_crit, up to 3 M) and Mg²⁺ tolerance (up to 50 mM) are layer C ranges. In fatty-acid experiments, salt at this level destroys vesicles |

## 3. Reading

- **The model points to cool, slowly flushed, enclosed water bodies with long temperature cycles.** These are lakes or ponds with seasonal or multi-day cycles, not vents or hot pools.
- **The closest described setting is a closed-basin lake,** cooler than the hot-spring scenario and fed slowly enough that molecules live for years. Carbonate-rich lakes (Toner & Catling 2020) match the residence time, the near-neutral to slightly acidic pH and the cycling.
- **Two conflicts must be stated in the paper:**
  1. Salt. The top points are as salty as sea water, where measured fatty-acid vesicles fail (Milshteyn et al. 2018). Either fresher versions of these environments also reach heredity, which the maps do not yet show, or the result depends on membranes more salt-tolerant than fatty acids.
  2. Acidity. Three of the four points need a membrane with pKa_app near 6, at the edge of the modelled range and below measured fatty-acid and alcohol mixtures (from about pH 6.5; Jordan et al. 2019).
- **The cold preference agrees in direction with work on ice** (Attwater et al. 2010), but the model has no freezing, so this is a hint, not a test.

## Open items

- Check whether fresher variants of the top environments (salinity below 0.1 M) keep heredity, using the existing map data (CPU) or a small set of confirmation runs, if the owner wants it.
- Read the full text of Attwater et al. (2010) for temperatures and the eutectic concentrations.
- Verify volume and pages of Toner & Catling (2020).

## References

- Attwater, J., Wochner, A., Pinheiro, V. B., Coulson, A. & Holliger, P. (2010). Ice as a protocellular medium for RNA replication. *Nature Communications* 1, 76. doi:10.1038/ncomms1076. Checked 29 October 2026 against the RePEc record and the MRC LMB group page (abstract only).
- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452. (Checked for the [condition review](condition_review.md).)
- Jordan, S. F., Rammu, H., Zheludev, I. N., Hartley, A. M., Maréchal, A. & Lane, N. (2019). Promotion of protocell self-assembly from mixed amphiphiles at the origin of life. *Nature Ecology & Evolution* 3, 1705–1714. doi:10.1038/s41559-019-1015-y. Checked 28 October 2026.
- Kelley, D. S. et al. (2005). A serpentinite-hosted ecosystem: the Lost City hydrothermal field. *Science* 307(5714), 1428–1434. Checked 26 October 2026.
- Martin, W., Baross, J., Kelley, D. & Russell, M. J. (2008). Hydrothermal vents and the origin of life. *Nature Reviews Microbiology* 6, 805–814. (Checked for the condition review.)
- Milshteyn, D., Damer, B., Havig, J. & Deamer, D. (2018). Amphiphilic compounds assemble into membranous vesicles in hydrothermal hot spring water but not in seawater. *Life* 8(2), 11. doi:10.3390/life8020011. Checked 29 October 2026 in PubMed Central (PMC6027054): waters, ion data, results and stated causes.
- Toner, J. D. & Catling, D. C. (2020). A carbonate-rich lake solution to the phosphate problem of the origin of life. *PNAS*, doi:10.1073/pnas.1916109117 (approved 27 November 2019). Checked 29 October 2026 against the authors' PDF: phosphate, pH, temperature range and wet–dry cycles. Volume and pages are not yet verified.
