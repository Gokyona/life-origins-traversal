# Stage 4c-4 specification: amphiphiles and vesicles (level E3)

**Status:** specification, drafted 27 October 2026, 21:55 Toronto time · **Plan:** [stage4_plan.md](stage4_plan.md), step 4c-4 · **Builds on:** [4c-3](stage4c3_spec.md) (template copying, L2), with all of 4c-1 to 4c-3 unchanged

The fourth layer adds membranes to the chemistry of 4c-3. Amphiphiles come from the feedstock. Above a critical concentration they assemble into vesicles, and the vesicles enclose whatever chains are nearby when they form. They grow by taking up amphiphiles, divide when gentle shear or cycling acts on grown vesicles, and break down where ions, pH or heat destroy membranes. No rule says "be a cell": a vesicle is a physical container, and what it holds is chance.

The coupling between content and vesicle growth (L3, 4c-5) needs one physical law in this layer: osmotic growth (Chen, Roberts & Szostak 2004). 4c-5 then adds nothing new and only detects.

## 1. Amphiphiles

- **Supply.** The amphiphile-precursor fraction of the feedstock (layer A variable 9, decoded but unused so far) is now used. Precursors become amphiphile monomers A at the synthesis rate E, the same energy-driven step as for the other monomers. A is a single pool per site. Amphiphile chain length and mixture enter only through layer C constants (below).
- **Losses.** Free A washes out at 1/τ and mixes between sites like the monomer pools.
- **Assembly.** Membrane amphiphile in a site is B. Free A above the critical vesicle concentration CVC goes into membrane at rate k_as (fast, 10³ /h) times (c_A − CVC). Membrane returns to free A at the same rate when c_A < CVC. At steady state, without other losses, c_A = CVC and B = total − CVC.
- **CVC.** CVC = CVC₀ · exp(γ_T (T − 25 °C)), where CVC₀ is layer C, 0.01–100 mM, log, and γ_T is layer C, 0–0.05 /°C. Wet–dry cycling concentrates A in the dry phase by the same evaporation factor as other species.

## 2. When membranes are stable

Membranes form vesicles only where the conditions allow bilayers. Outside this window, assembled amphiphile is non-compartment material (oil or micelles): it counts toward B for the mass balance but forms no vesicles, and existing vesicles break down at rate k_break (layer C, 0.1–10 /h, log). The window is defined by:

- **pH:** |pH − pKa_app| ≤ w_pH. Fatty-acid bilayers form near the apparent pKa of the acid in the membrane. pKa_app is layer C, 6–10; w_pH is layer C, 0.3–2.
- **Divalent ions:** Mg²⁺ ≤ Mg_crit. Mg²⁺ destabilises fatty-acid vesicles (fatty-acid salts crystallise). A myristoleic acid / glycerol monomyristoleate 2:1 mixture tolerates up to 4 mM (Mansy & Szostak 2008), and citrate-chelated Mg²⁺ allows RNA copying inside vesicles (Adamala & Szostak 2013). Mg_crit is layer C, 0.5–50 mM, log, standing for the mixture and any chelation.
- **Salinity:** salinity ≤ sal_crit, layer C, 0.03–3 M, log.
- **Temperature:** the hottest phase is ≤ T_mem, layer C, 40–100 °C. Mixed fatty-acid membranes are stable over a wide range (Mansy & Szostak 2008).

## 3. Vesicles

Each site holds up to NV = 8 vesicles. Vesicle j has a membrane amount a_j (amphiphile molecules) and an internal volume v_j = v₀ (a_j / a₀)^(3/2), where a₀ is the membrane of a newly formed vesicle. a₀ and v₀ are fixed so that a newborn vesicle encloses a fraction φ₀ of the site (layer C, 10⁻³–0.3, log). Every chain slot gets a compartment label: 0 = outside, j = inside vesicle j.

- **Formation.** If the window holds, B grows and a vesicle slot is free, new vesicles form at rate k_form (B_free / a₀), with k_form 10 /h. B_free is the membrane amphiphile not yet in vesicles. On formation, each outside chain is enclosed with probability φ_j = v_j / V_site.
- **Rehydration.** When a dry phase ends, every vesicle breaks open and re-forms with the same membrane, re-enclosing chains with probability φ_j × f_rehyd. f_rehyd is layer C, 1–10, and stands for the encapsulation by lamellae on rehydration (Damer & Deamer 2020). The chains a vesicle held before are kept inside; only new ones are added.
- **Growth.** Vesicle j takes up membrane from B_free, and from other vesicles by exchange, at a rate proportional to a_j · (1 + β · Π_j). Π_j is the osmotic load: the enclosed chain mass (letters) divided by v_j, relative to the site's outside chain mass per volume. β is layer C, 0–10. With β = 0, growth is by area only. With β > 0, vesicles with more enclosed chains take membrane from emptier ones, as observed by Chen, Roberts & Szostak (2004). Membrane exchange between vesicles runs at rate k_ex (layer C, 10⁻³–1 /h, log), with net flow towards higher a_j(1 + βΠ_j).
- **Division.** A vesicle with a_j ≥ 2a₀ divides at rate k_div = s_shear · mix_factor + s_cycle · (wet–dry switches per hour). Gentle shear divides grown, filamentous fatty-acid vesicles while they keep their contents (Zhu & Szostak 2009). s_shear and s_cycle are layer C, 0.01–10, log. mix_factor is the layer A mixing variable already used for monomer mixing. Each enclosed chain goes to either daughter with probability ½, and paired strands go together. Membrane splits in half. A daughter needs a free vesicle slot; if none is free, the smaller daughter's contents are released outside.
- **Losses.** Vesicles wash out at 1/τ with their contents. They break down at rate k_break outside the stability window, releasing their contents outside.

**Permeability.** Enclosed chains grow from the site's activated monomers times f_perm (layer C, 0.01–1). Activated nucleotides cross fatty-acid membranes within seconds to minutes, while tetranucleotides are retained for hours (Mansy & Szostak 2008). So chains of length ≥ 4 do not leave a vesicle except by breakdown, washout or division. Copies and catalysed ligations happen only between chains with the same label: a copy inherits its template's label.

New layer C constants (12): CVC₀, γ_T, pKa_app, w_pH, Mg_crit, sal_crit, T_mem, φ₀, f_rehyd, β, s_shear, s_cycle. Two more are fixed for now: k_ex and f_perm, set at the geometric middle of their ranges, with a sensitivity check later. The model then has 54 dimensions.

## 4. Detectors

- **E3 loose:** vesicles persist. In at least half of the sites, at least one vesicle exists at every one of the last 8 checkpoints, and the vesicle population spans ≥ 5 vesicle turnover times (1 / (1/τ + k_break·[outside window] + division-limited loss)).
- **E3 strict:** E3 loose, and the vesicles hold chains: at the last checkpoint, at least half of the sites have a vesicle that encloses at least one chain of length ≥ L_H.
- E3 is linked to nothing: membranes need amphiphiles, not polymers. Overlaps with E1, E2, L1 and L2 are reported.
- **Vesicle lineage record (for L3 in 4c-5).** For each vesicle: id, parent id, generation, birth time, membrane at birth and at division, and the enclosed chain mass and enclosed lineage sequences (as in the L2 detector) at each checkpoint. This adds no new detector here.

Composition fields added: vesicles per site, membrane share of amphiphile, mean vesicle size, enclosed share of chain mass, divisions per vesicle lifetime, maximum generation, share of the L2 lineage mass enclosed, and which came first (E3 or L2).

## 5. Validation before any map

1. **No amphiphiles:** with the amphiphile fraction 0, the run is bit-identical to chem4c3 for the same seed (the layers below are unchanged).
2. **CVC:** in a closed site with a fixed amphiphile amount and no other losses, B equals total − CVC above the CVC and 0 below it, within sampling error, across a scan of amounts.
3. **Stability window:** vesicles appear only inside the window. Moving pH, Mg²⁺, salinity or temperature across a bound breaks existing vesicles at rate k_break.
4. **Encapsulation statistics:** the number of chains enclosed at formation is binomial with p = φ_j (χ² test).
5. **Division:** each daughter's content follows a binomial ½ split (χ² test), duplexes stay together, and membrane is conserved through division.
6. **Osmotic competition (Chen et al. 2004, qualitative):** with β > 0, in a site with half the vesicles loaded with chains and half empty, loaded vesicles gain membrane and empty ones lose it. With β = 0 there is no difference (control).
7. **Growth–division cycle:** with amphiphile supply and shear, vesicle number rises and settles at a steady state where divisions balance losses. Without shear and cycling, vesicles grow but do not divide.
8. **Power cut:** kill and resume gives bit-identical results.
9. **Memory pre-flight** at the worst case (box2 points, all vesicle slots used), under the standing rule.

## 6. Map

A joint map of E1, E2, L1, L2_long and E3 (loose and strict) over the 54 dimensions, COMPASS only, with ground truth 1000 × 16 and a budget of 32k, as in 4c-3. Then a focused map of the region where L2 and E3 can occur together. That region is the overlap of the 4c-3 heredity region with the membrane stability window, and it is the input for 4c-5.

## 7. What this layer cannot show

- The membrane laws are coarse. A stability window and an exponential CVC stand in for real phase diagrams. Ranges stay wide (layer C) until measured values for specific amphiphile mixtures are checked.
- Copying in 4c-3 does not depend on Mg²⁺, while membranes do. The known conflict between Mg²⁺-dependent copying and Mg²⁺-sensitive membranes (Adamala & Szostak 2013) is therefore only half modelled. A Mg²⁺ dependence of k_t would be a 4c-3 revision; it is listed for after this layer.
- Membranes do not catalyse or protect anything chemically: enclosure changes only washout, mixing, permeability and who can copy whom.

## References

- Adamala, K. & Szostak, J. W. (2013). Nonenzymatic template-directed RNA synthesis inside model protocells. *Science* 342(6162), 1098–1100. (Checked for the [condition review](condition_review.md); pages from the Szostak lab's reference list.)
- Chen, I. A., Roberts, R. W. & Szostak, J. W. (2004). The emergence of competition between model protocells. *Science* 305(5689), 1474–1476. (Checked on 27 October 2026 in the Caltech authors repository record.)
- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452.
- Mansy, S. S. & Szostak, J. W. (2008). Thermostability of model protocell membranes. *Proceedings of the National Academy of Sciences USA* 105(36), 13351–13355. (Checked on 27 October 2026 in PubMed Central PMC2528865: Mg²⁺ tolerance up to 4 mM for the 2:1 MA/GMM mixture; activated nucleotides equilibrate within 30 s at 90 °C; tetranucleotides retained over 1 h.)
- Zhu, T. F. & Szostak, J. W. (2009). Coupled growth and division of model protocell membranes. *Journal of the American Chemical Society* 131(15), 5705–5713. (Checked on 27 October 2026 in the Szostak lab's reference list and an author-group page; issue from the latter.)
