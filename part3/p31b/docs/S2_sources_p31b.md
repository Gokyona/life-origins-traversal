# P3-1b §2 source extraction (read 30 October 2026, 22:15–22:30 Toronto time)

## §2.1 Copying accuracy: measured mismatch / misincorporation data
**Fang, Jia, Xing & Szostak (2025), PNAS 122(21) e2505720122.** Main text and SI, local copies (`lot-p3/literature`, md5s in README_fang2025.md), searched for mismatch, misincorporat, fidelity, error, and the self-pairs A:A, G:G, C:C, U:U, s2U:s2U. Text read from the Europe PMC JATS copy.
- **No measured misincorporation rate, mismatch-extension rate or fidelity anywhere in the main text or SI.**
- The only statements are qualitative:
  - s2U forms a self-pair "energetically equivalent to a canonical A:U base pair". At equal A it is outcompeted by A:s2U, so "excessive misincorporation of s2U is avoided" (citing their ref. 30).
  - Fidelity is "a particular concern" that needs deep-sequencing experiments not yet done.
  - Ref. 16, Heuberger et al. (2015) JACS 137:2769: s2U "enhances the rate and fidelity" of nonenzymatic primer extension. It is cited only.
- So: no measured ε for AS2U in Fang 2025.

**Sources of part3_p30_sources.md:**
- The one fidelity number is Prywes, Blain, Del Frate & Szostak (2016) eLife 5:e17756: "average fidelity 98%". That is the modern alphabet (A, U, G, C) copied with helper oligomers, 100–200 mM Mg²⁺, pH 8–9, as recorded in P3-0.
- It is an average over all mismatch types for the 4-letter alphabet. It is not a per-pair k_mis/k_match for a 2-letter GC or AU alphabet.
- No other mismatch data in those sources.

**Candidates outside the spec's source list (not read; for the reviewer):**
- Fang ref. 30 (s2U:s2U self-pair crystal structure / competition study).
- Heuberger 2015 (s2U fidelity).
- Prywes 2016 per-mismatch data, if its figures give them.

## §2.2 Salt correction: originals read
1. **SantaLucia (1998) PNAS 95:1460–1465** (PMC19045, read in the PMC HTML; equations read from the PMC equation images M13, M14). DNA.
   - Eq. 7: ΔG°37(oligomer, [Na+]) = ΔG°37(unified oligomer, 1 M NaCl) − 0.114 × N × ln[Na+].
   - Eq. 8: ΔS°(oligomer, [Na+]) = ΔS°(unified oligomer, 1 M NaCl) + 0.368 × N × ln[Na+] (e.u., cal mol⁻¹ K⁻¹).
   - N = total phosphates in the duplex / 2 (8-bp duplex without terminal phosphates: 7). ΔH° assumed salt-independent.
   - Fit to 26 duplexes in 0.01–0.3 M NaCl, < 17 bp. SD 0.60 kcal/mol for ΔG; slope SD 0.033. No Mg²⁺ term.
2. **Chen & Znosko (2013) Biochemistry 52:7477** (PMC3859436, author manuscript, read in PMC HTML; equations from the MathJax source). RNA.
   - 18 duplexes in 71, 121, 221, 621 mM and 1.021 M Na+.
   - The SantaLucia form refitted to RNA: ΔG°37(2) = ΔG°37(1) − 0.1016 × N × ln([Na+]₂/[Na+]₁).
   - Recommended (eq. 26, 0.14 kcal/mol): ΔG°37(2) = ΔG°37(1) + (0.324 f_GC − 0.468) ln([Na+]₂/[Na+]₁) + 0.133 (ln²[Na+]₂ − ln²[Na+]₁). This has no length term.
   - The reference is 1.021 M Na+ (standard melting buffer).
3. **Arteaga, Adams, Meyer, Richardson, Hoener & Znosko (2023) Biophys J 122:565** (PMC9941723, author manuscript, read in PMC HTML; equations from the MathJax source). RNA, Mg²⁺ alone (0.5, 1.5, 3.0, 10.0 mM, no monovalent ions), 17 duplexes.
   - Best ΔG form (eq. 12, 0.30 kcal/mol): ΔG°37(Mg) = ΔG°37(1.021 M Na+) + (0.078 f_GC + 0.337) ln([Mg²⁺]/[Na+]) + 0.066 (ln²[Mg²⁺] − ln²[Na+]).
   - Length-dependent alternative (eq. 10, "modified SantaLucia", 0.48 kcal/mol): ΔG°37(2) = ΔG°37(1) − 0.00796 × (N_bp − 1) ln([Mg²⁺]/[Na+]).
   - The authors say mixtures of Mg²⁺ and monovalent ions were not studied for RNA ("Systematic thermodynamic data ... are needed").
4. **Owczarzy et al. (2008) Biochemistry 47:5336** (Mg²⁺ / mixed-ion rule, DNA): **NOT READ.**
   - Not open access. pubs.acs.org shows a bot-verification page, which I did not try to get past. There is no PMC copy.
   - Per the brief, the Mg²⁺/monovalent mixing-rule sub-item stops here.
