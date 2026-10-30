# Side-product inhibition hypothesis (cold-gap residual, metric A) — analysis only, nothing applied

30 Oct 2026, Toronto. Script `inhibit_fit.py`, output `inhibit_fit.json`. Data: Aitken et al. 2023 only (main + SI;
curves digitised by eye, about +-0.2 mM). Literature: one search found no study of inhibition of nonenzymatic primer
extension by AppA, Passerini products or spent reagent (only enzymatic pyrophosphate inhibition).

## Evidence
| Observation (Table 2: >= +1 at 24 h, RT) | Reading |
|---|---|
| Acetaldehyde 47%, glycolaldehyde 41% vs MBA 65%, although Fig. S3 shows as much or more A*A (peak ~4.0 / 3.3 vs 2.9 mM) and far more impurities (~4.6 / 8.0 vs 1.2 mM at 400 min). Without inhibition the copying integral is 1.19 / 0.87 x MBA | **for** (strongest) |
| More nucleotide (15 mM each) 56% < standard 65%, although activation should rise (AppA ~ [N]^2) | for (weak) |
| More catalyst (75 mM) 50%; AppA formation "correlated with the concentration of 1MeI" | for (weak) |
| Fig. 7D acetaldehyde: faster for 72 h, then stalls at 40% vs 60% (MBA); authors suspect "denaturing organic solutes and inhibitory byproducts such as AppA" | for (qualitative) |
| Less MeNC/MBA (100 mM) 47% (Fig. S5A: 32% vs 62% at -3 C) | against a pure-inhibition view: activation dominates |
| 2AI 5 / 15 mM: 49 / 50% (optimum at 10) | ambiguous (A*A vs *A exchange) |
| pH 7.5 / 7.0: 40 / 13% | activation chemistry, not informative |
| -3 C: impurities near zero (Fig. 5), AppA 1.3 mM after 96 h with 40 mM AMP (Fig. 6) | consistent: side products suppressed in the cold |

## Minimal term and fit
r_ext -> r_ext / (1 + I / K_I), I = accumulated side products (mM). Fit on the three RT aldehyde yields (K_I and rate scale):
**K_I = 4.0 mM** (fitted 0.65 / 0.54 / 0.35 vs 0.65 / 0.47 / 0.41; RMS 0.06).
Same K_I, protocols A, B (side products per addition: RT 1.5 mM, -3 C 0.16 mM per 20 mM nucleotide), C (no side products):

| E_ext | A >= +1 24 h, -3 C (obs 0.87) | B >= +7 5 d (obs 0.60) | C pre-activated, +3 C (obs 0.08) |
|---|---|---|---|
| 0 | 0.87 | 0.98 | 0.56 |
| 10 | 0.75 | 0.82 | 0.29 |
| 15 | 0.68 | 0.65 | 0.19 |
| 18.4 | 0.63 | 0.52 | 0.13 |
| 20 | 0.61 | 0.46 | 0.11 |
K_I x0.5 / x2 at E_ext 18.4: A 0.69 / 0.60, B 0.79 / 0.28.

## Conclusion
Inhibition by side products is supported by the aldehyde series. With K_I ~ 4 mM it closes metric B together with C at
E_ext ~ 15-20 kJ/mol. It does NOT close metric A: at E_ext 18-20 the -3 C prediction is 0.61-0.63 vs 0.87, a remaining
per-letter rate factor of ~2.0 (2.55 without inhibition). No single E_ext in 0-20 fits A, B and C with this term.
The A residual stays unexplained by the data at hand (possible causes, untested: first-step-specific effects of the
monomer *C / bridged C*G mix in the CMP+GMP system, which differs from the AMP NMR system used for all A*A curves).
