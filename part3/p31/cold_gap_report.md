# P3-1 test 6 follow-up: the cold-benefit gap (Aitken et al. 2023)

30 October 2026, Toronto time. CPU only (`cold_gap.py`, output `cold_gap.json`). Paper and SI: `lot-p3/literature/` (MD5 sums in `MD5SUMS_aitken2023.txt`, matching the repository's checksums).

## 1. Experimental facts (from the paper and SI)

| Item | Value | Where |
|---|---|---|
| Activation (Fig. 5) | AMP 20 mM, 2AI 10 mM, 2AT/DCI/4-aminopyrimidine 50 mM, MgCl2 30 mM, HEPES 200 mM pH 8, MBA + MeNC 200 mM once, 10% D2O; 31P NMR for 400 min; RT and −3 °C | SI p. 3; Fig. 5 |
| A*A at RT | peak ≈2.9 mM at ≈80 min, falling to ≈1.5 mM at 400 min; *A rises to ≈7.3 mM; AMP plateaus at ≈7.5 mM within ≈2 h | Fig. 5 (read by eye) |
| A*A at −3 °C | rises more slowly to ≈3.4 mM by ≈250 min, flat to 400 min; *A ≈2.6 mM; AMP still falling at 400 min (≈10 mM) | Fig. 5 |
| Authors' explanation | "off-pathway reactions are suppressed at low temperatures"; A*A "significantly more resistant to hydrolysis, persisting at >3 mM over 10 h"; AppA and Passerini products suppressed | main text p. 16147 |
| −3 °C is liquid | "without any freezing or eutectic phase concentration; freezing was suppressed due to the high concentrations of salt and organic solutes" | p. 16146 |
| ≥+1 at 24 h (Table 2) | standard (RT) 65%; −3 °C 87% (MBA), 85% (acetaldehyde). CMP + GMP 10 mM each, 2AI 10 mM, catalysts 50 mM, Mg 30 mM, MBA/MeNC 200 mM once; primer 2 µM, template 3 µM, overhang 5′-CCG-3′ | Table 2; SI p. 6 |
| ≥+7 after 5 days | −3 °C 60% (Fig. 7C: "62% after 6 days"); RT 7% (Fig. S4). NMPs 10 mM each, six helper trimers 500 µM each, 2AI 30 mM, catalysts 60 mM, Mg 30 mM, MBA/MeNC 200 mM added every 24 h; template overhang 5′-GCGCCUCAU | Fig. 7C, S4; SI p. 8 |
| −3 °C system stalls after ~5 d | attributed to denaturing organic solutes and inhibitory byproducts (AppA) | p. 16148 |
| Re-activation at −3 °C (Fig. 6) | AMP 40 mM, 2AI 30 mM; A*A cycles between 1.2 and 6.0 mM; AppA 1.3 mM after 96 h | Fig. 6 |
| Pre-activated control | *N 10 mM, *NNN 500 µM, Mg 30 mM, 24 h: ≥+7 42% at RT vs 8% at **+3 °C**: cold slows copying here | Fig. 7A,B |

Per step: ≥+1 at 87/65 = 1.34; 1.34^7 = 7.8, close to the ≥+7 ratio 60/7 = 8.6. The data are consistent with a per-step rate gain compounding over 7 letters.

## 2. Matched reproduction in the model

Protocol A: first added letter C on G, C\*C on G, K_m 0.39 mM, k_max 39 /h. Protocol B: the same rate law as a homogeneous per-letter Poisson rate. One rate scale s is fitted per protocol so that the RT yield equals the observed value; −3 °C is then predicted. The A\*A supply is the measured Fig. 5 curve (per 10 mM nucleotide = half the 20 mM AMP curve), superposed for every reagent addition, with exponential tails after 400 min:
- RT: 0.124 /h (from Fig. 5);
- −3 °C: ≤ 0.057 /h (Fig. 6, an upper bound because activation is still running), with 0.02 /h as the lower case.

| Scenario | ≥+1 at 24 h, −3 °C (obs 0.87) | ≥+7 at 5 d, −3 °C (obs 0.60) |
|---|---|---|
| RT A\*A, extension term only (E_ext 18.4) | 0.41 | 0.002 |
| Formation only (cold first 400 min, RT tail, E_ext 0) | 0.74 | 0.17 |
| Persistence only (RT first 400 min, cold tail, E_ext 0) | 0.72 | 0.17 |
| Cold A\*A, E_ext 18.4 (from the pre-activated control) | 0.54–0.56 | 0.02–0.04 |
| Cold A\*A, E_ext 0 | 0.78–0.80 | 0.31–0.41 |
| Cold A\*A, E_ext −10 | 0.89–0.91 | 0.69–0.78 |
| Cold A\*A, E_ext 18.4, K_m with ΔH_bind −50 / −100 kJ/mol | 0.63–0.65 | 0.06–0.07 |

E_ext scan (cold A\*A tail 0.057 / 0.02):

| E_ext (kJ/mol) | ≥+1 | ≥+7 |
|---|---|---|
| +20 | 0.51 / 0.54 | 0.02 / 0.03 |
| +10 | 0.65 / 0.67 | 0.09 / 0.13 |
| 0 | 0.78 / 0.80 | 0.31 / 0.41 |
| −10 | 0.89 / 0.91 | 0.69 / 0.78 |
| −20 | 0.96 / 0.97 | 0.94 / 0.97 |

Both metrics are matched at an apparent in-situ E_ext of about −8 to −10 kJ/mol.

A 4-parameter ODE fit of the Fig. 5 curves (reagent loss, activation rate, selectivity, A\*A loss) is also in `cold_gap.json`. It is poorly identified (RMS 0.35–0.48 mM; the cold A\*A loss runs to its bound), so I use it only qualitatively. It agrees that A\*A persistence is the activation-side term that matters, and that the activation rate, reagent loss and selectivity are not large on their own.

## 3. Decomposition, and what is missing

1. **Activated-block persistence and formation (activation side).** In the model these are E_act,h and the T-independent k_act·a. With the measured curves, persistence and formation each give about half of the observed +1 benefit (0.72–0.74 vs 0.87; 0.17 vs 0.60 for +7). Apparent E of A\*A loss: 0.124 /h at 22 °C vs ≤0.057 (likely ~0.02) at −3 °C, i.e. **≥21, likely ~50 kJ/mol**. That overlaps the low end of the spec's 50–100 range; the measured lower bound of ~21 lies below it.
2. **Reagent.** chem_p31 has a constant agent supply a, with no reagent loss and no re-addition schedule. In the experiment the activation burst ends within ~2 h at RT but runs for >7 h at −3 °C. A constant a cannot represent this, and with a constant a the time-matched cold benefit depends on calibration (test 6: ≤2.6×).
3. **Extension.** The pre-activated control gives an apparent **E_ext ≈ +18 kJ/mol** (Poisson rates 6.2 vs 3.7 for 42% vs 8% ≥+7 at 24 h, 22 vs 3 °C). That is the low end of the spec range 20–80. The in-situ data need ≈ −9 kJ/mol. The ≈27 kJ/mol difference is a cold-favoured factor of ≈2.8 per step that **none of the model's temperature terms provides**. Candidates, all absent from the model:
   - suppression of off-pathway chemistry that consumes reagent and nucleotides (AppA, Passerini products);
   - inhibition of copying by those byproducts and solutes at RT;
   - in protocol B, in-situ activation and persistence of the helper trimers.
4. **K_m(T) via binding enthalpy.** Small effect, because c ≈ 1–3 mM is already above K_m (0.39 mM). ΔH_bind −100 kJ/mol adds only ≈0.1 to the +1 yield.
5. **Duplex stability.** In protocol A the primer is a long, stable duplex, so there is no effect. In B, helper-trimer binding is cold-favoured, but the pre-activated control (helpers present) is still slower in the cold.

## 4. Literature values

| Quantity | Value | Source |
|---|---|---|
| A\*A apparent loss, in-situ mixture | 0.124 /h (22 °C); ≤0.057 /h (−3 °C) | Aitken 2023, Fig. 5, 6 (read by eye) |
| A\*A in copying buffer | t½ ≈ 2 h (A\*U, A\*C); 245 min (C\*C) at RT | Ding et al. 2022 (as in the P3-1 spec) |
| Activated monomer hydrolysis | 4.0e-3 /h (2AI-activated C, primer-extension conditions, 25 °C) | Kaddour/Walton kinetic model, PMC6340129 |
| Extension, apparent E | +18 kJ/mol (22 → 3 °C, pre-activated, helpers) | derived here from Aitken 2023 Fig. 7A,B |
| ΔH of bridged-dinucleotide binding, E_a of 2AI-imidazolide hydrolysis | **not found** (searched; no multi-temperature measurements located) | — |

## 5. Proposed minimal model change (not applied to chem_p31.py)

1. **Reagent pool with loss.** Replace the constant k_act·a by an agent pool a(t) supplied at a layer-A rate (or as periodic additions), lost at k_r(T) = k_r,RT·exp(−E_r/R(1/T − 1/296)). Activation proceeds at k_act·a(t), and a fraction y(T) of activations leads to side products (AppA) instead of blocks. Anchor: Fig. 5 (burst ~2 h at RT, >7 h at −3 °C) for k_r,RT ≈ 0.5–1 /h and E_r ≈ 40–70 kJ/mol (rough). This is the main structural gap.
2. **Widen E_act,h to 20–100 kJ/mol.** Aitken's A\*A data give ≥21 kJ/mol.
3. **Widen E_ext to −20…80 kJ/mol, with ~0 to +20 as the measured band** (+18 pre-activated; ≈ −9 effective in situ). State that a negative effective value stands for the unmodelled in-situ factors (item 3 of section 3); it is not a physical activation energy.
4. **Temperature floor.** Lower the T range to −10 °C, with liquid water only above a freezing point T_f = −3.72 °C × salinity (M, ideal NaCl, i = 2): 0 to −4.5 °C over the model's salinity range. Aitken's −3 °C solution stayed liquid because of ~0.8–1 osmolal organic solutes, which the model does not represent. Also allow the wet-phase floor (T − 10d) to reach T_f instead of 0 °C. Without this, cold-activation regimes cannot be sampled at all.

Item 1 is the one that changes the science; items 2 and 3 are range changes; item 4 changes the environment space. None is applied.
