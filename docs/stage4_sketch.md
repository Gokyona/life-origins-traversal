# Stage 4 sketch: from environmental conditions to life

**Status:** early sketch, 25 October 2026. It is not yet a plan, and no runs have been done. Stage 3 will be finished first.

## Why

The stage 2–3 models start from vesicles that already contain replicators. Their conditions mix two kinds of quantity. Some are environmental (food supply, degradation as a stand-in for dilution, leak as membrane permeability, spatial versus well-mixed). Others are properties of the chemical system itself (division size S, parasite advantage σ, gene types K, error rate μ, membrane synthesis km, heritable variation η). Those maps therefore answer "which chemical systems become alive", not "which environments give rise to life". The second question is the one the project set out to answer.

The original environmental conditions are:

| | condition |
|---|---|
| P1 | energy flux |
| P2 | temperature |
| P3 | food (monomer) supply |
| P4 | dilution |
| P5 | number of monomer types |
| P6 | catalytic probability |
| P7 | polymerisation / hydrolysis ratio |
| P8 | copying error rate |
| P9 | amphiphile fraction |
| P10 | membrane permeability |
| P11 | environmental cycling (wet–dry, hot–cold) |
| P12 | mixing strength |

## Idea: two layers

Layer A maps the environment onto model quantities. Layer B is the stage 3 protocell model, or a bottom-up replacement where a mapping is not enough. Only P1–P12 are inputs. Quantities such as S, km and μ become intermediate results.

| condition | how it enters | type |
|---|---|---|
| P1 energy flux | monomer activation: food production and the condensation side of polymerisation | mapping (scenario constant for yield) |
| P2 temperature | Arrhenius factors on replication, degradation, membrane turnover and leak; error rises with T; vesicles fail above the lipid phase transition | mapping |
| P3 food supply | sets food directly | mapping |
| P4 dilution | adds to degradation of free molecules and gives a cell washout rate (new term) | mapping |
| P5 monomer types | molecular diversity in the RAF layer; template alphabet (information per position vs per-position error) | partly mapping, partly bottom-up |
| P6 catalytic probability | whether and how fast a metabolism exists: the stage 1e RAF model inside each cell | bottom-up |
| P7 polymerisation / hydrolysis | chain-length distribution sets template length L and so genome error 1 − (1 − ε)^L; whether templates emerge at all | partly mapping, partly bottom-up |
| P8 copying error rate | per-position ε → μ through L | mapping |
| P9 amphiphile fraction | whether vesicles form at all (critical vesicle concentration), membrane growth and vesicle size | bottom-up (self-assembly) |
| P10 membrane permeability | leak of templates and uptake of food | mapping |
| P11 environmental cycling | periodic rates (as in stage 1b); cycle-driven vesicle fusion and division | partly mapping, partly bottom-up |
| P12 mixing | spatial structure: neighbourhood, exchange between cells, well-mixed vs grid (Péclet / Damköhler numbers) | mapping |

## Possible staging

- **4a:** use the mapping layer only, over the mappable conditions, with the bottom-up conditions (P5, P6, P7, P9, the division part of P11) held at fixed favourable values. This gives an environment-driven map quickly.
- **4b:** replace the fixed values with bottom-up layers one at a time, starting with P6 through the existing RAF model.

## Caveat

Every mapping brings unknown constants (activation energies, yields). The environment map will therefore be a family of maps under stated assumptions. COMPASS can treat the uncertain constants as extra conditions, which makes the sensitivity explicit. All physical relations and parameter values will be checked against published sources before use.
