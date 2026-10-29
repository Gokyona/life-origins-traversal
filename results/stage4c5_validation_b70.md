# Stage 4c-5 validation: coupling, L3 and L4 detection (Intel Arc Pro B70)

**Status:** validation done 29 October 2026, 04:00 Toronto time; report: [stage4c5_coupling_b70.md](stage4c5_coupling_b70.md). **6 of 8 pre-registered tests pass; tests 3 and 6 fail.** Under the stop rule (spec §3), no map was started. Test 6b, pre-registered after the diagnosis of test 6, passes. · **Specification:** [stage4c5_spec.md](../docs/stage4c5_spec.md), with the revisions and implementation choices in its §6 · **Code:** chem4c5.py (md5 c747dc68), validate4c5.py (4dc1d04f), kill_resume_4c5.sh (958ac576); kept on the benchmark machine and to be added to [stage4/](../stage4/README.md) · **Data:** validate4c5.json (md5 f80f6e5e)

## Before the GPU validation

The spec was revised three times, and the test fixtures once, before any GPU result. Each change came from a CPU test or diagnostic, and each is recorded in spec §6:

1. **Released chains are enclosed again.** In wet–dry cycles, vesicles re-enclose released chains. So original molecules stay in a vesicle family without being copied. A paired random-copy comparison did not remove this effect. Now only chains copied after the window starts count.
2. **The first fixture was closed.** It had no hydrolysis and no washout, so all chain slots were full before the window started. The fixtures now use the full #3709 environment.
3. **Ancestors too far back.** The ancestors were taken about 290 vesicle generations before the end, and copy errors had changed the exact sequence. Inheritance is now judged against a nearby ancestor, with at least 90% sequence identity.

**CPU diagnostic at four real joint points** (2 seeds each, screen fidelity):
- L3 loose: 8 of 8 runs.
- L3 strict: 6 of 8. #3709, engine12538 and engine15233 pass in 2 of 2. engine11230 passes in 0 of 2, with site shares of 0.45 and 0.39 against the threshold of 0.5.
- L4: 0 of 8. Heritability is strong everywhere, but ρ(w, z) is negative at every point.

## Tests

| Test | Result | Detail |
|---|---|---|
| 8. Memory pre-flight (worst case, batch 256) | pass | peak 15.4 GiB VRAM (cap 27.9), host 13.4 GiB (scope 24), growth +0.5% over 5 batches; 1,549–1,595 s per batch; record 422 MiB per batch |
| 1. Bit-identity with chem4c4 | pass | 9 of 9 cases identical (3 seeds at #3709, engine12538, and #3709 without amphiphiles) |
| 2. Record consistency | pass | 0 errors over 78,000 and 123,000 vesicle ids; membrane is conserved exactly through division |
| 3. L3 fixture (β = 5 against β = 0) | **fail** | β = 5: coupling detected in 16 of 16 runs (permutation p = 0.001 in every run), but ρ(g, Π) reached the pre-registered 0.3 in only 3 of 16 (median about 0.27). β = 0: L3 loose in 0 of 16 |
| 4. Inheritance fixture (copying against k_t = 0) | pass | copying: paired inheritance 16 of 16 (site share 0.50–0.81, random-copy share 0). k_t = 0: 0 of 16 |
| 5. False positives at real points (β = 0, 80 runs) | pass | L3 loose or L4 in 0 of 80. Only 31 of the 80 runs passed the L2/E3 gate (none at #3709, C and engine12538) |
| 6. L4 fixture (β = 5 against β = 0) | **fail** | β = 5: L4 in 3 of 16 (criterion ≥ 14), while L3 strict and L3 loose held in 16 of 16. ρ(w, z) median about 0.21. β = 0: 0 of 16 |
| 7. Kill and resume | pass | runs, engine selections, results and accuracy curves identical. The short test run had no L3 or L4 hit, so stored records were not compared (0 of 0 files) |

## Test 6 diagnosis and test 6b

An exact replay of test 6, bit-identical to the saved runs (replay_test6.json, md5 95736855), shows the following:
- heritability passed in 16 of 16 runs;
- every run had at least 17 ancestors;
- selection beyond drift failed in 13 of 16.

The ancestors were the vesicles alive at checkpoint 0 (375 h). By then, the copying founders had already displaced the non-copying ones. All final descendants came from copying founders, so the window saw only variation within one lineage.

Test 6b was written into the spec before it ran (§6). It takes the founders at t = 0. The statistic is the share of final descendants from the copying founders, compared with the maximum of the same share in 16 runs at β = 0 (validate4c5_6b.json, md5 abae6552).

| | Share from copying founders | Pass |
|---|---|---|
| Null, β = 0 (16 runs) | 0.42–0.58, maximum 0.576 | – |
| β = 5 (16 runs) | **1.0 in all 16** (heritability 16 of 16) | **16 of 16** (needed ≥ 14) |
| Second β = 0 set (16 runs, new seeds) | 0.37–0.59 | 1 of 16 (allowed ≤ 2) |

**Test 6b passes.** It shows that selection of a heritable, copied content among vesicle families is detectable when the observation starts before selection has run its course. It does not validate the power of the map's L4 statistic (ρ(w, z) among the vesicles at the window start). That statistic sees only variation within the population present at the window start. Its power there remains untested, and L4 results from it are reported with that caveat.

## Reading

- **The L3 detectors work as intended:**
  - Coupling is detected whenever it is built in, and never without it (tests 3 and 5).
  - Inheritance is detected only when copying happens (test 4).
  - The failure of test 3 lies in the effect-size threshold of 0.3, which was a guess, not in detection.
- **The L4 detector has too little power in its own positive control** (test 6). Heritability and selection beyond drift are both required, and in 13 of 16 runs one of them was not significant. Which one is being checked by an exact replay of test 6, with the per-run details saved.
- **The random-copy control passes L3 loose but not L3 strict.** Osmotic load counts letters whatever their sequence. So L3 loose measures growth coupling through mass, and L3 strict is what carries heredity.

## Next

The next step is decided with the owner. The options are:
- run the point study and focused map for L3 only, with L4 reported as not validated;
- or revise the L4 test with a new pre-registered criterion and validate it again first.
