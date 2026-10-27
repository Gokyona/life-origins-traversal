# Stage 4 code

Code for [stage 4](../docs/stage4_plan.md), the bottom-up environment search. Each file is a byte-for-byte copy of the version that ran on the benchmark machine (Intel Arc Pro B70), checked by md5.

| file | purpose | md5 |
|---|---|---|
| `s4a_scaling.py` | 4a: COMPASS at 20–30 dimensions ([report](../results/stage4a_scaling_b70.md)) | 99afa784… |
| `chem4c1.py` | 4c-1: polymer chemistry, detectors E1 and E2 ([specification](../docs/stage4c1_spec.md)); version with the Poisson fix used for the published map | 19da968712aa4e669a3c9f41c31f64a1 |
| `run4c1_map.py` | 4c-1 map: ground truth, Sobol and COMPASS, importance, partial dependence, slices, composition | dd78ab510e77e390d0ad1788f636980a |
| `ckpt.py` | atomic saves for crash-safe runs | 66757e944fdb7128768ce9db4acb0aff |
| `thermal.py` | GPU thermal guard (read-only: pauses the job, never changes fans or limits), GPU lock (one GPU job at a time) and CPU temperature check; 4c-1 ran with the earlier version 4e2ae639 (in git history) | 3cebc44c4dec70327cdd9696ccf07f8e |
| `validate4c1.py` | 4c-1 validation tests (a)–(g) | 36ea4eb3b3780aa88889f6e01409f500 |
| `compare_runs.py` | kill-and-resume test: compares an interrupted and an uninterrupted run | b17c20ac33c740f0f5ca35c7752c93ad |
| `kill_resume_test.sh` | runs that test | 91f5a7bacc7be48b28e0797e43aa5a32 |
| `resume_all.sh` | after a power cut, (re)starts the watchdog in its memory-limited scope; the watchdog then resumes the queue one job at a time (run by hand; the 4c runs used earlier versions 15b7efa5 and 07a4a602, in git history) | ee47bc0dfdeb25db60872c29dc2264f7 |
| `watchdog.sh` | restarts a queued job that died or stalled, then starts the next; each job in its own memory-limited scope (version with the 4c-3 map in the queue; earlier version c1cf591c in git history) | c7f9af1818e68666aaf459751f27028c |
| `temp_watch.sh` | overnight temperature watch (read-only on hardware): pauses the stage-4 jobs at 90 °C GPU or CPU, resumes below 75 °C | df348296eeaf3ae11dfcd40871697957 |
| `preflight_mem.py` | memory pre-flight run before every long job: peak GPU and host memory, growth and speed per batch size on worst-case point sets | 2e70b90740438c21547f7e0e157b9c19 |
| `report4c1.py` | text report of a 4c-1 map result | fdcf029620071b05e9ee4b0c7f3d40f3 |
| `sens4c1_hydrolysis.py` | hydrolysis-law sensitivity check (pH and Mg²⁺ floors) | 9aba9c99439f74ed3a7d6b27ade6a436 |
| `compare_gt4c1.py` | compares the original and rerun ground truth | 4001af6f50721e5a813aa0e207787f99 |
| `chem4c2.py` | 4c-2: catalysis and self-sustaining networks (level L1) on top of the 4c-1 chemistry ([specification](../docs/stage4c2_spec.md)) | f5eefabca737fd78c163bff96807a302 |
| `cat4c2.py` | 4c-2 catalysed ligation among the abundant species of a site, RAF of a site | 13a08af5aad16b5b0acee0ac8029f8f8 |
| `raf4c2.py` | 4c-2 fixed-hash catalysis rule and sparse maxRAF | ecfdb020747bccac1ec07843b6e65a72 |
| `run4c2_map.py` | 4c-2 map (E1, E2, L1; COMPASS only), optional focus box and detector settings | d45d1f57345333aa1ccdbb2ca4c699b6 |
| `validate4c2_t1.py` | 4c-2 validation test 1: RAF threshold | 767d2d31c1956c02295d49d231fa2780 |
| `validate4c2.py` | 4c-2 validation tests 2–4 | 915ac103b7a2bb71d72763cf1e9c7b57 |
| `cost4c2.py` | 4c-2 cost per run | 1845bba37d05088f490d82ea42644566 |
| `sens4c2_detector.py` | 4c-2 detector-sensitivity check (RAF sites, refresh, catalysis start) | 7441be7cb2c15781db1cf3128bc97b1c |
| `chem4c3.py` | 4c-3: template copying and heredity (level L2) on a two-letter alphabet ([specification](../docs/stage4c3_spec.md)) | cab66933cc8fcd93f5c5b3c78f780952 |
| `cat4c3.py` | copy of `cat4c2.py` with deterministic tie-breaking among abundant species (needed for bit-identical kill-and-resume) | 8dba0b4fd493c30c74b46704df97f3b6 |
| `validate4c3.py` | 4c-3 validation tests 1–5 and 7–9 ([results](../results/stage4c3_validation_b70.md)) | 89c05093dd39949c9a3d03453c21ce70 |
| `scan4c3.py` | 4c-3 two-stage scans for a natural L2 reference (full space, `--box` and `--box2` regions) | 9ad2afead5dea884a76607ef29ebbd48 |
| `sens4c3.py` | 4c-3 check whether the layer-C bounds rule L2 out: processivity table q(l) and k_t / primer-length sensitivity at the box reference | 1bb216847409d6c439c0c830022e0adf |
| `s2stats4c3.py` | reproduction number R (long copies per template per lifetime) at the box2 stage-2 points | b63b2247049d3b3a6fea9bdd7556de85 |
| `sensph4c3.py` | acid-branch check (hydrolysis pH floor 5) and copy-length cap check at the acidic L2 points | b53b2013fd2d770e41705f0cdef7b9f8 |
| `make_ref4c3.py` | writes a natural L2 reference point from the box2 scan for validation test 9 | 60b5a5069825bc4a7cb7a2b62ac0e663 |
| `run4c3_map.py` | 4c-3 map (E1, E2, L1, L2; COMPASS only, batched) | 426a8091ecc745e51f7163775ec6865b |
| `compare_runs4c3.py` | 4c-3 kill-and-resume comparison (runs bit-identical, refitted map analysis within 0.002) | 80ea617c45e7c49ca2b818dfb8b75782 |
| `kill_resume_test4c3.sh` | runs that test | db545f96525e76c75f18d32c37ddf287 |

**Results**

| file | content | md5 |
|---|---|---|
| `map4c1_result.json` | 4c-1 map, published version (ground truth rerun with the fixed code) | c6358fc9afd3880542cb7baa93826e62 |
| `map4c1_result_original.json` | 4c-1 map as first run (first 1,056 ground-truth environments with the code before the Poisson fix) | dd368e1582539a2a49be8072117244e7 |
| `sens4c1_result.json` | hydrolysis-law sensitivity | ba13573fd9851d25fb32fad1d45e63c8 |
| `validate4c1_final.json` | validation tests (a)–(f), final detector code | 0d9e71be8f9788d1cbba14ce5de5e8bb |
| `validate4c1.json` | earlier validation run, includes cost test (g) | 2ee04becd6a7300dba1391298f6da833 |
| `validate4c2_t1.json` | 4c-2 test 1 | ea072f73528de93fac49d5de7f94f66d |
| `validate4c2.json` | 4c-2 tests 2–4 | 8f3f3f2362dab1fe9241c731fff4a603 |
| `cost4c2.json` | 4c-2 cost | 80e84c9797ceaafa91ddc31341d9e5e3 |
| `map4c2b_result.json` | 4c-2 map over the full 37-dimensional space | 62d1b9860442446f1d81e6ffdcbf6688 |
| `sens4c2_result.json` | 4c-2 detector sensitivity; adopted arm c | 4267b7c1ee4d4bd28cb0b73c52134c1a |
| `map4c2_focus_result.json` | 4c-2 focused map (p0 0.01–1, κ 10⁴–10⁶, L_cat 3–5) | 106757c697c10dcb921b78746a098733 |

`run4c1_map.py`, `run4c2_map.py` and `run4c3_map.py` import `bench1b.py` and `surrogate.py` from [stage1b](../stage1b/); set `LOT_STAGE1B` to that folder.
