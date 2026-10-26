# Stage 4 code

Code for [stage 4](../docs/stage4_plan.md), the bottom-up environment search. Each file is a byte-for-byte copy of the version that ran on the benchmark machine (Intel Arc Pro B70), checked by md5.

| file | purpose | md5 |
|---|---|---|
| `s4a_scaling.py` | 4a: COMPASS at 20–30 dimensions ([report](../results/stage4a_scaling_b70.md)) | 99afa784… |
| `chem4c1.py` | 4c-1: polymer chemistry, detectors E1 and E2 ([specification](../docs/stage4c1_spec.md)); version with the Poisson fix used for the published map | 19da968712aa4e669a3c9f41c31f64a1 |
| `run4c1_map.py` | 4c-1 map: ground truth, Sobol and COMPASS, importance, partial dependence, slices, composition | dd78ab510e77e390d0ad1788f636980a |
| `ckpt.py` | atomic saves for crash-safe runs | 66757e944fdb7128768ce9db4acb0aff |
| `thermal.py` | GPU thermal guard (read-only: pauses the job, never changes fans or limits) | 4e2ae639395d9a1cf65387e2cf807f91 |
| `validate4c1.py` | 4c-1 validation tests (a)–(g) | 36ea4eb3b3780aa88889f6e01409f500 |
| `compare_runs.py` | kill-and-resume test: compares an interrupted and an uninterrupted run | b17c20ac33c740f0f5ca35c7752c93ad |
| `kill_resume_test.sh` | runs that test | 91f5a7bacc7be48b28e0797e43aa5a32 |
| `resume_all.sh` | restarts every unfinished job with `--resume` after a power cut (run by hand) | 15b7efa5f3db74f9ba369f3a6a534501 |
| `watchdog.sh` | restarts a queued job that died or stalled, then starts the next | 11e03771ae80c6581d9b5cff82058f39 |
| `report4c1.py` | text report of a 4c-1 map result | fdcf029620071b05e9ee4b0c7f3d40f3 |
| `sens4c1_hydrolysis.py` | hydrolysis-law sensitivity check (pH and Mg²⁺ floors) | 9aba9c99439f74ed3a7d6b27ade6a436 |
| `compare_gt4c1.py` | compares the original and rerun ground truth | 4001af6f50721e5a813aa0e207787f99 |

**Results**

| file | content | md5 |
|---|---|---|
| `map4c1_result.json` | 4c-1 map, published version (ground truth rerun with the fixed code) | c6358fc9afd3880542cb7baa93826e62 |
| `map4c1_result_original.json` | 4c-1 map as first run (first 1,056 ground-truth environments with the code before the Poisson fix) | dd368e1582539a2a49be8072117244e7 |
| `sens4c1_result.json` | hydrolysis-law sensitivity | ba13573fd9851d25fb32fad1d45e63c8 |
| `validate4c1_final.json` | validation tests (a)–(f), final detector code | 0d9e71be8f9788d1cbba14ce5de5e8bb |
| `validate4c1.json` | earlier validation run, includes cost test (g) | 2ee04becd6a7300dba1391298f6da833 |

`run4c1_map.py` imports `bench1b.py` and `surrogate.py` from [stage1b](../stage1b/); set `LOT_STAGE1B` to that folder.
