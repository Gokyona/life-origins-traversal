# Stage 4 code

Code for [stage 4](../docs/stage4_plan.md), the bottom-up environment search. Each file is a byte-for-byte copy of the version that ran on the benchmark machine (Intel Arc Pro B70), checked by md5.

| file | purpose | md5 |
|---|---|---|
| `s4a_scaling.py` | 4a: COMPASS at 20–30 dimensions ([report](../results/stage4a_scaling_b70.md)) | 99afa784… |
| `chem4c1.py` | 4c-1: polymer chemistry, detectors E1 and E2 ([specification](../docs/stage4c1_spec.md)) | 0c9d1264faaa4ba3de0ab51ac163225f |
| `run4c1_map.py` | 4c-1 map: ground truth, Sobol and COMPASS, importance, partial dependence, slices, composition | dd78ab510e77e390d0ad1788f636980a |
| `ckpt.py` | atomic saves for crash-safe runs | 66757e944fdb7128768ce9db4acb0aff |
| `thermal.py` | GPU thermal guard (read-only: pauses the job, never changes fans or limits) | 4e2ae639395d9a1cf65387e2cf807f91 |
| `validate4c1.py` | 4c-1 validation tests (a)–(g) | 36ea4eb3b3780aa88889f6e01409f500 |
| `compare_runs.py` | kill-and-resume test: compares an interrupted and an uninterrupted run | b17c20ac33c740f0f5ca35c7752c93ad |
| `kill_resume_test.sh` | runs that test | 91f5a7bacc7be48b28e0797e43aa5a32 |
| `resume_all.sh` | restarts every unfinished job with `--resume` after a power cut (run by hand) | 97065bbf43171188bf14618407f685c7 |

`run4c1_map.py` imports `bench1b.py` and `surrogate.py` from [stage1b](../stage1b/); set `LOT_STAGE1B` to that folder.
