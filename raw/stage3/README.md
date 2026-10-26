# Raw results of stage 3

Full outputs from the benchmark machine (Intel Arc Pro B70), uploaded 25 October 2026 and checked against the md5 sums recorded on that machine. Compact result files used in the reports are in [stage3/](../../stage3/): `maps_l4_compact.json`, `maps_sp_compact.json`, `combos_sp_summary.json`, `s3e_joint_results.json` (= s3e_joint.json), `s3e_straddle_results.json` (= s3e_straddle.json).

| file | step | bytes | md5 |
|---|---|---|---|
| `combos_sp.csv.gz` | 3 side (combinations) | 3381110 | `09e2c621551e1e0101b732c906bb4660` |
| `combos_sp.log` | 3 side (combinations) | 862 | `d135cbabc5bc6e49dea06734688565d0` |
| `combos_sp_L3L4.csv.gz` | 3 side (combinations) | 530060 | `b8369fd718d34008dcd4e1463eee99e0` |
| `combos_sp_summary.json` | 3 side (combinations) | 8595 | `279769a05ad55c14ab1ca604f241dff2` |
| `maps_l4.json` | 3c | 160731 | `9e96a095218577fc84ff4b2e104ae574` |
| `maps_l4.log` | 3c | 1301 | `c6ca6ed232aa123ff64a0325681d4def` |
| `maps_sp.json` | 3d | 331171 | `49eb6916f6ea1ffdcfb7024457fb5282` |
| `maps_sp.log` | 3d | 2705 | `b6c8f46c1635668dc565732d34761446` |
| `s3a.json` | 3a | 6031 | `5e4cc75413b6c96e237c9068b95ee2c1` |
| `s3a.log` | 3a | 3243 | `451127e4073768a2cf3ec5b703ed3f1e` |
| `s3b.json` | 3b | 15976 | `421467728b5afed8198d6aaef82c82c2` |
| `s3b.log` | 3b | 5896 | `5e7e51c0c18ca5def706916f9a2fe395` |
| `s3e_joint.log` | 3e | 3161 | `bdcac761d1701eba1bd58881c548bc07` |
| `s3e_joint_test.json` | 3e | 1467 | `a4d30853cacb82e5290e04e14eb51dd5` |
| `s3e_joint_test.log` | 3e | 816 | `a753afdbf5ddc257982238c183c869ef` |
| `s3e_straddle.log` | 3e | 4159 | `0d7b490b5e694b7deaaaa2aa4f6ce76d` |
| `s3e_straddle_test.json` | 3e | 909 | `d1d595c3ca397fb1d577ec9a85ffcd6a` |
| `s3e_straddle_test.log` | 3e | 1429 | `42a3dc970fa855b66fede033f9c1274e` |

`combos_sp.csv.gz`: one row per combination (786,432 rows): the ten condition values, P(L1)…P(L4) from the COMPASS maps, and the highest level with P > 0.5. `combos_sp_L3L4.csv.gz`: only the 81,721 rows that reach L3 or L4. The conditions are system properties of the stage 3d model, not environmental conditions; see [the combinations report](../../results/stage3_combos_system_properties.md).

Raw outputs of stages 1 and 2 are still on the benchmark machine and will be added.
