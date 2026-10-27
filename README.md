# life-origins-traversal

**Mapping the conditions under which life can emerge — by ordered, non-repeating traversal of condition space instead of waiting for chance.**

> Status: **Part I (stages 0–3): building and validating the COMPASS engine on test models**, complete ([stage 3 paper](paper/COMPASS_stage3.md)). **Part II (stage 4 on): the search over environmental conditions.** Code, results and write-ups are published here as they are produced, including negative results.


## The idea

Most origin-of-life simulations pick one set of conditions, run for a long time, and hope that self-replicating, self-sustaining chemistry appears. If it does not, someone guesses a new set of conditions and tries again. This wastes compute on hopeless conditions, silently re-tests equivalent ones, and never tells us *where* life is possible.

This project asks a different question:

> **Under which conditions does life emerge with high probability, and under which does it never emerge?**

Instead of tracing one lucky history, we try to draw the **phase diagram of life**: the boundary in condition space between chemistry that stays chemistry and chemistry that becomes self-replicating, self-maintaining and compartmentalised. Every set of conditions is tested at most once, in an order designed so that each test removes as much uncertainty as possible.

## Method in brief

1. **Condition space.** Physical parameters are reduced to dimensionless groups (e.g. Damköhler number, catalytic density, Eigen error number μL), merging equivalent conditions into one. After a [review of the condition list](docs/condition_review.md), stage 4 uses 13 primary environmental variables as inputs (temperature, pressure, pH, redox conditions, ultraviolet radiation, ion composition, feedstock composition and supply, dilution and flow, mixing and geometry, concentration mechanisms, wet–dry and hot–cold cycling, mineral surfaces, available time). Copying error rate, catalytic probability, polymerisation/hydrolysis ratio and membrane permeability are computed from them, not set directly.
2. **Hard constraints first.** Regions ruled out by thermodynamics or information theory are removed before any simulation.
3. **Ordered testing rules.**
   - Canonical key + permanent archive: a condition is never simulated twice.
   - Sobol low-discrepancy order: any prefix of the test sequence covers the space evenly.
   - Monotone inference: where "more is never worse" holds, one result settles a whole region (with spot checks).
   - Multi-fidelity funnel: cheap models screen, expensive models confirm.
   - Boundary refinement / active learning: compute is spent where life and non-life meet.
   - Sequential testing over enumerated seeds: outcomes are stochastic, so we estimate a probability, not a yes/no.
   - Knockout experiments: extract the minimal sufficient set of conditions.
4. **Life detector.** Each run is scored automatically on five levels: L0 plain chemistry → L1 self-maintaining metabolism (autocatalytic set) → L2 self-replication (heredity) → L3 protocell (membrane + metabolism + heredity coupled) → L4 Darwinian evolution.

## Repository layout

```
stage1/          minimal model + ordered testing engine + benchmark
  model.py       stochastic protocell-chemostat model, 6 conditions, GPU-batched (PyTorch)
  engine.py      ordered condition-space testing engine
  bench.py       ordered engine vs random search, accuracy against brute-force ground truth
docs/            write-ups (theory, method proposal, results)
results/         benchmark outputs (JSON + logs)
```

## Quick start

Requires Python 3.10+, PyTorch, NumPy and SciPy. The device is picked automatically (Intel XPU → CUDA → CPU).

```bash
pip install torch numpy scipy
cd stage1
# smoke test (about a minute on CPU)
python bench.py --budget 4000 --reps 1 --n-test 300 --gt-seeds 16 --out smoke.json
# full benchmark
python bench.py --budget 128000 --reps 5 --n-test 4000 --gt-seeds 128 --out results.json
```

## Roadmap

### Part I: the engine (stages 0–3)

Stages 0–3 build COMPASS and test it on a series of increasingly hard model systems. The conditions in these models are test conditions: several of them (division size, parasite advantage, gene types) are properties of the chemical system, not of the environment. Their maps validate the engine and the life detector; they do not answer which environments give rise to life.

- [x] Stage 0 — minimal stochastic model with a built-in error threshold, an origin event and a distractor dimension
- [x] Stage 1, first benchmark — ordered engine did not beat random search (see results)
- [x] Stage 1b — revised refinement beats random search on a harder benchmark (see results)
- [x] Stage 1c — confirmed with a fair map-maker: same accuracy with about 8× fewer simulations (see results)
- [x] Stage 1d — engine made 4× faster; soft selection becomes the default; late plateau remains (see results)
- [x] Stage 1e — the engine on an established model (Kauffman autocatalytic sets, RAF): about 2× fewer simulations (see results)
- [x] Stage 2 — dividing protocells: stochastic corrector population, membrane coupling, multi-fidelity test, maps for life-detector levels L1–L3 ([plan](docs/stage2_plan.md), [paper](paper/COMPASS_stage2.md))
- [x] Stage 3 — analysis tools and cross-checks ([plan](docs/stage3_plan.md), [paper](paper/COMPASS_stage3.md))
  - [x] 3a knockouts · [x] 3b minimal sufficient conditions · [x] 3c Darwinian evolution (L4) · [x] 3d independent spatial model
  - [x] 3e engine: joint nested map, comparison with a standard active-learning method (see results)
  - [x] 3f stage 3 paper, closing Part I and stating its limits

### Part II: the environment search (stage 4 on)

- [x] [Review of the condition list](docs/condition_review.md): the original 12 conditions were neither sufficient nor independent; replaced by 13 primary environmental variables (approved 25 October 2026)
- [ ] Stage 4 — fully bottom-up artificial chemistry (monomers → polymers → catalysts and templates → vesicles → protocells, nothing living put in) driven only by the 13 environmental variables, rate laws calibrated against experiments, then COMPASS maps ([plan](docs/stage4_plan.md), [sketch](docs/stage4_sketch.md))
- [ ] Stage 5 — beyond the stage 4 scope: other solvents and planetary or cosmic environments, possibly life not based on carbon
- [ ] Volunteer computing — hand out non-overlapping Sobol index ranges as work units so anyone can contribute GPU time

## Paper

**[COMPASS: ordered, non-repeating mapping of the conditions for life in stochastic origin-of-life models](paper/COMPASS_stage1.md)** (stage 1 report, 25 October 2026). The engine design, the stage 1 benchmarks, and a check against the known threshold for autocatalytic sets.

**[Mapping the conditions for dividing protocells with COMPASS](paper/COMPASS_stage2.md)** (stage 2 report, 25 October 2026). Protocell populations after the stochastic corrector and chemoton ideas, COMPASS benchmarks (about 8× and 4× fewer simulations), a partial negative result for multi-fidelity screening, and maps of three levels of life.

**[What protocells need, when they evolve, and what the COMPASS engine can and cannot tell us](paper/COMPASS_stage3.md)** (stage 3 report, closing Part I, 25 October 2026). Necessary mechanisms, minimal sufficient conditions, Darwinian evolution with a neutral-marker drift control, a cross-check on a spatial model, joint nested maps, a comparison with a standard active-learning rule, and the limits of Part I.

## Results so far

| stage | problem | outcome | report |
|---|---|---|---|
| 1 | 6 conditions, kNN map-maker | **Negative:** the ordered engine did not beat random search; its refinement step chased intrinsically random regions | [stage1_benchmark_b70.md](results/stage1_benchmark_b70.md) |
| 1b | 10 conditions, shared neural map-maker, disagreement-driven refinement | **First positive result:** +1.5 points overall and +8 to +12 points on hard conditions vs random and Sobol, from 4k runs on (5 repetitions). Caveat: baselines plateau, so the map-maker was limiting them; corrected in 1c | [stage1b_benchmark_b70.md](results/stage1b_benchmark_b70.md) |
| 1c | same as 1b, map-maker trained with early stopping (fair to all methods) | **Confirmed, smaller effect:** the engine reaches with 32k runs the accuracy the best baseline needs 256k for (0.989 / 0.947 all / hard), about **8× fewer simulations**; at 256k the methods nearly tie. Corrects the size of the 1b gap | [stage1c_benchmark_b70.md](results/stage1c_benchmark_b70.md) |
| 1d | same as 1c; batched map-maker training; late-round variants | **Engine 4× faster** (about 2 min per repetition). Soft selection (sample by disagreement²) is best: hard set 0.874 vs 0.849 at 8k and the most stable at 256k (0.960 ± 0.002). All engines plateau near 0.96 after 64k | [stage1d_benchmark_b70.md](results/stage1d_benchmark_b70.md) |
| 1e | Kauffman / Hordijk–Steel autocatalytic sets (RAF), 7 conditions, ~1 ms per run; implementation reproduces the known threshold | **Advantage transfers, smaller:** the engine matches Sobol sampling at half the budget (hard 0.955 at 32k vs 0.955 at 64k), about **2× fewer simulations**; the saving depends on how thin the life boundary is | [stage1e_raf_benchmark_b70.md](results/stage1e_raf_benchmark_b70.md) |
| 2a | Stochastic corrector protocell population (Szathmáry & Demeter 1987), 8 conditions, ~0.8 ms per population; reproduces the known behaviour | **About 8× fewer simulations:** COMPASS at 8k equals random and Sobol at 64k (0.987); half the errors at equal budget | [stage2a_scm_benchmark_b70.md](results/stage2a_scm_benchmark_b70.md) |
| 2b | Coupled protocells: membrane made inside the cell triggers division, 10 conditions, level L3 (homeostasis), ~0.9 ms per population | **About 4× fewer simulations:** COMPASS at 16k equals random and Sobol at 64k (0.977 vs 0.975–0.976); about 40% fewer errors at equal budget; hard set 0.79 vs 0.67–0.70 | [stage2b_protocell_benchmark_b70.md](results/stage2b_protocell_benchmark_b70.md) |
| 2d | Multi-fidelity COMPASS: cheap short simulations (7.8× cheaper, 97% agreement) for screening, five variants at equal cost | **Partial negative:** cheap screening helps only at small budgets (8k: 0.969 vs 0.964); from 16k COMPASS on full simulations is as good or better (64k: 0.986 vs 0.983–0.985); default unchanged | [stage2d_multifidelity_b70.md](results/stage2d_multifidelity_b70.md) |
| 2e | Condition maps for life levels L1, L2, L3 on the coupled protocell model, 32k simulations per level | COMPASS beats Sobol at every level (L1 0.969 vs 0.952, L2 0.989 vs 0.978, L3 0.980 vs 0.971); levels nested; thresholds match direct simulation; optimal cell size between parasite control and gene loss made quantitative | [stage2e_maps_b70.md](results/stage2e_maps_b70.md) |
| 3a | Mechanism knockouts of the coupled protocell model (equal split, no compartments, no selection, no membrane coupling, no decay, one gene) | Compartments and selection between cells are essential (no compartments: L1 −45%, L3 impossible; no selection: L2 −11%); random sharing matters only for large cells with strong parasites | [stage3a_knockouts_b70.md](results/stage3a_knockouts_b70.md) |
| 3b | Minimal sufficient conditions: smallest boxes in condition space with ≥ 95% life, confirmed by direct simulation | L1: 4 conditions held (17% of space); L2: 5 (6.9%); L3: 5, with narrow windows for membrane synthesis and cell size (1.9%); parasite advantage never needs holding | [stage3b_minimal_conditions_b70.md](results/stage3b_minimal_conditions_b70.md) |
| 3c | Darwinian evolution (L4): heritable quality q with 8 neutral markers as an internal drift control; map with COMPASS | L4 in 11% of space, all inside L3 (64% of L3 evolves); needs heritable variation η ≥ 0.02, parasite mutation μ < 0.04 and enough food; less food needs more variation; false positives 1%; COMPASS 0.876 vs Sobol 0.815 on hard points | [stage3c_evolution_b70.md](results/stage3c_evolution_b70.md) |
| 3d | Cross-check on an independent spatial model (16×8 grid, local replacement, leak between neighbours), 12 conditions | All 3a–3c conclusions hold; with almost no leak L2/L3 regions about a quarter smaller than well-mixed; compartments must be nearly sealed (leak > 0.7–0.8% per step destroys heredity and evolution); under strong parasites the largest cell that still evolves is about halved; COMPASS beats Sobol (L3 errors 1.7% vs 3.3%) | [stage3d_spatial_b70.md](results/stage3d_spatial_b70.md) |
| 3e | Engine: one joint map for L1–L4 with nesting built in; COMPASS's selection rule against the straddle rule (Bryan et al. 2005) on the 2b model | Joint map: zero contradictions between levels (separate maps 1.1–1.4%), same accuracy, 3× faster training; adopted. Straddle ties with COMPASS at every budget (64k: 0.987 vs 0.987; Sobol 0.976); all active rules need about 4× fewer simulations than Sobol, so COMPASS's value is the whole loop, not its scoring formula | [stage3e_engine_b70.md](results/stage3e_engine_b70.md) |
| 3 side | List of all 786,432 combinations of a 4-level grid over the spatial model's 10 conditions, scored by the COMPASS maps (system properties, not environmental conditions) | 10% of combinations reach L3/L4; direct check agrees 93–98% per level; most robust evolving systems: small cells, one gene type, low error, weak parasites, plenty of food | [stage3_combos_system_properties.md](results/stage3_combos_system_properties.md) |
| 4a | Engine at 20 and 30 dimensions (test model padded with irrelevant or weakly relevant inputs), before the environment search | COMPASS still beats Sobol: about 3.6–4.0× fewer simulations with irrelevant extra dimensions, 2.8–2.9× when many dimensions matter a little (10 dimensions: 4.4×); engine overhead does not grow with dimension | [stage4a_scaling_b70.md](results/stage4a_scaling_b70.md) |
| 4c-1 | First bottom-up layer: feedstock → activated monomers → chains, over 29 dimensions (13 environmental variables plus unmeasured constants), 64k simulations per method | Activated monomers reach working concentration in 44% of the space, persistent chains of ≥ 8 units form in 16%; residence time, temperature (cold better), feedstock and energy decide, minerals and salinity hardly matter; chain composition follows the feedstock, hot water favours peptide-like chains; COMPASS saves about 2× (E1) and 1.6× (E2) | [stage4c1_polymers_b70.md](results/stage4c1_polymers_b70.md) |
| 4c-2 | Catalysis and self-sustaining networks (L1) over 37 dimensions: fixed random catalysis rule, RAF detection with persistence and a paired catalysis-off control; full map and a focused map | Self-sustaining networks never appear in the full space; with highly catalytic short chains they appear in 0.3% of runs, favoured by long residence time, feedstock, UV and short catalysts; networks are small (median 27 reactions), nucleotide-like, raise their products about 1.5× with only 1.5% catalysed joins, and mostly arise without long chains; polymer and mineral catalysis do not help long chains; survived a power cut | [stage4c2_catalysis_b70.md](results/stage4c2_catalysis_b70.md) |
| 4c-3 *(in progress)* | Template copying and heredity (L2) over 42 dimensions, 2-letter alphabet: primed copying with errors, duplexes that must melt, lineage detector against two random-copy controls | Validation done: copy fidelity, strand separation, parasite takeover (short templates win 99.4%) and an error threshold (heredity in 15/16 runs at 1.7% errors per letter, none at 5%) reproduced; several drafted rules that blocked heredity by construction were found and fixed; under the final rules heredity needs temperature cycles faster than molecule lifetime; map pending | [stage4c3_validation_b70.md](results/stage4c3_validation_b70.md) |

## Contributing

Ideas, critiques and compute are welcome. Open an issue to discuss the model, the testing logic or the life detector. A volunteer-computing client is planned; until then, running `bench.py` with different seeds and posting the output in an issue already helps.

## License

Code: MIT. Documents and figures: CC BY 4.0.
