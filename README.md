# life-origins-traversal

**Mapping the conditions under which life can emerge — by ordered, non-repeating traversal of condition space instead of waiting for chance.**

> Status: early research (stage 0/1). Code, results and write-ups are published here as they are produced, including negative results.


## The idea

Most origin-of-life simulations pick one set of conditions, run for a long time, and hope that self-replicating, self-sustaining chemistry appears. If it does not, someone guesses a new set of conditions and tries again. This wastes compute on hopeless conditions, silently re-tests equivalent ones, and never tells us *where* life is possible.

This project asks a different question:

> **Under which conditions does life emerge with high probability, and under which does it never emerge?**

Instead of tracing one lucky history, we try to draw the **phase diagram of life**: the boundary in condition space between chemistry that stays chemistry and chemistry that becomes self-replicating, self-maintaining and compartmentalised. Every set of conditions is tested at most once, in an order designed so that each test removes as much uncertainty as possible.

## Method in brief

1. **Condition space.** Physical parameters (energy flux, temperature, food supply, dilution, catalytic probability, polymerisation/hydrolysis ratio, copying error rate, amphiphile fraction, membrane permeability, environmental cycling, mixing) are reduced to dimensionless groups (e.g. Damköhler number, catalytic density, Eigen error number μL), merging equivalent conditions into one.
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

- [x] Stage 0 — minimal stochastic model with a built-in error threshold, an origin event and a distractor dimension
- [x] Stage 1, first benchmark — ordered engine did not beat random search (see results)
- [x] Stage 1b — revised refinement beats random search on a harder benchmark (see results)
- [x] Stage 1c — confirmed with a fair map-maker: same accuracy with about 8× fewer simulations (see results)
- [x] Stage 1d — engine made 4× faster; soft selection becomes the default; late plateau remains (see results)
- [x] Stage 1e — the engine on an established model (Kauffman autocatalytic sets, RAF): about 2× fewer simulations (see results)
- [ ] Stage 2 — dividing protocells (stochastic corrector population, then membrane coupling and space); maps for life-detector levels L1–L3 ([plan](docs/stage2_plan.md))
- [ ] Stage 3 — knockout experiments, minimal sufficient conditions, cross-check on an independent model
- [ ] Volunteer computing — hand out non-overlapping Sobol index ranges as work units so anyone can contribute GPU time

## Paper

**[COMPASS: ordered, non-repeating mapping of the conditions for life in stochastic origin-of-life models](paper/COMPASS_stage1.md)** (stage 1 report, 25 October 2026). The engine design, all benchmarks so far, and a check against the known threshold for autocatalytic sets.

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

## Contributing

Ideas, critiques and compute are welcome. Open an issue to discuss the model, the testing logic or the life detector. A volunteer-computing client is planned; until then, running `bench.py` with different seeds and posting the output in an issue already helps.

## License

Code: MIT. Documents and figures: CC BY 4.0.
