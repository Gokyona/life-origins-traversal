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
- [ ] Stage 1 — prove (or refute) that ordered testing beats random search on a problem with brute-force ground truth
- [ ] Stage 2 — spatial protocell model (membranes, growth, division) and a full condition map
- [ ] Stage 3 — knockout experiments, minimal sufficient conditions, cross-check on an independent model
- [ ] Volunteer computing — hand out non-overlapping Sobol index ranges as work units so anyone can contribute GPU time

## Honest early findings

On a small CPU pilot run, the first engine (several seeds per condition) did **not** beat random search: repeating seeds on the same condition is wasteful, while one seed per condition plus smoothing across neighbouring conditions works better. A revised engine following that lesson is slightly ahead on the hardest (stochastic) conditions. A full benchmark on an Intel Arc Pro B70 is running; results will be posted under `results/`.

## Contributing

Ideas, critiques and compute are welcome. Open an issue to discuss the model, the testing logic or the life detector. A volunteer-computing client is planned; until then, running `bench.py` with different seeds and posting the output in an issue already helps.

## License

Code: MIT. Documents and figures: CC BY 4.0.
