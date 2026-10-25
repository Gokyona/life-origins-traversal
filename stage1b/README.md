# Stage 1b: harder benchmark, shared neural map-maker, disagreement-driven refinement

Follows the negative stage 1 result ([results/stage1_benchmark_b70.md](../results/stage1_benchmark_b70.md)).

- `model10.py`   10-condition protocell chemostat (parasites, hydrolysis, wet/dry cycling, 2 distractors); ~11% of conditions produce life
- `surrogate.py` small neural network fitted to (condition, successes, runs); shared by all methods; bootstrap ensemble for the engine
- `bench1b.py`   random / Sobol baselines vs `engine3` (one run per condition, Sobol coverage, refinement where bootstrap networks disagree)

```bash
python bench1b.py --budget 256000 --reps 5 --n-test 8000 --gt-seeds 128 --out results_1b.json
```

## Stage 1c re-check

Stage 1b's baselines did not improve with more data, because the 600-step fit underfits large datasets. `--surrogate es` replaces both the evaluation map-maker and the engine's ensemble with early-stopping fits (training length set by held-out loss; 256 hidden units for evaluation, 128 for the engine's ensemble):

```bash
python bench1b.py --surrogate es --budget 256000 --reps 5 --n-test 8000 --gt-seeds 128 --out results_1c.json
```
