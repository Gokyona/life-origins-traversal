# Stage 1e: the engine on an established model (Kauffman autocatalytic sets)

`raf.py` is the binary polymer model of Kauffman (1986) in the formulation of Hordijk & Steel (2004), with an exact maxRAF test batched on the GPU. "Life" here means level L1: a self-sustaining autocatalytic (RAF) network exists. Seven conditions: maximum polymer length n (5–10), catalysis level f (0.1–10, log scale), food length t (1–3), fraction of reactions present (0.3–1), fraction of molecules able to catalyse (0.1–1), and two distractors.

Check against theory first (P(RAF) vs f; the 50% point should sit near f ≈ 1–2):

```bash
python raf.py --device xpu
```

Benchmark (reuses `../stage1b/bench1b.py`; about 1 ms per instance on the B70):

```bash
python bench_raf.py --device xpu --surrogate es --fast --budget 64000 --reps 5 --n-test 4000 --gt-seeds 64 \
  --methods random_1seed,sobol_1seed,engine3,engine4_soft --out results_raf.json
```
