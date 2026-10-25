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

## Design check against the literature

- Molecules, food set (strings up to length t) and one reversible reaction per split match the binary polymer model of Hordijk & Steel; the number of reactions is (n−2)·2^(n+1)+4.
- Detection is the exact Hordijk–Steel maxRAF reduction (closure of the food set, then remove reactions that are uncatalysed or unreachable, repeat).
- Theory check (t = 2, all reactions present, every molecule a potential catalyst, 200–400 instances per cell): P(RAF) = 0.5 at f ≈ 1.1 (n = 6), 1.25 (n = 8), 1.3 (n = 10), and the transition gets steeper with n, in line with the reported 1–2 reactions per molecule.
- Catalysts per reaction are Poisson (approximating binomial) with at most 10 slots; truncation probability is below 2·10⁻⁴.

Things to keep in mind when reading results:

- `c` (fraction of molecules able to catalyse, with f held fixed) has no measurable effect: at n = 8, f = 1.3, P = 0.74 / 0.72 / 0.71 / 0.76 for c = 1 / 0.5 / 0.2 / 0.1. In practice it is a third distractor.
- `dens` is partly confounded with f (f counts all possible reactions) but also has its own effect: with f·dens fixed at 1.3, P = 0.76 / 0.50 / 0.23 for dens = 1 / 0.65 / 0.4.
- Food length matters strongly: P = 0.10 / 0.69 / 0.91 for t = 1 / 2 / 3.
- A RAF is a structural property: it says a network could sustain itself, not that its kinetics would. Pure RAF networks also have limited evolvability (Vasas et al. 2010, 2012), so this model corresponds to level L1 of our life detector only.
