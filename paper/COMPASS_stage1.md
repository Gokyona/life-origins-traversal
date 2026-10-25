# COMPASS: ordered, non-repeating mapping of the conditions for life in stochastic origin-of-life models

**Stage 1 report** · Gokyona Wahiro (pen name), independent researcher, Toronto · 25 October 2026

Code and data: <https://github.com/Gokyona/life-origins-traversal> · Text: CC BY 4.0 · Code: MIT

---

## Abstract

Simulations of the origin of life usually run a model for a long time and wait for a rare event. We take a different route. We treat "does life emerge under these conditions?" as a question about a map: a model's condition space is divided into a region where self-sustaining, self-replicating chemistry appears and a region where it never does, and we want the boundary between them. We present **COMPASS** (Condition-space Ordered Mapping by Prioritised Active Sampling of Stochastic systems), an engine that tests conditions in an ordered sequence that never repeats a condition. It starts with low-discrepancy (Sobol) coverage, runs each condition once, builds a shared neural map of the success probability, and spends most of its budget where an ensemble of bootstrap maps disagrees. That is where data are missing, not where the outcome is random by nature. None of these parts is new on its own. Our contribution is to combine them for this question, together with a benchmark protocol that compares every method under an identical map-maker against brute-force ground truth, and an honest record of what failed along the way.

On a 10-condition stochastic protocell model, COMPASS reaches with 32,000 simulations the accuracy that the best baseline (Sobol sampling) needs 256,000 for: about 8× fewer simulations. On Kauffman's autocatalytic sets in the Hordijk–Steel binary polymer model, our GPU implementation first reproduces the published catalysis threshold for the appearance of RAF sets (to within 2% for maximum polymer length n = 7–10, the range covered by the published fit). There, COMPASS needs about 2× fewer simulations than Sobol sampling for the same accuracy, and the threshold read off its map with 32,000 simulations lies within 7% of the directly simulated value on average across 15 slices of the condition space. The saving depends on how thin the boundary is, and the engine's own compute cost matters when simulations are cheap. We describe these limits, and how the engine will be used on expensive protocell models and in volunteer computing.

## 1. Introduction

How the first living system arose from non-living chemistry is still open. Computational work has mostly followed one of two styles: long stochastic simulations of a chosen model, hoping to observe an origin event, or theory that identifies conditions such as the error threshold for replication (Eigen 1971) or the level of catalysis at which autocatalytic sets appear (Kauffman 1986; Hordijk & Steel 2004; Mossel & Steel 2005; Hordijk, Kauffman & Steel 2011).

Our starting point is a simple observation. For a given model, the useful result is often not one trajectory but the **set of conditions** under which life-like behaviour appears, a "phase diagram of life". Mapping it by brute force is costly: the condition space has many dimensions, each simulation is random, and the interesting conditions lie in a thin region near the boundary. We therefore ask how a computer should *order* its tests so that every simulation adds as much as possible to the map and no condition is ever tested twice.

This is an active-learning problem for a stochastic binary outcome. The tools are well known: low-discrepancy designs (Sobol 1967; Owen 1998), query by committee (Seung, Opper & Sompolinsky 1992), ensembles as uncertainty estimates (Lakshminarayanan, Pritzel & Blundell 2017), separating model uncertainty from intrinsic randomness (Houlsby et al. 2011), and active learning of threshold boundaries and level sets (Bryan et al. 2005; Gotovos et al. 2013). For a survey see Settles (2009). We do not claim a new algorithm. What we add is:

1. an engine, COMPASS, that combines these parts for origin-of-life condition spaces, with design choices learned from failed versions;
2. a benchmark protocol in which all methods share the same map-maker, are scored against brute-force ground truth on held-out conditions, and are repeated five times, so that differences come from *which conditions were simulated*;
3. a test on an established model with a known answer, Kauffman's autocatalytic sets, where we check both the simulator and the map against published results;
4. a complete record, including two negative or corrected results.

## 2. The COMPASS engine

### 2.1 Problem

A model has d conditions, scaled to the unit cube u ∈ [0,1]^d. One simulation at u returns a binary outcome (life or no life) and is random, so each condition has a success probability p(u). The goal is to classify every condition as life (p > 0.5) or no life with as few simulations as possible. Accuracy is measured on held-out conditions whose p is known from brute force. We report it on all conditions and on the **hard set**: conditions with 0.1 ≤ p ≤ 0.4 or 0.6 ≤ p ≤ 0.9, which are random but still decidable.

### 2.2 Components

**R1. Archive with canonical keys (never repeat).** Each tested condition is quantised to a canonical key (a 4096-level grid per dimension) and stored. A candidate whose key is already in the archive is never simulated again.

**R2. Ordered coverage.** Candidates come from one scrambled Sobol sequence of 2^19 points. The engine starts with the first 2,048 points, and in every later round 30% of the batch continues along the sequence. This keeps coverage of the whole space and protects against a wrong map.

**R3. One run per condition.** Every condition is simulated once. Repeated runs per condition were consistently wasteful in our tests (Section 5.1): the map-maker averages neighbouring single runs better than repetition does.

**R4. A shared map-maker.** A small neural network (two hidden layers of 256 units, SiLU) is trained on (condition, outcome) pairs with a binary cross-entropy loss to give p̂(u). Training stops early on a 10% held-out split, keeping the best weights, so training length adapts to the amount of data (Section 5.2 shows why this matters).

**R5. Refinement where knowledge is missing.** Five networks (128 units), each trained on a bootstrap resample, are evaluated on the whole Sobol pool. Most of each batch (70%) goes to untested candidates where the five predictions disagree most. Disagreement is high where more data would change the map; it is low where the outcome is simply random (p near 0.5 in a well-sampled region). An earlier rule that targeted p̂(1−p̂) chased such intrinsically random regions and lost to random search (Section 5.1).

**R6. Soft selection.** Instead of taking the strict top of the disagreement ranking, refinement points are drawn with probability proportional to disagreement squared. This spreads them along the boundary rather than piling them into the single most uncertain spot (Section 5.3).

**R7. Growing batches on the GPU.** Each round adds max(1,024, N/3) points, where N is the number of points so far, so the number of retraining rounds grows only logarithmically with the budget (12–17 rounds here). Simulations run as one batched tensor computation on the GPU. The five bootstrap networks are trained together as one batched network, which made the engine 4× faster (Section 5.3).

### 2.3 Algorithm

```
pool    <- scrambled Sobol sequence (2^19 points); archive <- {}
test(first 2048 pool points)                      # R1-R3, one run each
while runs < budget:
    batch <- max(1024, runs/3)
    take 30% of batch: next untested pool points in sequence order      # R2
    fit 5 bootstrap networks on all (u, outcome)                         # R4-R5
    score each untested pool point by the std of the 5 predictions
    draw 70% of batch with probability ~ score^2                         # R6
    simulate new points once each; add their keys to the archive         # R1, R3
final map <- early-stopped network on all data                           # R4
```

## 3. Benchmark protocol

- **Ground truth.** 4,000–8,000 random held-out conditions, each simulated 64–128 times by brute force on the GPU.
- **Methods.** `random_1seed` (uniform random, one run each), `sobol_1seed` (Sobol sequence, one run each), `sobol_8seed` (eight runs each), and COMPASS variants (`engine3`: strict top-k refinement; `engine4_soft`: the default, soft selection).
- **Scoring.** Each method produces an ordered list of simulated conditions. At fixed budgets (2k, 4k, ..., 256k simulations) the same map-maker (R4) is fitted to the prefix of that list and scored on the held-out set. All methods therefore differ only in *which* conditions they simulated.
- **Repetitions.** 5 independent repetitions per method; tables give mean ± standard deviation.
- **Hardware.** One Intel Arc Pro B70 (32 GB) with PyTorch 2.14 (XPU backend), a 6-core desktop CPU.

## 4. Test problems

**Protocell chemostat (6 and 10 conditions).** A well-mixed flow reactor with food F, a master replicator R and mutants or parasites M, simulated by Poisson tau-leaping (Gillespie 2001) for 800 steps. R appears spontaneously from food at a low rate, copies itself with fidelity q = e^(−μL), so the Eigen error threshold is built in (Eigen 1971), and is lost to dilution and temperature-dependent hydrolysis. Mutants copy with relative fitness σ and act as parasites; the 10-condition version adds wet–dry cycling and two distractor conditions with no effect. A run counts as life only if R ends at least 30 copies above five times the abiotic background that spontaneous formation alone would sustain. About 19% (6-D) and 11% (10-D) of conditions produce life.

**Autocatalytic sets (RAF model, 7 conditions).** Kauffman's binary polymer model in the formulation of Hordijk & Steel (2004). Molecules are binary strings up to length n; each split of a string ab into a + b is one reversible reaction (ligation and cleavage); the food set is all strings up to length t; each molecule catalyses each reaction independently with probability p, and the catalysis level f = p × (number of reactions) is the mean number of reactions a molecule catalyses. "Life" means the network contains a RAF set: a reflexively autocatalytic, food-generated set of reactions. This corresponds to self-sustaining metabolism, the first level of our life detector. RAF existence is decided exactly with the Hordijk–Steel maxRAF reduction, which we batch on the GPU. Conditions: n (5–10), f (0.1–10, log scale), t (1–3), fraction of possible reactions present (0.3–1), fraction of molecules able to catalyse (0.1–1, with f held fixed), and two distractors. One simulation costs about 1 ms, roughly 1,000× more than the chemostat.

## 5. Results

### 5.1 First engine: a negative result (6 conditions)

The first engine used repeated runs with a sequential test, a kNN map and refinement by p̂(1−p̂). It did **not** beat random search (budget 128k, 5 repetitions). The best methods at every budget were the simplest: one run per condition, at random or along a Sobol sequence (hard-set accuracy at 128k: `sobol_1seed` 0.892, `random_1seed` 0.885). The refining engines were worst, with a dip to 0.68 on the hard set at 4k runs. Two lessons shaped COMPASS: repeated runs per condition waste budget (R3), and uncertainty scores of the form p(1−p) chase intrinsically random regions (R5). Monotonicity-based inference, also tried, saved runs but did not improve accuracy and is not part of COMPASS. Full tables: [stage1_benchmark_b70.md](../results/stage1_benchmark_b70.md).

### 5.2 Disagreement-driven refinement and a corrected first estimate (10 conditions)

With disagreement-driven refinement and a shared neural map trained for a fixed 600 steps, the engine first appeared to win by 8–12 points on the hard set at every budget (stage 1b). The baselines, however, did not improve at all from 2k to 256k runs. Diagnosis showed that the fixed-length fit underfits large datasets. With early stopping (R4), every method improves a lot (e.g. `random_1seed` at 256k: hard set 0.702 → 0.944), and the engine's advantage becomes a saving in simulations rather than a permanent gap (Table 1). [Stage 1b](../results/stage1b_benchmark_b70.md), [stage 1c](../results/stage1c_benchmark_b70.md).

**Table 1.** 10-condition protocell model, accuracy on the 666 hard conditions (of 8,000 held out), mean ± s.d. over 5 repetitions. Early-stopping map-maker; values from the stage 1d run with the batched map-maker.

| method | 4k | 8k | 16k | 32k | 64k | 128k | 256k |
|---|---|---|---|---|---|---|---|
| sobol_1seed | 0.807 ± 0.021 | 0.826 ± 0.010 | 0.843 ± 0.012 | 0.866 ± 0.014 | 0.903 ± 0.011 | 0.920 ± 0.010 | 0.942 ± 0.009 |
| COMPASS, strict top-k | 0.811 ± 0.018 | 0.849 ± 0.014 | 0.916 ± 0.010 | 0.942 ± 0.011 | 0.959 ± 0.004 | 0.961 ± 0.006 | 0.953 ± 0.004 |
| **COMPASS (default)** | **0.830 ± 0.020** | **0.874 ± 0.012** | **0.918 ± 0.006** | 0.938 ± 0.008 | 0.956 ± 0.009 | 0.962 ± 0.006 | **0.960 ± 0.002** |

On all 8,000 conditions, COMPASS at 32k (0.988–0.989) matches `sobol_1seed` at 256k (0.989), and on the hard set 0.938–0.942 at 32k matches 0.942 at 256k: **about 8× fewer simulations for the same map.** An independent earlier run (stage 1c, unbatched map-maker) gave the same picture: 0.989 / 0.947 at 32k vs 0.990 / 0.948 at 256k (all / hard).

### 5.3 Engine speed and soft selection

Profiling showed that engine time went almost entirely into training small networks, not into simulation: 256k chemostat simulations take under a second on the GPU. Each training step was bound by kernel-launch overhead (2.3–3.6 ms per step regardless of data size). Training the five bootstrap networks as one batched network (R7) cut the engine's time per repetition from about 510 s to about 128 s with unchanged accuracy; bfloat16 and `torch.compile` gave no gain. Of two fixes tested for the plateau after 64k runs, soft selection (R6) was the clearer improvement: +1.9 and +2.5 points on the hard set at 4k and 8k, and the smallest spread at 256k (Table 1). The plateau near 0.96 remained for all variants, which points to the map-maker or the reference answers rather than to point selection. [Stage 1d](../results/stage1d_benchmark_b70.md).

### 5.4 Autocatalytic sets: checking the simulator, the engine and the map against known results

#### 5.4.1 The simulator reproduces the published threshold

Hordijk, Kauffman & Steel (2011) report that, with food strings up to length 2, the catalysis level at which RAF sets appear with probability 0.5 grows only slowly with n. A linear fit to their simulations for n = 7–20 gives f_A(n) = 1.0970 + 0.0189 n, i.e. "between one and two reactions" per molecule. We measured the same threshold by direct simulation for t = 2, all reactions present and every molecule a potential catalyst: 400 instances at each of 49 log-spaced values of f between 0.5 and 10, with the 50% crossing and a 95% interval from the crossings of p ± 2 s.e.

**Table 2.** Catalysis level f* at which P(RAF) = 0.5, food up to length 2.

| n | this work, direct simulation (95% interval) | Hordijk, Kauffman & Steel (2011), f_A(n) | difference |
|---|---|---|---|
| 5 | 0.50 (< 0.50–0.63) | outside fitted range (n = 7–20) | — |
| 6 | 1.05 (0.97–1.11) | outside fitted range | — |
| 7 | 1.21 (1.16–1.24) | 1.229 | −1.6% |
| 8 | 1.25 (1.23–1.27) | 1.248 | −0.1% |
| 9 | 1.28 (1.26–1.29) | 1.267 | +0.6% |
| 10 | 1.31 (1.29–1.32) | 1.286 | +1.6% |

Within the range where the published fit applies, our simulator agrees with it to within 2%, and the threshold rises slowly with n as reported (slope 0.0189 per unit of n in their fit). For n ≤ 6 the threshold is lower; the published fit was not made there, and small networks behave differently. The transition also sharpens with n (at n = 10, P(RAF) is 0.01 at f = 1.06, 0.28 at f = 1.28 and 0.98 at f = 1.54), in line with the literature.

The checks also showed two properties of our condition space worth stating. The fraction of molecules able to catalyse (with f held fixed) has no measurable effect (at n = 8, f = 1.3: P = 0.74, 0.72, 0.71, 0.76 for fractions 1, 0.5, 0.2, 0.1), so in practice it is a third distractor. The fraction of reactions present partly overlaps with f, but also has its own effect: with f × fraction fixed at 1.3, P = 0.76, 0.50, 0.23 for fractions 1, 0.65, 0.4.

#### 5.4.2 Benchmark on the RAF model

Ground truth: 4,000 held-out conditions × 64 simulations. Life fraction 0.418; hard set 932 conditions (23% of the space, versus 8% in the protocell model).

**Table 3.** RAF model, accuracy (mean ± s.d., 5 repetitions).

| method | | 4k | 8k | 16k | 32k | 64k |
|---|---|---|---|---|---|---|
| random_1seed | all | 0.941 ± 0.004 | 0.948 ± 0.002 | 0.958 ± 0.006 | 0.966 ± 0.002 | 0.970 ± 0.003 |
| sobol_1seed | all | 0.945 ± 0.003 | 0.952 ± 0.004 | 0.957 ± 0.003 | 0.966 ± 0.002 | 0.970 ± 0.003 |
| **COMPASS (default)** | all | **0.949 ± 0.003** | **0.959 ± 0.003** | **0.966 ± 0.002** | **0.970 ± 0.002** | **0.973 ± 0.001** |
| random_1seed | hard | 0.859 ± 0.020 | 0.880 ± 0.008 | 0.914 ± 0.019 | 0.939 ± 0.007 | 0.952 ± 0.008 |
| sobol_1seed | hard | 0.869 ± 0.008 | 0.891 ± 0.012 | 0.912 ± 0.009 | 0.939 ± 0.009 | 0.955 ± 0.007 |
| **COMPASS (default)** | hard | **0.876 ± 0.012** | **0.914 ± 0.010** | **0.937 ± 0.005** | **0.955 ± 0.006** | **0.964 ± 0.002** |

COMPASS at 16k matches Sobol at 32k (0.966 / 0.937 vs 0.966 / 0.939, all / hard), and COMPASS at 32k matches Sobol at 64k (0.970 / 0.955 vs 0.970 / 0.955): **about 2× fewer simulations.** The smaller saving than in Section 5.2 is consistent with the shape of the problem. The uncertain region is broad, and the boundary is spread along f, t and the reaction fraction rather than forming a thin shell, so concentrating tests near it pays off less. [Stage 1e](../results/stage1e_raf_benchmark_b70.md).

#### 5.4.3 Does the engine's map recover the known threshold?

The benchmark measures accuracy over the whole space. A sharper test of reliability is whether the map the engine builds gives back a known quantity. For each slice of the condition space (fixed n and t, all reactions present, every molecule a potential catalyst, distractors at 0.5) we read off the catalysis level f* where the map predicts P(RAF) = 0.5, and compare it with the direct-simulation value from Section 5.4.1. Maps were built from the first 8k, 16k and 32k simulations of each method (3 repetitions). Fifteen slices have a crossing inside the grid (t = 1 and 2 for n = 5–10; t = 3 for n = 8–10).

**Table 4.** Threshold f* read off the maps, food up to length 2 (mean ± s.d. over 3 repetitions).

| n | direct simulation (95% interval) | COMPASS 16k | COMPASS 32k | Sobol 16k | Sobol 32k |
|---|---|---|---|---|---|
| 6 | 1.047 (0.968–1.105) | 1.046 ± 0.054 | 1.059 ± 0.059 | 1.060 ± 0.027 | 1.056 ± 0.050 |
| 7 | 1.209 (1.160–1.242) | 1.171 ± 0.027 | 1.185 ± 0.039 | 1.167 ± 0.022 | 1.171 ± 0.001 |
| 8 | 1.247 (1.226–1.270) | 1.220 ± 0.034 | 1.235 ± 0.032 | 1.231 ± 0.046 | 1.227 ± 0.017 |
| 9 | 1.275 (1.255–1.294) | 1.256 ± 0.048 | 1.262 ± 0.028 | 1.272 ± 0.074 | 1.261 ± 0.021 |
| 10 | 1.306 (1.292–1.320) | 1.285 ± 0.059 | 1.280 ± 0.021 | 1.306 ± 0.097 | 1.280 ± 0.019 |

**Table 5.** Mean absolute error of the map threshold over all 15 slices and 3 repetitions, as |log10 f*_map − log10 f*_direct| and the equivalent relative error in f.

| method | 8k | 16k | 32k |
|---|---|---|---|
| Sobol sampling | 0.052 (12.6%) | 0.051 (12.5%) | 0.036 (8.7%) |
| **COMPASS (default)** | 0.054 (13.1%) | **0.037 (9.0%)** | **0.029 (6.8%)** |

With 32k simulations, COMPASS's map places the RAF threshold within 7% of the direct value on average across all slices. At 16k it is as close as Sobol sampling at 32k, the same factor of about 2 seen in accuracy. At 8k the two are equal. On the t = 2 slices, both maps at 32k fall inside or just below the 95% interval of the direct value, and both recover the published result: a threshold between 1 and 2 reactions per molecule that rises slowly with n. The largest errors occur where the threshold lies near the edge of the explored range (t = 3, small n) and, for both methods, the maps at low budgets tend to place the threshold slightly too low. The engine was never told about the threshold or about n, t and f; it found the boundary from single yes/no simulations. The same slices also show how strongly the threshold depends on the food set: with only monomers as food (t = 1) the direct threshold rises from 1.37 (n = 5) to 2.29 (n = 10), while with food up to length 3 small networks already form RAF sets below f = 0.5.

## 6. Discussion

**When COMPASS helps.** The saving was about 8× on a model with a thin life region and about 2× on a model with a broad uncertain region. It was largest at small and medium budgets, and the methods converge when simulations are plentiful. The engine is therefore most useful exactly where it is needed: expensive simulations and limited budgets.

**Cost of the engine itself.** COMPASS spends compute on training its maps: about 2 minutes per run of the benchmark on our GPU, versus 15–90 s for the baselines. With the chemostat (microseconds per simulation) this overhead dominates. With the RAF model (1 ms) it is about equal to the simulation cost. For models where one simulation takes seconds or more, it becomes negligible: at 1 s per simulation, 256k baseline simulations would take about 71 hours and 32k COMPASS simulations about 9.

**Reliability.** We consider a map reliable when (i) the simulator reproduces known results (Table 2), (ii) the map gives back a known quantity (Tables 4–5), (iii) accuracy is measured on held-out conditions with brute-force ground truth, and (iv) results are repeated. All four hold here. The accuracy ceilings near 0.96–0.97 come partly from the map-maker's capacity and partly from noise in the 64–128-run reference answers. We did not separate the two.

**Limitations.** The protocell models are our own toy models. The RAF criterion is structural: it says a network could sustain itself, not that its kinetics would, and pure RAF networks have limited evolvability (Vasas, Szathmáry & Santos 2010; Vasas et al. 2012). All results come from one GPU and one map-maker architecture. We compared against random and Sobol sampling only, not yet against other active-learning methods such as BALD or level-set estimation on the same benchmark. The engine's settings (30/70 split, batch schedule) were chosen by hand.

## 7. Outlook

1. **Expensive protocell models (stage 2).** Compartments that grow and divide, with competing replicators and parasites, in the spirit of the stochastic corrector model. Each simulation will cost orders of magnitude more, which is where fewer, better-placed tests matter most.
2. **Head-to-head comparison** with BALD, level-set estimation and Gaussian-process methods on the same benchmark.
3. **Minimal sufficient conditions.** Once a map exists, knock out conditions one at a time to find which are necessary for life to appear.
4. **Volunteer computing.** Because the Sobol sequence is deterministic and every condition has a canonical key, ranges of Sobol indices can be handed out as non-overlapping work units, checked by redundant computation, and merged into one archive.
5. **Beyond the origin of life.** The same engine applies to any stochastic simulator with a yes/no outcome and an expensive simulation: phase diagrams in materials science, epidemic thresholds, or reliability limits in engineering.

## 8. Reproducibility

All code is in the repository. Results in this report come from the following commands on an Intel Arc Pro B70; per-stage reports give the exact commits.

```bash
# protocell model, stages 1c/1d (stage1b/)
python bench1b.py --device xpu --surrogate es --fast --budget 256000 --reps 5 --n-test 8000 --gt-seeds 128 \
  --methods sobol_1seed,engine3,engine4_anneal,engine4_soft,engine4_both --out results_1d.json
# RAF threshold check and benchmark (stage1e_raf/)
python raf.py --device xpu
python bench_raf.py --device xpu --surrogate es --fast --budget 64000 --reps 5 --n-test 4000 --gt-seeds 64 \
  --methods random_1seed,sobol_1seed,engine3,engine4_soft --out results_raf.json
```

The threshold analysis of Section 5.4.3: `python stage1e_raf/threshold_check.py --device xpu` (defaults: 3 repetitions, budgets up to 32k, 400 direct runs per grid point).

## Acknowledgements

Code, analyses and text were developed with the help of an AI assistant (Claude, Anthropic). All simulations ran on the author's own hardware.

## References

- Bryan, B., Nichol, R. C., Genovese, C. R., Schneider, J., Miller, C. J. & Wasserman, L. (2005). Active learning for identifying function threshold boundaries. *Advances in Neural Information Processing Systems* 18. <https://proceedings.neurips.cc/paper/2005/hash/8e930496927757aac0dbd2438cb3f4f6-Abstract.html>
- Eigen, M. (1971). Selforganization of matter and the evolution of biological macromolecules. *Die Naturwissenschaften* 58(10), 465–522.
- Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems. *Journal of Chemical Physics* 115, 1716. <https://doi.org/10.1063/1.1378322>
- Gotovos, A., Casati, N., Hitz, G. & Krause, A. (2013). Active learning for level set estimation. In *Proceedings of the 23rd International Joint Conference on Artificial Intelligence (IJCAI)*, 1344–1349. <https://www.ijcai.org/Abstract/13/202>
- Hordijk, W. & Steel, M. (2004). Detecting autocatalytic, self-sustaining sets in chemical reaction systems. *Journal of Theoretical Biology* 227(4), 451–461. <https://doi.org/10.1016/j.jtbi.2003.11.020>
- Hordijk, W., Kauffman, S. A. & Steel, M. (2011). Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. *International Journal of Molecular Sciences* 12(5), 3085–3101. <https://doi.org/10.3390/ijms12053085>
- Hordijk, W. & Steel, M. (2017). Chasing the tail: the emergence of autocatalytic networks. *BioSystems* 152, 1–10. <https://www.sciencedirect.com/science/article/abs/pii/S030326471630274X>
- Houlsby, N., Huszár, F., Ghahramani, Z. & Lengyel, M. (2011). Bayesian active learning for classification and preference learning. arXiv:1112.5745. <https://arxiv.org/abs/1112.5745>
- Kauffman, S. A. (1986). Autocatalytic sets of proteins. *Journal of Theoretical Biology* 119, 1–24. <https://www.sciencedirect.com/science/article/abs/pii/S0022519386800479>
- Lakshminarayanan, B., Pritzel, A. & Blundell, C. (2017). Simple and scalable predictive uncertainty estimation using deep ensembles. *Advances in Neural Information Processing Systems* 30. <https://proceedings.neurips.cc/paper/2017/hash/9ef2ed4b7fd2c810847ffa5fa85bce38-Abstract.html>
- Mossel, E. & Steel, M. (2005). Random biochemical networks: the probability of self-sustaining autocatalysis. *Journal of Theoretical Biology* 233, 327–336. <https://www.sciencedirect.com/science/article/abs/pii/S0022519304004825> (preprint: arXiv:q-bio/0406044)
- Owen, A. B. (1998). Scrambling Sobol' and Niederreiter–Xing points. *Journal of Complexity* 14(4), 466–489.
- Settles, B. (2009). Active learning literature survey. Computer Sciences Technical Report 1648, University of Wisconsin–Madison. <https://research.cs.wisc.edu/techreports/2009/TR1648.pdf>
- Seung, H. S., Opper, M. & Sompolinsky, H. (1992). Query by committee. In *Proceedings of the 5th Annual Workshop on Computational Learning Theory (COLT)*, 287–294. <https://doi.org/10.1145/130385.130417>
- Sobol', I. M. (1967). On the distribution of points in a cube and the approximate evaluation of integrals. *USSR Computational Mathematics and Mathematical Physics* 7(4), 86–112. <https://doi.org/10.1016/0041-5553(67)90144-9>
- Vasas, V., Szathmáry, E. & Santos, M. (2010). Lack of evolvability in self-sustaining autocatalytic networks constraints metabolism-first scenarios for the origin of life. *Proceedings of the National Academy of Sciences USA* 107(4), 1470–1475. <https://pmc.ncbi.nlm.nih.gov/articles/PMC2824406>
- Vasas, V., Fernando, C., Santos, M., Kauffman, S. & Szathmáry, E. (2012). Evolution before genes. *Biology Direct* 7, 1. <https://doi.org/10.1186/1745-6150-7-1>

### How references were checked

Each entry was checked on 25 October 2026 against the publisher's page, the authors' own copy, or the proceedings page: titles, authors, venue, volume and pages. Where a source page could not be opened (publisher robot blocks), details were cross-checked against at least two independent reference lists, and a link is given only to a page that was confirmed to belong to that paper. The quantitative claim taken from the literature, the fitted threshold f_A(n) = 1.0970 + 0.0189 n for n = 7–20 with food up to length 2, was read from Hordijk, Kauffman & Steel (2011) directly.
