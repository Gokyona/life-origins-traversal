# Stage 4d-seq: how the environment must change along the ladder

**Date:** 29 October 2026, Toronto time · **Plan:** [stage4_plan.md](../docs/stage4_plan.md), step 4d-seq (descriptive, no new simulation) · **Code:** analyze4dseq.py (md5 f0ba428c), on the benchmark machine · **Data:** result4dseq.json (95477997), tables_4dseq.txt (5f927879), overview and per-variable plots; to be added to the repository

## What was done

The existing map runs were pooled, and for each level the conditions under which it is reached were read off. The maps are 4c-1, 4c-2 (b and focus), 4c-3 (full, box2 reps 0 and 1), 4c-4 (full, box3, box3 v2) and the 4c-5 point study. Layer A, the plan's 13 environmental variables, is coded as 19 dimensions, and all 19 are reported. Levels are nested: E1; E2; L1 (with E1); L2_long; E3 strict; joint (E3 strict and L2_long); L3 strict.

**Sample caveats:**
- Ground-truth runs are used wherever a level has at least 10 positive ground-truth runs: E1 (63,535 runs), E2 (18,385), L1 (85) and E3 strict (304).
- L2_long (14 runs), joint (7) and L3 (298 runs at the 10 point groups of 4c-5) rest on COMPASS hits and confirmed points only. They are not a random sample.
- L1's 85 ground-truth runs come from only a few points.
- The "all feasible" baseline mixes full and focused maps.

## Medians of positive runs per level

| Variable | E1 | E2 | L1 | L2_long | E3 strict | joint | L3 |
|---|---|---|---|---|---|---|---|
| Temperature (°C) | 54.8 | 33.8 | 59 | 8.9 | 20.7 | 4.5 | 7.8 |
| pH | 6.84 | 6.48 | 5.64 | 5.17 | 7.32 | 5.85 | 6.73 |
| Feedstock supply (relative) | 0.21 | 0.8 | 3.7 | 5.2 | 3.1 | 4.5 | 6.3 |
| Nucleotide-like share of feedstock | 0.55 | 0.56 | 0.82 | 0.77 | 0.88 | 0.79 | 0.67 |
| Residence time τ (h) | 2,500 | 7,400 | 12,000 | 28,000 | 12,000 | 30,000 | 29,000 |
| Cycling depth | 0.57 | 0.62 | 0.28 | 0.98 | 0.95 | 1.0 | 0.97 |
| Cycle period (h) | 28 | 26 | 102 | 60 | 44 | 236 | 118 |
| Thermal-gradient factor | 3.4×10⁴ | 2.3×10⁴ | 3.6×10³ | 2.5×10⁶ | 3.4×10⁴ | 3.1×10⁶ | 5.1×10⁵ |
| H₂ (relative) | 0.0025 | 0.0051 | 0.041 | 0.10 | 0.0024 | 0.11 | 0.029 |
| UV (relative) | 0.025 | 0.040 | 0.062 | 0.0021 | 0.0024 | 0.0011 | 0.0011 |
| Mg²⁺ (M) | 0.0027 | 0.0025 | 0.0051 | 0.0008 | 0.0017 | 0.0012 | 0.0005 |
| Salinity (M) | 0.61 | 0.58 | 0.33 | 0.73 | 0.62 | 0.83 | 0.74 |
| Mixing between sites | 0.032 | 0.034 | 0.0062 | 0.028 | 0.46 | 0.0014 | 0.087 |
| Mineral availability | 0.030 | 0.028 | 0.16 | 0.048 | 0.026 | 0.068 | 0.046 |
| Available time (h) | 8.2×10⁶ | 9.9×10⁶ | 4.2×10⁸ | 1.5×10⁸ | 1.4×10⁷ | 2.3×10⁸ | 7.3×10⁶ |

Pressure, evaporation factor, amphiphile share and mineral type barely change. The 10th and 90th percentiles are in tables_4dseq.txt.

## Changes from one level to the next

**No conflicts.** For every pair of consecutive levels and every variable, the 10–90% ranges overlap (114 comparisons). The ladder is a gradual narrowing, not a sequence of contradictory environments.

The largest shifts of the median:
- **E1 → E2:** more feedstock, longer residence time, cooler water (55 → 34 °C).
- **E2 → L1:** more H₂, more minerals, more available time, weaker thermal gradients. These rest on few points.
- **L1 → L2_long (the largest step):**
  - cold water (59 → 9 °C) with deep cycling (depth 0.28 → 0.98);
  - strong thermal gradients (3.6×10³ → 2.5×10⁶);
  - little UV;
  - residence times of years.
- **L2_long → E3 strict:** these are different sets of runs, since vesicles alone are found mostly in milder, better-mixed sites (mixing 0.03 → 0.46, 9 → 21 °C).
- **E3 strict → joint:** mixing reverses (0.46 → 0.0014). The joint level needs nearly isolated sites, plus strong thermal gradients and long available time.
- **joint → L3:** little change. L3 is measured at the joint points themselves.

## Reading

- The top of the ladder asks for an environment that is:
  - cold and deeply cycling;
  - slowly flushed;
  - poorly mixed between sites;
  - shielded from UV;
  - supplied with nucleotide-like feedstock.
- No step needs conditions that exclude the previous level. In this model, one stable environment of this kind could carry chemistry through all the levels reached. A change of environment between levels is not required, though it is not ruled out.
- **Caveat:** for the rare levels, the numbers describe the points found, not the whole region where each level occurs. The focused maps were also chosen from earlier maps.
- **Carry-over is not simulated:** whether the products of one level survive a change of conditions to the next.
