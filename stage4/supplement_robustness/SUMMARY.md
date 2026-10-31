# P3-1b summary (stopped early on the owner's decision, 31 October 2026, 05:25 Toronto time)

Spec: part3_p31b_spec.md (commit e02a9fe, md5 91d0dbc9). All implementation choices and the reviewer decisions (Reviews 1–10) are in docs/IMPLEMENTATION_P31B.md. Final code: chem_p31b 44f30af0, chem4c4b 61d9a4f4, chem4c5b 6abf1b2b. Base code chem_p31 cfc1e3db is unchanged.

## What was done and what was not

| Item | Status |
|---|---|
| §1 MULTI_START implementation | done (switch; off = chem_p31 bit for bit) |
| §1 tests 1–4 + 2e (sampler) | PASS (validate_p31b.json d22f26dc) |
| §1 test 5 kill/resume | PASS (5 + 5 batches, 160 runs, IDENTICAL) |
| §1 test 6 memory pre-flight, batch 256 | PASS (worst peak reserved 14.43 GiB of the 22.3 GiB limit; host 4.82 GiB; growth 0.5%; 2.1–4.7 s/run with controls) |
| Test 1 reruns after every code change (s3 switches, R6 sampler, Review 7, Review 9) | all PASS (0 differing fields) |
| R6 sampler (Cornish–Fisher start, exact inversion) | PASS: 117,440,512 draws per module identical to the previous sampler; 2e(ii) at λ = 0.02…1e4 pass; max 1 search step |
| §2.1 measured mismatch data | none in Fang 2025 (main + SI) or the P3-0 sources; only Prywes 2016 "average fidelity 98%" for the 4-letter alphabet |
| §2.2 salt | sources read: SantaLucia 1998, Chen & Znosko 2013, Arteaga 2023. Owczarzy 2008 not readable, so the Mg/Na mixing rule is a recorded gap. SALT_CORR decided (Review 4) but **not implemented or validated** before the stop |
| §3 L2 ablations at the 7 confirmed L2 points | done: base, A0, A1, A2, A2a, A4, A6 + diagnostics A2b–A2e, records, A5 at 3 points |
| A5 at engine435 / engine5704 | not testable at affordable cost (~105 h each); A5 part 2 (rep1_1926, p1558) not run |
| chem4c5b port (coex) | validation (a) bit-identity PASS, (b) 2a/2d/2e PASS; (c) pre-flight under A2c FAILED (reserved 25.3 / 26.9 GiB, allocated 11–12 GiB; chem4c4 has no cache release). Fix approved (Review 10); its gate was running when the work stopped |
| §3 coex groups (E3, joint, L3), A3, coex detector perturbation | **not run** |
| L2 detector perturbation (L_H, persistence, margin) | **not run** (replays queued, stopped) |

## L2 results (7 points × 16 seeds, FID (16, 256, 600), both random-copy controls, common seed 83100 for all variants)

Published L2_strict/16: p1558 16, p3709 16, engine435 16, engine5704 12, engine10653 16, rep1_1926 16, rep1_9823 16. Base vs published, Fisher exact two-sided: every point p ≥ 0.46, none flagged.

Variants:
- A0: POISSON_INV.
- A1: growth cap 16.
- A2: MULTI_START (+ slot fix + inversion sampler).
- A2a: slot-rounding fix only.
- A2b: base with random slot-cut order.
- A2c: A2 without the half-period melt floor.
- A2d: A2 with slot-index cut order.
- A2e: A2 without the floor and with slot-index order.
- A4: paired strands hydrolyse 10× more slowly.
- A5: phase-resolved melting, dt ≤ period/4.
- A6: chains ≥ 4 wash out at half rate.

### L2_strict (of 16)

| variant | p1558 | p3709 | engine435 | engine5704 | engine10653 | engine_rep1_1926 | engine_rep1_9823 |
|---|---|---|---|---|---|---|---|
| base | 15 | 15 | 14 | 9 | 16 | 16 | 16 |
| A0 | 16 | 16 | 15 | 11 | 16 | 16 | 16 |
| A1 | 15 | 16 | 15 | 11 | 16 | 16 | 16 |
| A2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A2a | 15 | 15 | 14 | 9 | 16 | 16 | 16 |
| A2b | 1 | 16 | 14 | 0 | 16 | 16 | 16 |
| A2c | 1 | 0 | 0 | 0 | 0 | 16 | 5 |
| A2d | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| A2e | 0 | 0 | 0 | 0 | 3 | 16 | 4 |
| A4 | 15 | 15 | 14 | 9 | 16 | 16 | 16 |
| A5_part1 | – | 16 | – | – | 16 | – | 16 |
| A6 | 16 | 16 | 16 | 16 | 16 | 16 | 16 |

### Copies per run (mean)

| variant | p1558 | p3709 | engine435 | engine5704 | engine10653 | engine_rep1_1926 | engine_rep1_9823 |
|---|---|---|---|---|---|---|---|
| base | 37,478 | 40,053 | 37,097 | 38,905 | 37,619 | 39,019 | 34,946 |
| A0 | 37,417 | 40,140 | 36,960 | 39,070 | 37,599 | 38,997 | 34,945 |
| A1 | 37,346 | 40,132 | 37,015 | 38,714 | 37,519 | 39,139 | 34,913 |
| A2 | 9,010 | 4,810 | 33,324 | 35,040 | 1,014 | 5,835 | 2,089 |
| A2a | 37,478 | 40,053 | 37,097 | 38,905 | 37,619 | 39,019 | 34,946 |
| A2b | 37,440 | 40,110 | 36,942 | 38,964 | 37,470 | 39,030 | 34,931 |
| A2c | 34,094 | 35,974 | 33,656 | 35,352 | 34,170 | 35,379 | 31,665 |
| A2d | 9,004 | 4,810 | 33,357 | 35,011 | 992 | 5,810 | 2,084 |
| A2e | 34,142 | 36,026 | 33,661 | 35,532 | 34,108 | 35,359 | 31,683 |
| A4 | 37,478 | 40,053 | 37,097 | 38,905 | 37,619 | 39,019 | 34,946 |
| A5_part1 | – | 40,706 | – | – | 38,970 | – | 35,738 |
| A6 | 17,779 | 19,838 | 17,396 | 19,063 | 17,597 | 19,002 | 17,188 |

### L2_share_last5 (mean)

| variant | p1558 | p3709 | engine435 | engine5704 | engine10653 | engine_rep1_1926 | engine_rep1_9823 |
|---|---|---|---|---|---|---|---|
| base | 0.423 | 0.570 | 0.556 | 0.312 | 0.673 | 0.791 | 0.707 |
| A0 | 0.451 | 0.599 | 0.528 | 0.335 | 0.663 | 0.770 | 0.696 |
| A1 | 0.448 | 0.558 | 0.511 | 0.373 | 0.654 | 0.747 | 0.675 |
| A2 | 0.012 | 0.002 | 0.151 | 0.062 | 0.000 | 0.004 | 0.000 |
| A2a | 0.423 | 0.570 | 0.556 | 0.312 | 0.673 | 0.791 | 0.707 |
| A2b | 0.175 | 0.336 | 0.239 | 0.068 | 0.349 | 0.441 | 0.299 |
| A2c | 0.164 | 0.149 | 0.002 | 0.055 | 0.177 | 0.320 | 0.217 |
| A2d | 0.047 | 0.010 | 0.159 | 0.056 | 0.000 | 0.022 | 0.001 |
| A2e | 0.160 | 0.142 | 0.002 | 0.064 | 0.196 | 0.349 | 0.194 |
| A4 | 0.423 | 0.570 | 0.556 | 0.312 | 0.673 | 0.791 | 0.707 |
| A5_part1 | – | 0.633 | – | – | 0.788 | – | 0.734 |
| A6 | 0.517 | 0.687 | 0.606 | 0.480 | 0.686 | 0.784 | 0.731 |

### Copy concentration by template lineage: Gini / top-10% share (mean)

| variant | p1558 | p3709 | engine435 | engine5704 | engine10653 | engine_rep1_1926 | engine_rep1_9823 |
|---|---|---|---|---|---|---|---|
| base_conc | 0.86 / 0.83 | 0.89 / 0.86 | 0.87 / 0.86 | 0.86 / 0.81 | 0.89 / 0.88 | 0.91 / 0.91 | 0.90 / 0.90 |
| A2b_conc | 0.59 / 0.44 | 0.72 / 0.59 | 0.57 / 0.41 | 0.47 / 0.32 | 0.60 / 0.44 | 0.64 / 0.48 | 0.63 / 0.47 |
| A2c_conc | 0.54 / 0.38 | 0.53 / 0.37 | 0.26 / 0.20 | 0.46 / 0.31 | 0.53 / 0.37 | 0.59 / 0.43 | 0.56 / 0.40 |

### A4 check: bond breaks per run in base (all strands / paired strands)

| point | breaks | paired breaks |
|---|---|---|
| p1558 | 3,210 | 0 |
| p3709 | 281 | 0 |
| engine435 | 3,209 | 0 |
| engine5704 | 1,760 | 0 |
| engine10653 | 3,004 | 0 |
| engine_rep1_1926 | 1,511 | 0 |
| engine_rep1_9823 | 1,104 | 0 |
Identity checks:
- A0_r6 ≡ A0 and A2_r6 ≡ A2 (new sampler).
- base_rec ≡ base and A2_rec ≡ A2 (Review 7 records).
- base/A2b/A2c ≡ *_conc (concentration record).
In every case, all shared stored fields are identical.

## Preliminary conclusions (L2 only; coexistence levels untested)

1. **Robust at all 7 points**:
   - A0 (sampler);
   - A1 (growth cap 4 → 16);
   - A6 (size-dependent washout: L2 everywhere, with half the copies);
   - A2a and A4: identical to base in every stored field. The slot-rounding cut never fires at these points (0 copies lost), and paired strands never break (0 paired-strand breaks in all 112 base runs).
2. **L2 depends on the one-start copy rule (A2)**: 0/16 at all 7 points. The decomposition:
   - Copying is entirely slot-limited at these points. In base, 88–99% of copy tries (p3709 88%, the others 96–99%) are cut for lack of free slots; under A2, every template-step is slot-limited.
   - **Slot-order priority (numerical artefact):** base keeps copies in slot-index order, which gives the top 10% of template lineages 81–91% of all copies. A random cut order (A2b), with copy counts unchanged, loses L2 at 2/7 points (p1558, engine5704) and halves the L2 share at the other 5. The top-10% share falls to 32–59%.
   - **Melt-wait floor:** it explains the copy collapse at the long-period points under A2 (e.g. engine10653: 1,014 vs 37,619 copies per run), but not the loss of L2. Without the floor (A2c), copy counts return to base level, yet L2 is recovered only at rep1_1926 (16) and partly at rep1_9823 (5).
   - **Slot-index priority under multi-start (A2e)** stays near A2c (3/16 at engine10653, 16 at rep1_1926, 4 at rep1_9823, 0 elsewhere). So priority order alone does not restore L2 under multi-start.
   - **Phase-resolved melting under the one-start rule (A5)** keeps L2 at the 3 long-period points tested (16/16).
3. **Reading:** the published L2 at these points arises in a slot-limited regime. There, which templates get copied is decided by the allocation rule, not by copying rates. With fair allocation (random order, or multi-start shared in proportion to rate), L2 survives at most at 2–5 of 7 points.
   - This is a dependence on numerical assumptions of the copy step. It should enter the Paper 2 assumptions table as "L2 depends on slot allocation / one-start rule".
   - Not tested: whether it also holds for E3, joint and L3 at the coexistence points.
4. **Sampler:** torch.poisson on the Intel XPU is biased and overdispersed at small means. At λ = 0.02 the pooled mean is −1.1% (z −8.6), per-seed means vary by ±5–8%, and the χ² p values are as small as 1e-287. On CPU it is clean. The published Stage 4 / P3-1 runs use it (initiation counts and similar). A0 at the L2 points changes nothing qualitatively (L2 16/16/15/11/16/16/16 vs base 15/15/14/9/16/16/16). Its effect at the coexistence points is untested.

## Files

- results/l2_summary.csv: all variants × points (means over 16 runs).
- results/l2/runs.jsonl: every per-run output.
- results/validation/: tests 1–6, R6 gate log, chem4c5b validation (a)/(b) and the failed (c) log, queue logs.
- results/diagnostics/review1_diag: the test-2 decomposition (Review 1).
- code/: all code versions used (code/versions/ holds earlier chem_p31b, chem4c4b, chem4c5b and validate_p31b versions by md5 prefix).
- docs/: implementation table with all reviews, §2 source extraction, spec, R6 gate result.

Earlier outputs of the run killed on 30 October (lot-p3/p31b/void_2026-10-30/) were voided and are not included.
