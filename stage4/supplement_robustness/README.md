# Stage 4 supplementary study: how much does the heredity result depend on built-in rules?

**31 October 2026, Toronto time.** This study was planned as step P3-1b of Part III, after an external review of the model's assumptions. That is why file and code names carry "p31b".

The part that was run tests the Stage 4 results, not Part III. It is therefore published as a supplement to Stage 4 and [Paper 2](../../paper/COMPASS_stage4.md), which was revised with its findings on 31 October 2026 (abstract, Sections 6–8).

The owner stopped the work early, at 05:25, after the heredity (L2) part. Start with **[SUMMARY.md](SUMMARY.md)**, which has every table.

## Question

The heredity level L2 was confirmed at 7 environments in Stage 4. Does it still appear there if built-in simplifications of the model are changed one at a time?

The 7 points were each re-run with 16 seeds and both random-copy controls, under these changes:
- an exact Poisson sampler;
- a higher cap on untemplated growth;
- several copy starts per template and step;
- the order in which copies get free slots;
- protection of paired strands;
- phase-resolved melting;
- size-dependent washout.

## Findings (L2 only)

1. **Robust.** L2 does not depend on:
   - the Poisson sampler;
   - the untemplated-growth cap (4 → 16);
   - size-dependent washout;
   - phase-resolved melting (3 long-period points).

   Paired-strand protection and the slot-rounding fix change nothing, because neither event occurs at these points.
2. **Depends on the copy step.** Copying is slot-limited: 88–99% of copy attempts find no free slot. The engine fills free slots in slot-index order, so the top 10% of template lineages get 81–91% of all copies.
   - **Random slot order, same copy count:** L2 is lost at 2 of 7 points and its share halves at the other 5.
   - **Several copy starts per step, shared by rate:** L2 survives at no point with the melt wait, and fully at only 1 point without it.

   So the published L2 arises largely from a priority rule of the numerical scheme, not from copying rates.
3. **XPU Poisson bias.** `torch.poisson` on the Intel XPU is biased at small means: −1.1% at λ = 0.02. An exact sampler leaves L2 unchanged at these points.
4. **Not done:**
   - the coexistence groups (E3, joint, L3) and A3;
   - detector-threshold perturbation;
   - A5 at the remaining points;
   - the salt correction, which was decided but not implemented;
   - the 4c-5 port's memory pre-flight, which failed on cached memory. Its fix was approved but not run.

**Meaning for Stage 4.** Together with the two earlier limits in Paper 2 (copying 10⁹–10¹¹ times faster than measured chemistry, and no real base pair matching the heredity constants), the heredity results of Stage 4 should be read as properties of the model's numerical copying rule, not as evidence that heredity can start from scratch in real chemistry.

## Files

- **Plan:** [docs/part3_p31b_spec.md](../../docs/part3_p31b_spec.md) (copy in `docs/`).
- **Implementation choices and reviewer decisions:** [docs/IMPLEMENTATION_P31B.md](docs/IMPLEMENTATION_P31B.md), Reviews 1–10. Each was written before the results it governs.
- **Final code:**
  - `code/chem_p31b.py` 44f30af0;
  - `code/chem4c4b.py` 61d9a4f4;
  - `code/chem4c5b.py` 6abf1b2b.

  Earlier versions are in `code/versions/`. The Stage 4 modules chem4c4 and chem4c5 are published with Paper 2 and are not repeated here.
- **Results:**
  - `results/l2_summary.csv`, with every variant × point; `code/make_summary.py` rebuilds it from `results/l2/runs.jsonl`;
  - the validation outputs and logs.
- **Integrity:**
  - `MD5SUMS.txt` lists all files (`md5sum -c MD5SUMS.txt`).
  - The archive as transferred from the B360 was p31b_publish.tar.gz: 2,021,182 bytes, md5 abc7c5e300539d9c27c34a4c118b91ad.
  - This README replaces the archive's short README.
