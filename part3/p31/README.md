# Part III, P3-1: implementation and validation files

**30 October 2026, Toronto time.** Validation only: no map and no long run has been started. Maps begin after the owner reviews these results.

- **Report:** [REPORT_P31.md](REPORT_P31.md), with every test, the cold-acceptance table, implementation choices (a)–(n), cost and pre-flight, temperatures and science flags.
- **Specification:** [docs/part3_p31_spec.md](../../docs/part3_p31_spec.md). **Cold-benefit analysis:** [docs/part3_cold_gap.md](../../docs/part3_cold_gap.md) (copies in `ref/`, with [cold_gap_report.md](cold_gap_report.md) and [inhibit_report.md](inhibit_report.md)).
- **Final code:** `code/chem_p31.py`, md5 cfc1e3db797789fadda313d038d0c981. Earlier versions are kept for traceability. `code/chem4c3.py` is an unchanged reference copy of the Stage 4c-3 chemistry.
- **Integrity:** `MD5SUMS.txt` lists all files (`md5sum -c MD5SUMS.txt`). Archive as transferred from the B360: p31_publish.tar.gz, 389,719 bytes, md5 ae5098d8fdac065ff3808524a2b98670.
- **Literature:** `literature/` holds text extracts and MD5s of Aitken et al. (2023) main text and SI, used for the cold-acceptance test; no PDFs.

## Summary

All nine validation tests pass on the final code:

1. Duplex thermodynamics.
2. Extension kinetics.
3. Pulse integration.
4. Bit-identical Stage 4 regression.
5. Activation calibration.
6. Cold acceptance against Aitken et al. (2023).
7. Error threshold near ε·L_copy ≈ 0.85.
8. 0/64 null false positives.
9. Bit-identical kill and resume.

Cold acceptance fits the 5-day and pre-activated results with one parameter set, but it has two caveats. The 24-hour single-letter result keeps a residual of ×2.25 per letter, recorded as unexplained. The accepted parameter region also lies at the edge of the approved ranges.

The memory pre-flight passes at batch 512 with per-step cache release: 11.5 GiB VRAM, 41% of the cap.

**Science flag:** with measured duplex thermodynamics, primers of 2–3 letters rarely extend. In 640 random A/s2U environments (a diagnostic, not a map), 22 made any copy and none made a long copy.

`results/regress_4c3_stored.pkl` is a Python pickle of stored Stage 4 results used by test 4; load it only if you trust this repository.
