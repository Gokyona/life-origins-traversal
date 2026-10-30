# Gokyona Environment

Gokyona Environment is a local experiment platform built on the bottom-up chemistry of the life-origins-traversal
project (stage 4: monomers → polymers → templated copying → vesicles). It removes the 64-letter cap, lets 0/1
sequences evolve in an environment taken from the study's results, and shows the run live in a web page on your
own computer, including the simulated time and statistical tests for Darwinian evolution.

> **Model, not life.** Length beyond 64 letters is unvalidated. Paper result: without a sequence trait, selection
> stops after one sweep. The time scale follows the model's rate constants from the literature; it is not a
> calibrated Earth clock.

## Install

Python 3.10+. Install PyTorch for your hardware first (Intel XPU, NVIDIA CUDA, Apple MPS or CPU), then this package:

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/xpu     # Intel Arc (or cu124 / cpu; plain `pip install torch` on macOS)
pip install -e .
```

Dependencies: torch and numpy only (the web server is standard library).

## Quick start

```bash
gokyona-env            # or: python -m gokyona_env   -> opens http://localhost:9300
```

1. Optional, once: measure your hardware (writes `~/.gokyona-env/hardware.json`; the page then shows "your hardware
   allows max length N at ~X ms per step"):
   `python -m gokyona_env.engine.calibrate --sites 64 --slots 512`
2. In the page: pick a preset (e.g. **Darwin demo**), adjust sliders if you like, **Create**, then **Start**.
3. **Pause / Resume / Stop** at any time. A run checkpoints every 4 minutes and on pause/stop, and **Start** on an
   existing run resumes it.

Stop the page server with Ctrl+C. Nothing starts by itself (no autostart). One simulation at a time; on Linux with
systemd the engine runs in a memory-limited scope (`--mem-max`, default 24G). The server listens on 127.0.0.1 only.
All data live in `~/.gokyona-env` (or `--data-dir`).

UI details (options, panels, mock data for testing): see [README_UI.md](README_UI.md).

## What is simulated

The engine (`gokyona_env/engine/gk_chem.py`) is the stage-4 chemistry (chem4c4 v2 with the 4c-5 lineage ideas):
activated monomers from feedstock, polymerisation and hydrolysis, wet–dry cycling, catalysed ligation, templated
copying with errors (copies are the complementary strand and start at letter 0), duplexes, amphiphiles, vesicles that
grow, divide and wash out. Sequences are stored in W 64-bit words per chain, so the maximum length is 64·W.

**Validation:** with W = 1 the engine is **bit-identical** to the study's chem4c4 (every output and the final state,
4 of 4 test cases with copying and vesicles; `tests/identity_w1.py`). Longer chains (W > 1) use the same rules;
no bit is ever set beyond a chain's length (`tests/smoke_w.py`).

**Presets:** engine12538 and the 9c joint points of stage 4c-4 (heredity + vesicles), where the 5a point study found
coupled protocells (L3 strict). **Darwin demo** = engine12538 with the sequence trait on.

## The sequence trait (optional; off = the study's physics)

- **Landscape** (fixed seed): *functional domain* (default) — the first 24 letters, mean of a 6-mer table over its
  19 windows, z-scored against random sequences (so the strength `s` is in SD units); *epistatic* — the window value
  also depends on its position (rugged landscape); *additive mean* — the whole chain (reference; nearly flat).
- **Pair score** (default on): a chain's trait = mean of its own and its complement's domain, so the trait is
  inherited every generation (copies are the complementary strand).
- **Modes** (what the trait does):
  - *retention* — a vesicle's washout rate is multiplied by exp(−s·z̄) of its enclosed chains (**assumption**);
  - *division* — a vesicle's division rate × exp(s·z̄) (**assumption**);
  - *membrane* — a vesicle's membrane uptake × exp(s·z̄);
  - *stability* — a chain's bond breaking × exp(−s·z);
  - *self* / *shared replicase* — copy rate × exp(s·z) / × (1 + s·max(0, mean replicase z)).
- **Multi-start copying** (option): Poisson(λ·dt) copy starts per template and step instead of at most one.
- **Built-in control:** half of the sites use a landscape re-shuffled every step (same score distribution, no
  heritable link), so the Darwin panel compares real against control in the same run.

## Experimental options (version 1.0; all default off)

All of these are switched off unless a run's config sets them. With all of them off, the engine at one 64-bit word per
chain (W = 1) is bit-identical to chem4c4, the Stage-4 engine of paper 2 (`tests/identity_w1.py`); every option passes
`tests/options_invariants.py`. They are research options and are not validated beyond the smoke tests and screens
described in "What we found".

| Config key | What it does |
|---|---|
| `ration` | monomer and slot rationing among started copies: `random` (default for new runs from the page), `vesicle`, or `slot` (legacy chem4c4 order; configs without the key) |
| `census_years` | founder-tag census interval for clades by descent (default 100) |
| `trait.env_landscape` | environment-dependent landscape: two landscapes mixed by w(t); schedules `constant`, `season`, `trend`, `catastrophe`; `variable` T / pH / salinity; `gradient`, `band_amp`, `band_period_years`, `season_amp` (latitude) |
| `env_chem` / `trait.env_landscape.couple_chemistry` | the environment variable acts on the chemistry (phase table, membrane window, CVC) |
| `catastrophe`, `impact` | per-site and global vesicle-killing events; impact env and feed pulses |
| `migration_m`, `migration_mode` (`random`, `ring`, `line`), `half_pool_mix`, `local_pool_mix` | vesicle migration and pool mixing between sites |
| `uv_gradient` | per-site UV on the latitude axis (synthesis, activation, UV breaks) |
| `protected_domain`, `promoter_required`, `domain_break_guard` | break fragments without the start mark are unscored / uncopyable; breaks inside the domain suppressed |
| `k_max`, `no_reencapsulation`, `division_mode` (`content`), `division_n0`, `segregate`, `division_membrane_req` | single-genotype vesicles: cap per vesicle, no re-entry, content-triggered division, strand segregation |
| `division_trait_s`, `access_trait_s` | couple division threshold / monomer access to the vesicle trait (imposed couplings, for comparison) |
| `cat_motif`, `local_conc`, `salt_per_step` | motif-based heritable catalysis, mass action with internal vesicle concentrations, per-step re-hash control |
| `vesicle_min_size` | vesicles below a minimum membrane break down and release their contents |
| `growth_uncapped`, `copy_kinetic`, `monomer_competition` | time-step-convergent growth and copying (no 4-letter growth cap; copies draw monomers over time; rate-proportional monomer split) |
| `copy_physical`, `growth_physical`, `long_primers`, `mineral_oligo` | physical (literature-scale) elongation and growth rates, chance long primers, mineral-surface oligomer supply |
| `vesicle_stats` | per-vesicle statistics in the history (no effect on dynamics) |

## Recommended Darwin demo (tested 29 Oct 2026 on an Intel Arc Pro B70)

Preset **Darwin demo** (engine12538 environment with feedstock ×4 — above the study's feed range — and the trait
on: retention, s = 4, functional domain, epistatic, pair score), size 64 sites × 512 slots, **64 vesicle slots per
site**, L_MAX 1024. Measured: 111 ms per step, ~1,400 simulated years per wall-clock hour, 0.76 GiB GPU memory,
~4 GiB host memory (flat after warm-up), checkpoint 27 MB. In short tests the real half's mean trait rose
~0.5–1 SD above the control within the first ~10–20 simulated years ("selection on the trait"); a second shift
("sustained evolution") was not seen within 116 simulated years.

## Darwin panel: definitions (fixed before the platform's runs)

- **Selection on the trait** (one shift): real − control mean z > 3 × the SD of the control's report-to-report
  values (last 20 reports, ≥ 5) for ≥ 5 consecutive reports; the plateau = the mean difference over those reports.
- **Sustained evolution**: ≥ 2 successive shifts (each a new plateau passing the same rule above the previous one),
  or ≥ 2 successive lineage sweeps (a new dominant 24-letter lineage with ≥ 30% share whose score gain exceeds the
  99th percentile of the control half's jumps).
- **ρ(w, z)**: Spearman correlation between a vesicle's trait z and its number of descendants one report later,
  against a 199-permutation null (99th percentile), per half.
- **L4 lamp**: ρ(w, z) above its null99 in ≥ 3 of the last 5 reports (a single report beyond the 1% null is
  expected by chance about once in 100 reports), together with a selection status other than "no selection".
- The ladder lamps (E1 … L4) are live proxies of the study's batch detectors (L3 lamps are not computed live).

### Clades by descent (founder tags; added 29 Oct 2026)

24-letter lineage keys change with every copy error, so a clade that spreads while mutating never shows up as one
"lineage". Founder tags follow descent instead:

- **Census:** at the first step and then every `census_years` simulated years (default 100; config key
  `census_years`), every chain gets a new unique tag (census number in the high 32 bits). Copies and fragments
  inherit their parent's tag, migration moves it with the chain, new chains from monomers get none. Tags are saved in
  the checkpoint. With the tags the engine is otherwise unchanged (W = 1 identity still holds).
- **Clade** = all chains carrying one tag; per half the page shows the top clade's share of the half's tagged lineage
  chains, its sites, its mean z and the number of clades.
- **Clade event** (logged for both halves): a tag first reaches **≥ 30 %** of its half's tagged chains.
- **Counted clade sweep** (real half only), all three must hold:
  1. share ≥ 30 %;
  2. share > the **95th percentile of the control half's top-clade share** over all reports so far (needs ≥ 20
     control reports; before that the event is "pending" and is judged, with its share at the event, once there are
     20) — small isolated sites fix neutral clades by drift, and the control half measures how often;
  3. the clade's mean z (length-weighted, at the event) > the **real half's mean z at its census** (first report after
     that census) — a fitter clade, not just drift.
- **Sustained evolution** additionally = **≥ 2 counted clade sweeps**. The rule was fixed before any run used it.

## What we found (29 Oct 2026; Intel Arc Pro B70; md5 = first 8 hex digits of each result file)

All runs use the Darwin demo (engine12538 environment, feed x4, retention s = 4, functional domain, epistatic,
pair score). Criteria for the later tests were written to file before the runs (md5 given).

(The result files and run data named below are not part of this package.)

1. **Which trait can be selected at all** (rounds 1-4, short GPU runs; `test_results_2026-10-29/runs/s2_*`, `s3_*`,
   `s4*`, `scan_*`; model diagnostics `saturation_loss.json` d6b95803): copy starts are saturated at the paper's time
   step (lambda dt >= 10^6), monomers bind (78-99% of started copies dropped), 93-99.9% of lineage chains are lost with
   washed-out vesicles, divisions are triggered from outside. Only a trait on vesicle **retention** gave a response.
2. **Long run, slot-order rationing** (`darwin_long1`, 64 x 512, 2009 simulated years; history 7faf639d, status
   4e10af5b): one fast shift (real - control +0.41 at 5 years; plateau 0.569), then flat for 2000 years
   (250-year blocks 0.47-0.52); no lineage sweep (top lineage < 1.2%), L4 never persistent.
3. **Migration** (`darwin_mig1`, m = 1e-3, 2381 years; history 85dd1cb4, status f2d2edfa; plus a scan of m = 1e-4 ...
   1e-2, `runs/mig_*`): same picture - one shift (plateau 0.477), flat, no sweep, top lineage <= 1.2%; no hidden
   cloud of near-identical domains (largest Hamming-3 ball 1.0%; domain_clusters_2381y.txt 4746e082).
4. **Why it stops after one step** (CPU analyses on darwin_long1; analysis_first_shift.json f2e64b6b,
   fresh_domains.json 5f11f4b7, copy_bias.json 9aeecd0a): 79% of the variation in z lies within vesicles (invisible
   to vesicle washout); the landscape is far from exhausted (max z 3.48, best chain 2.48, mean 0.37); copy errors are
   negligible. The plateau is a balance: washout selection (+0.003 per step, chain level) against (a) copying biased
   to low-z templates (-0.002) and (b) bond-break fragments that start with an unselected, random domain (-0.001).
   The copy bias comes from **slot-order rationing**: almost all copies went to the first 20% of slots, and long
   (fitter) templates lost the monomer race.
5. **Fair rationing** (2 x 2 screen `fair_results.json` efc2d09b, criteria c46a1b8e; two 2-h runs at full size,
   criteria faa20215, results `fair_long_s1_results.json` b881022f and `fair_long_s2_results.json` 7b6d0a87):
   random-order rationing removes the copy bias (copy-template z = population z) and raises the stable plateau from
   ~0.50 to ~0.78 SD within ~250 years, replicated in two seeds (criteria a and d pass in both). After that the level
   is flat for ~2300 years (late slope +-0.0002 per 100 years); the rule's second shift (+0.11) lies inside the
   plateau band; no lineage or clade takeover (< 2%), L4 never persistent.
6. **Fragment tests** (screen3, criteria db332d7b, results `screen3_results.json` d353362a): suppressing breaks in
   the first 24 letters (guard; 5.5% of breaks) has no effect, as predicted; making break fragments unscored and
   uncopyable (P1) raises the z of marked chains (0.87) but not the population z over all chains (-0.02 / -0.09 vs
   fair rationing alone), and cuts copying of scored templates by ~83%. Negative.

7. **Physics-based heredity tests (30 Oct 2026):** at literature-scale copy rates (option `copy_physical`, ~1–2 letters per
   hour instead of ~1e9–1e11) no copying happens from short primers; templated ligation and chance long primers cannot start
   heredity from scratch either, also with stable pairing (T_m8 50 C), physical growth rates and a mineral oligomer supply.
   The platform's heredity therefore rests on the fast copying of the Stage-4 engine.

**In short:** in this chemistry, selection on a sequence trait produces one fast step to a stable plateau. Part of
the original plateau was an artefact of slot-order rationing (now fixed by default); the rest is set by fragment
regeneration and copy capacity. No option tested so far gives open-ended adaptation, sweeps or clades.
Run data of the long runs live in `~/.gokyona-env/runs/` on the machine that ran them; they are not part of the package.

## Files

| Path | What |
|---|---|
| `gokyona_env/server.py`, `gokyona_env/web/` | local web server (127.0.0.1:9300) and page |
| `gokyona_env/engine/gk_chem.py` | the chemistry (multi-word sequences; W = 1 equals chem4c4) |
| `gokyona_env/engine/gk_cat3.py`, `gk_cat4.py`, `gk_raf.py`, `words.py` | catalysis, RAF hashing, word arithmetic |
| `gokyona_env/engine/trait.py` | sequence → trait rule and the control landscape |
| `gokyona_env/engine/darwin.py` | live statistics, sweeps, lineage tree, ρ(w, z), trait-response rule, clade sweeps |
| `gokyona_env/engine/run.py` | one open-ended run: status, history, checkpoints, pause/resume/stop |
| `gokyona_env/engine/calibrate.py` | memory and step time against length; max length for your hardware |
| `gokyona_env/engine/devices.py`, `thermal.py` | device choice (xpu → cuda → mps → cpu), temperatures; optional GPU guard |
| `gokyona_env/data/presets.json`, `params.json` | presets and slider definitions |
| `tests/` | identity, smoke and invariant tests; `tests/research/` analysis scripts of the 29–30 Oct 2026 screens (not installed with the package; `identity_w1.py` needs `code/lot-s4/` of the paper repository) |

License: MIT (see LICENSE).
