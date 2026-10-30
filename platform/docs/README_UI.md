# Gokyona Environment: local web UI

A visual front end for the Gokyona origin-of-life simulation. **It is a model, not life.** The UI
always shows the engine's honesty banner and time note at the top.

## Install

You need Python 3.10 or newer. Install PyTorch for your hardware first, then this package:

```bash
python3 -m venv .venv && . .venv/bin/activate

# pick ONE torch build for your hardware (see https://pytorch.org/get-started/locally/):
pip install torch --index-url https://download.pytorch.org/whl/xpu     # Intel Arc / XPU
pip install torch --index-url https://download.pytorch.org/whl/cu124   # NVIDIA CUDA (match your driver)
pip install torch                                                      # macOS (Apple MPS), or default
pip install torch --index-url https://download.pytorch.org/whl/cpu     # CPU only

pip install -e .            # from the gokyona-env folder (installs numpy, and torch if still missing)
```

The web server itself uses only the Python standard library. torch is needed only by the engine.

## Quick start

```bash
gokyona-env                      # or: python -m gokyona_env
```

This opens http://localhost:9300 in your browser (add `--no-browser` to skip that). Options:

| option | default | meaning |
|---|---|---|
| `--data-dir` | `~/.gokyona-env` | runs live in `<data-dir>/runs/<name>/`; `hardware.json` is read from here |
| `--host` | `127.0.0.1` | bind address. Keep this on loopback: the API can start jobs |
| `--port` | `9300` | port |
| `--mem-max` | `24G` | memory limit for the systemd scope (used only where `systemd-run --user` exists) |
| `--no-browser` | off | do not open a browser |

Without installing, from this folder: `python3 -m gokyona_env` (any Python 3.10 or newer).

**Nothing starts by itself.** The engine runs only when you press **Start**. Starting is refused while
another job is active: any `gk-*` systemd scope, any `lot-*` unit other than `lot-tempwatch` (checked
only where systemctl exists), or any run whose `status.json` says running, starting or paused. The
"Start check" line under the buttons shows the reason. If an engine crashed and left a stale "running"
status, edit or delete that run's `status.json`.

### Testing without the engine

```bash
python3 tools/mock_run.py --data-dir ~/.gokyona-env --seconds 120
```

This writes a fake run `mock` with plausible numbers every second, including Darwin data. It obeys
Pause, Resume and Stop. Delete `<data-dir>/runs/mock` afterwards.

## Using it

1. **Run control.** Pick a preset; the environment sliders take the preset's values. Sliders marked
   ㏒ use a log scale, and the physical value and unit appear on the right. Moved sliders turn amber;
   only those are sent, and every other variable keeps its exact preset value. Set the sequence trait
   (enabled, mode self/shared, s 0–2, L_rep, k = 6 read-only, seed, fidelity link with s_f), plus
   L_MAX, sites and slots, and a run name (letters, digits, `_`, `-`, at most 40 characters).
   **Create** writes `config.json` but does not start the run. **Start** asks for confirmation and
   launches the engine. **Pause**, **Resume** and **Stop** write `control.json`; the engine acts on
   it at its next report.
2. **Run selector** (top bar): choose which run to watch. The pill shows its state.

## Panels

| panel | shows |
|---|---|
| Header | honesty banner, time note, hardware line (`perf.device`, `perf.hw_note`, or `hardware.json`) |
| Time | sim hours, days and years, progress to 1 sim year, cycles, cycle period, vesicle and molecular generations, speed (sim h per wall s), wall time, ETA to 1 sim year |
| Chain length | chain count, mean and max length, L_MAX, nucleotide and enclosed shares; mean and max over time (log y); current length histogram (log count) |
| Lineages | distinct sequences (≥ L_H) over time; top sequences as 0/1 colour strips (blue 0, amber 1, first 128 letters) with length, count, score and key |
| Vesicles | vesicles per site, generation max and mean, divisions, molecular generations; chart over time |
| Ladder levels | lamps E1, E2, L1, L2, E3, L3 loose, L3 strict, L4 (green reached, dim not reached, grey unknown) |
| Trait / selection | trait settings; score mean, max and dominant over time; ρ(w,z) over time against the current null99 |
| Machine | device, steps, step time, GPU temperature, VRAM against budget (bar); step time over wall time (from history, else sampled by the browser) |
| Darwin | big status light (red no selection, amber one sweep, green sustained evolution) with its note; real and control half cards (scores, ρ vs null99, L4); real vs control score chart (solid real, dashed control) with sweep triangles on the time axis (filled = beyond null); ρ(w,z) vs null99 for both halves; lineage tree (x = birth time, colour = score, size = share; real and control in separate bands) |

The UI polls status every 2 s and fetches history incrementally (`since=N`). A missing or null field
shows "—".

## API (for reference; see CONTRACT.md)

`GET /api/runs`, `/api/status?run=`, `/api/history?run=&since=`, `/api/presets`, `/api/params`,
`/api/check[?run=]`, `/api/hardware`; `POST /api/create`, `/api/start`, `/api/control`.
