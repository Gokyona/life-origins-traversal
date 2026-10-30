#!/usr/bin/env python3
"""Fake Gokyona run for testing the UI without the engine (stdlib only, no torch, no GPU).

Writes <data-dir>/runs/<name>/{config.json,status.json,history.jsonl} and updates them every
second with plausible numbers (incl. the v1.1 darwin fields, sweeps and a lineage tree).
Obeys control.json (pause / resume / stop) like the real engine.

    python3 tools/mock_run.py [--data-dir ~/.gokyona-env] [--name mock] [--seconds 0] [--interval 1]
    (--seconds 0 = run until Ctrl+C or a "stop" command)
"""
import argparse
import datetime
import hashlib
import json
import math
import os
import pathlib
import random
import time

BANNER = ("MOCK DATA - not a simulation. Model, not life. Length beyond 64 is unvalidated. "
          "Paper result: without a sequence trait, selection stops after one sweep.")
TIME_NOTE = "Time scale follows the model's rate constants from the literature; it is not a calibrated Earth clock."


def write_atomic(path, obj):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1))
    os.replace(tmp, path)


def key(*parts):
    return hashlib.md5("/".join(map(str, parts)).encode()).hexdigest()[:6]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(pathlib.Path.home() / ".gokyona-env"))
    ap.add_argument("--name", default="mock")
    ap.add_argument("--seconds", type=float, default=0.0)
    ap.add_argument("--interval", type=float, default=1.0)
    a = ap.parse_args()

    rng = random.Random(7)
    d = pathlib.Path(a.data_dir).expanduser() / "runs" / a.name
    d.mkdir(parents=True, exist_ok=True)
    for f in ("history.jsonl", "control.json"):
        (d / f).unlink(missing_ok=True)
    trait = {"enabled": True, "mode": "shared", "s": 1.0, "k": 6, "seed": 2026, "L_rep": 24,
             "fidelity_link": False, "s_f": 0.0}
    cfg = {"name": a.name, "preset": "darwin_demo", "env": {"T": 4.47, "pH": 7.58}, "trait": trait,
           "lmax": 256, "sites": 64, "slots": 512, "vram_budget_gb": 20, "report_every_s": 1,
           "checkpoint_every_s": 240, "created": datetime.datetime.now().replace(microsecond=0).isoformat()}
    write_atomic(d / "config.json", cfg)

    state, wall, sim_h, steps = "running", 0.0, 0.0, 0
    sweeps, nodes = [], []
    dom = {"real": None, "control": None}
    dom_score = {"real": 0.05, "control": 0.05}
    t0 = time.time()
    last_ctl = None
    with open(d / "history.jsonl", "a") as hist:
        while True:
            # ---- control
            try:
                ctl = json.loads((d / "control.json").read_text())
            except (OSError, ValueError):
                ctl = None
            if ctl and ctl != last_ctl:
                last_ctl = ctl
                state = {"pause": "paused", "resume": "running", "stop": "stopped"}.get(ctl.get("cmd"), state)
            if a.seconds and time.time() - t0 >= a.seconds:
                state = "finished"

            if state == "running":
                wall += a.interval
                sim_h += 12.0 + rng.random() * 4
                steps += 140
                # lineage events
                for half in ("real", "control"):
                    p = 0.10 if half == "real" else 0.05
                    if dom[half] is None or rng.random() < p:
                        parent = dom[half]
                        gain = rng.uniform(0.02, 0.12) if half == "real" else rng.uniform(-0.05, 0.05)
                        dom_score[half] = max(0.01, min(0.99, dom_score[half] + gain))
                        k = key(half, sim_h, rng.random())
                        nodes.append({"id": k, "parent": parent, "born_h": round(sim_h, 1),
                                      "score": round(dom_score[half], 3), "share": 0.05,
                                      "len": rng.randint(24, 96), "half": half})
                        dom[half] = k
                    # shares grow for the current dominant, decay for the others
                    for n in nodes:
                        if n["half"] != half:
                            continue
                        n["share"] = round(min(0.95, n["share"] * 1.25 + 0.02) if n["id"] == dom[half]
                                           else n["share"] * 0.85, 3)
                        if (n["id"] == dom[half] and n["share"] > 0.5
                                and not any(s["key"] == n["id"] for s in sweeps)):
                            sweeps.append({"sim_h": round(sim_h, 1), "half": half, "key": n["id"],
                                           "score": n["score"], "share": n["share"],
                                           "beyond_null": half == "real" or rng.random() < 0.1})
                nodes = nodes[-80:]

            tt = wall
            mean_len = 4 + 20 * (1 - math.exp(-tt / 30)) + rng.gauss(0, 0.3)
            max_len = min(256, int(16 + 6 * tt + rng.randint(0, 6)))
            distinct = int(50 + 25 * tt + rng.randint(0, 20))
            ves = 7.6 + 0.5 * math.sin(tt / 5) + rng.gauss(0, 0.05)
            gen_max = int(3 * tt)
            sm_r = 0.02 + 0.5 * dom_score["real"] * (1 - math.exp(-tt / 20))
            sm_c = 0.02 + 0.05 * rng.random()
            sd_r, sd_c = dom_score["real"], dom_score["control"] * 0.3
            rho_r = min(0.9, 0.05 + 0.02 * tt + rng.gauss(0, 0.02))
            rho_c = rng.gauss(0, 0.05)
            n99_r, n99_c = 0.22 + rng.gauss(0, 0.005), 0.23 + rng.gauss(0, 0.005)
            real_sw = [s for s in sweeps if s["half"] == "real" and s["beyond_null"]]
            dstatus = ("sustained evolution" if len(real_sw) >= 2 else "one sweep" if real_sw else "no selection")
            edges = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
            counts = [int(max(0, 9000 * math.exp(-i * 0.9) * (1 + tt / 40) + rng.randint(0, 30))) for i in range(len(edges) - 1)]
            top = []
            for i in range(8):
                L = rng.choice([32, 48, 64, 96, 128, 160])
                top.append({"len": L, "count": int(200 / (i + 1)) + rng.randint(0, 5),
                            "score": round(max(0, dom_score["real"] - 0.05 * i), 3),
                            "seq": "".join(rng.choice("01") for _ in range(L)), "key": key("top", i, int(tt / 5))})
            status = {
                "name": a.name, "state": state, "error": None, "preset": "darwin_demo", "env": cfg["env"], "trait": trait,
                "lmax": 256, "lmax_auto_max": 8192, "W": 4, "sites": 64, "slots": 512,
                "time": {"sim_h": round(sim_h, 2), "sim_days": round(sim_h / 24, 3), "sim_years": round(sim_h / 8766, 5),
                         "cycles": round(sim_h / 63.4, 2), "period_h": 63.4, "wall_s": wall,
                         "speed_sim_h_per_wall_s": 14.0, "eta_wall_h_to_1y": round(max(0.0, (8766 - sim_h) / 14 / 3600), 3)},
                "perf": {"steps": steps, "step_ms": round(85 + rng.gauss(0, 4), 1), "gpu_temp_c": int(60 + tt / 10 + rng.randint(0, 2)),
                         "vram_gib": round(7.3 + rng.random() * 0.2, 2), "vram_budget_gib": 20,
                         "device": "mock (no GPU)", "hw_note": "mock: your hardware allows max length 4096 at ~12 steps/s"},
                "chains": {"n": 21000 + int(30 * tt), "mean_len": round(mean_len, 2), "max_len": max_len,
                           "len_hist": {"edges": edges, "counts": counts}, "nuc_share": 0.7, "enclosed_share": 0.95},
                "lineages": {"distinct_ge_LH": distinct, "top": top},
                "vesicles": {"per_site": round(ves, 2), "gen_max": gen_max, "gen_mean": round(gen_max * 0.52, 1),
                             "divisions": int(400 * tt), "molecular_generations": round(tt * 7.0, 1)},
                "levels": {"E1": True, "E2": True, "L1": tt > 3, "L2": tt > 6, "E3": True, "L3_loose": tt > 8,
                           "L3_strict": tt > 12, "L4": None if tt < 4 else len(real_sw) >= 2},
                "trait_stats": {"score_mean": round(sm_r, 4), "score_max": round(min(1, sd_r + 0.1), 3),
                                "score_dominant": round(sd_r, 3), "rho_wz": round(rho_r, 3), "null99": round(n99_r, 3)},
                "darwin": {
                    "status": dstatus,
                    "status_note": "sustained = >= 2 successive sweeps by a fitter lineage, each beyond the neutral null",
                    "real": {"score_mean": round(sm_r, 4), "score_max": round(min(1, sd_r + 0.1), 3), "score_dominant": round(sd_r, 3),
                             "rho_wz": round(rho_r, 3), "null99": round(n99_r, 3), "L4": len(real_sw) >= 2},
                    "control": {"score_mean": round(sm_c, 4), "score_max": round(min(1, sd_c + 0.3), 3), "score_dominant": round(sd_c, 3),
                                "rho_wz": round(rho_c, 3), "null99": round(n99_c, 3), "L4": False},
                    "sweeps": sweeps,
                    "tree": {"nodes": nodes},
                },
                "banner": BANNER, "time_note": TIME_NOTE,
            }
            write_atomic(d / "status.json", status)
            if state == "running":
                hist.write(json.dumps({
                    "wall_s": wall, "sim_h": round(sim_h, 2), "mean_len": round(mean_len, 2), "max_len": max_len,
                    "distinct": distinct, "vesicles": round(ves, 2), "gen_max": gen_max,
                    "score_mean": round(sm_r, 4), "score_max": round(min(1, sd_r + 0.1), 3), "score_dominant": round(sd_r, 3),
                    "rho_wz": round(rho_r, 3), "n_chains": 21000 + int(30 * tt),
                    "score_mean_real": round(sm_r, 4), "score_mean_ctl": round(sm_c, 4),
                    "score_dom_real": round(sd_r, 3), "score_dom_ctl": round(sd_c, 3),
                    "rho_wz_real": round(rho_r, 3), "rho_wz_ctl": round(rho_c, 3),
                    "null99_real": round(n99_r, 3), "null99_ctl": round(n99_c, 3)}) + "\n")
                hist.flush()
            if state in ("stopped", "finished"):
                print(f"mock run {a.name}: {state}", flush=True)
                return
            time.sleep(a.interval)


if __name__ == "__main__":
    main()
