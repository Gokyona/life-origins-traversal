"""Memory pre-flight for P3-1 (copy of lot-s4/preflight_mem.py, md5 2e70b907, adapted to chem_p31 and a
pair chemistry; standing rule: peak <= 75% of the caps, no growth). Point sets: "random" (Sobol, full P3-1
space), "dense" (worst case for templates and run length: cold base 0-5 C, high feed, concentration factor and
activating agent, long available time, slow hydrolysis, many short pulses resolved or not) and "extreme"
(every dimension at 0.02 or 0.98). The 4c-3 description follows.
Memory pre-flight for stage-4 GPU jobs (standing rule from the owner, 27 Oct 2026: every new
or changed program passes this before any long run).

Runs map-sized batches of chem4c3 (batch size and fidelity as the job uses them) with both paired
controls on worst-case point sets, and records per batch: peak XPU memory allocated and reserved
(torch.xpu.max_memory_*), process RSS after the batch, and the scope's cgroup memory.current /
memory.peak. Point sets: "random" (Sobol, full space), "box2" (the scan region with the most
templates per site), "extreme" (every dimension at 0.02 or 0.98, at random). 
For each batch size in --sizes: 4 batches (random, box2, extreme, random again). Pass for a size:
peak XPU reserved and host peak (cgroup memory.peak) <= --limit (default 0.80) of the caps
(XPU: chem4c3.XPU_FRACTION x card memory; host: the scope's MemoryMax), and memory.current
after the second random batch within 10% of the first (no growth). The chosen batch size is the
passing one with the lowest s/run on the random batches.
    python preflight_mem.py --out preflight_map4c3.json --sizes 512,768,1024,1536 [--fid screen]
"""
import argparse, json, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_p31 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402

dev = torch.device("xpu")


def cgroup_file(name):
    rel = open("/proc/self/cgroup").read().strip().split("::")[-1]
    p = f"/sys/fs/cgroup{rel}/{name}"
    try:
        v = open(p).read().strip()
        return None if v == "max" else int(v)
    except OSError:
        return None


def rss():
    for ln in open("/proc/self/status"):
        if ln.startswith("VmRSS:"):
            return int(ln.split()[1]) * 1024
    return 0


def dense(n, seed):
    """Worst case for templates, copies and run length (u-space intervals)."""
    X = qmc.Sobol(d=chem.DIM_P31, scramble=True, seed=777 + seed).random(n).astype(np.float32)
    ix = {k: chem.NAMES_P31.index(k) for k in chem.NAMES_P31}
    box = {"T": (0.0, 0.1), "feed": (0.8, 1.0), "x_nuc": (0.8, 1.0), "thermo": (0.6, 1.0), "act_a": (0.8, 1.0),
           "avail_time": (0.8, 1.0), "pH": (0.2, 0.5), "tau": (0.8, 1.0), "UV": (0.0, 0.2), "k_UV": (0.0, 0.2),
           "T_p": (0.7, 1.0), "tau_p": (0.0, 0.5), "f_p": (0.3, 1.0), "k_c": (0.6, 1.0), "L_t": (0.0, 0.3),
           "k_act": (0.5, 1.0), "Mg": (0.7, 1.0), "evap": (0.5, 1.0), "depth": (0.0, 1.0), "period": (0.0, 0.5)}
    for k, (lo, hi) in box.items():
        X[:, ix[k]] = lo + (hi - lo) * X[:, ix[k]]
    return X


def point_set(name, n, k):
    if name == "random":
        return qmc.Sobol(d=chem.DIM_P31, scramble=True, seed=99 + k).random(n).astype(np.float32)
    if name == "dense":
        return dense(n, k)
    rng = np.random.default_rng(5 + k)
    return np.where(rng.random((n, chem.DIM_P31)) < 0.5, 0.02, 0.98).astype(np.float32)


def run_batch(name, X, b, fid, guard, gib):
    torch.xpu.empty_cache()
    torch.xpu.reset_peak_memory_stats()
    torch.xpu.synchronize()
    t = time.time()
    u = torch.tensor(X, dtype=torch.float32, device=dev)
    g = torch.Generator(device=dev)
    g.manual_seed(1000 + b)
    o = chem.simulate(u, g, fid, record=True)
    chem.paired_control(u, 1000 + b, fid, o)
    chem.paired_control_l2(u, 1000 + b, fid, o)
    torch.xpu.synchronize()
    sec = time.time() - t
    del u, o
    torch.xpu.empty_cache()
    guard.check()
    return {"set": name, "points": len(X), "seconds": round(sec, 1), "s_per_run": round(sec / len(X), 4),
            "peak_xpu_alloc_GiB": gib(torch.xpu.max_memory_allocated()),
            "peak_xpu_reserved_GiB": gib(torch.xpu.max_memory_reserved()),
            "xpu_reserved_after_GiB": gib(torch.xpu.memory_reserved()),
            "rss_after_GiB": gib(rss()), "cgroup_current_GiB": gib(cgroup_file("memory.current")),
            "cgroup_peak_GiB": gib(cgroup_file("memory.peak"))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "preflight_mem.json"))
    ap.add_argument("--sizes", default="512", help="batch sizes to test, comma separated")
    ap.add_argument("--fid", default="screen")
    ap.add_argument("--mem-max-gib", type=float, default=24.0)
    ap.add_argument("--limit", type=float, default=0.75, help="pass: peak <= limit x cap")
    ap.add_argument("--pair", default="AS2U")
    ap.add_argument("--job", default="")
    ap.add_argument("--sets", default="random,dense,extreme,random,random,random",
                    help="point sets per size, in order (growth is checked between the first and last)")
    a = ap.parse_args()
    chem.PAIR = a.pair
    guard = thermal.Guard(job=f"preflight {a.job}")
    cap_xpu = chem.XPU_FRACTION * torch.xpu.get_device_properties(0).total_memory
    gib = lambda x: None if x is None else round(x / 2 ** 30, 3)
    res = {"job": a.job, "fidelity": a.fid, "xpu_fraction": chem.XPU_FRACTION, "cap_xpu_GiB": gib(cap_xpu),
           "mem_max_GiB": a.mem_max_gib, "limit": a.limit,
           "pair": a.pair, "chem_p31_md5": os.popen(f"md5sum {os.path.join(HERE, 'chem_p31.py')}").read().split()[0],
           "sizes": {}}
    b = 0
    for B in [int(x) for x in a.sizes.split(",")]:
        rows = []
        try:
            seen = {}
            for name in a.sets.split(","):
                k = seen.get(name, 0)
                seen[name] = k + 1
                rows.append(run_batch(name, point_set(name, B, k), b, a.fid, guard, gib))
                b += 1
                print(B, rows[-1], flush=True)
                res["sizes"][str(B)] = {"batches": rows}
                ckpt.atomic_json(res, a.out)
        except torch.OutOfMemoryError as e:
            res["sizes"][str(B)] = {"batches": rows, "oom": str(e)[:300]}
            ckpt.atomic_json(res, a.out)
            print(B, "OOM", flush=True)
            torch.xpu.empty_cache()
            break
        worst_res = max(x["peak_xpu_reserved_GiB"] for x in rows)
        worst_host = max((x["cgroup_peak_GiB"] or x["rss_after_GiB"]) for x in rows)
        # growth: host memory after each later batch against the first (P3-1: >= 5 batches after the first)
        r0 = rows[0]
        growth = (max(x["cgroup_current_GiB"] for x in rows[1:]) / r0["cgroup_current_GiB"] - 1) \
            if r0["cgroup_current_GiB"] and len(rows) > 1 else None
        d = res["sizes"][str(B)]
        d.update({"worst_peak_xpu_reserved_GiB": worst_res, "worst_host_peak_GiB": worst_host,
                  "s_per_run_random": round(sum(x["s_per_run"] for x in rows if x["set"] == "random")
                                            / max(1, sum(x["set"] == "random" for x in rows)), 4)
                                      if any(x["set"] == "random" for x in rows) else None,
                  "s_per_run_worst": max(x["s_per_run"] for x in rows), "host_growth_random": growth,
                  "pass": bool(worst_res <= a.limit * cap_xpu / 2 ** 30 and worst_host <= a.limit * a.mem_max_gib
                               and (growth or 0.0) <= 0.10)})
        ckpt.atomic_json(res, a.out)
        print("SIZE", B, {k: d[k] for k in ("worst_peak_xpu_reserved_GiB", "worst_host_peak_GiB", "s_per_run_random",
                                            "s_per_run_worst", "host_growth_random", "pass")}, flush=True)
    ok = [int(k) for k, v in res["sizes"].items() if v.get("pass")]
    if ok:
        best = min(ok, key=lambda k: res["sizes"][str(k)]["s_per_run_random"])
        res["chosen_batch"] = best
    res["pass"] = bool(ok)
    res["thermal"] = guard.summary()
    ckpt.atomic_json(res, a.out)
    print("PREFLIGHT", "PASS" if res["pass"] else "FAIL", "chosen batch", res.get("chosen_batch"), flush=True)


if __name__ == "__main__":
    main()
