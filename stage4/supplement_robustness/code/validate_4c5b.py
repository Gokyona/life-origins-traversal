"""P3-1b Review 4 (R4.3): validation of chem4c5b / chem4c4b (the P3-1b switches ported to the 4c-5 chemistry).
  a  switches off: chem4c5b (record on) == lot-s4 chem4c5 (record on) on every output tensor of the reference and the
     final chain state, seeds 101-103 at #3709, engine12538 and #3709 with amph 0 (the 4c-5 test-1 set), FID (16, 256, 600)
  b  tests 2a, 2d, 2e of validate_p31b on chem4c4b. 2a/2d use the low-rate 4c-3-rule fixture of P3-1b test 2c (the P3-1 GC
     fixture does not exist in the 4c-4 chemistry), amph 0: 2a copy tries ON vs OFF, 2d post-slot copies ON vs OFF with
     SLOT_ROUND_FIX on, both within 3%, |z| < 3; 2e the sampler checks on chem4c4b._poisson_inv
  c  pre-flight at the 4c-5 worst case (validate4c5 test 8 set: 9a sets + engine12538, window open, slots full, record on,
     paired controls and paired_inherit), MULTI_START on, batch 256 (128 if 256 fails), 6 batches: peak XPU reserved
     <= 80% of the cap, host <= 80% of 24 GiB, no growth (<= 10%) after the first batch
    python validate_4c5b.py out.json [tests, e.g. abc]
"""
import json, math, os, pickle, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
S4 = os.path.join(HERE, "..", "..", "lot-s4")
sys.path.insert(0, HERE)
import chem4c5b as chemB  # noqa: E402  (p31b copies first)
import chem4c4b  # noqa: E402
import ckpt, thermal  # noqa: E402
sys.path.insert(1, S4)
import importlib.util  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


dev = torch.device("xpu")
guard = None
FID = (16, 256, 600)


def G():
    global guard
    if guard is None:
        guard = thermal.Guard(job="p31b chem4c5b validation")
    return guard


def gen(sd):
    g = torch.Generator(device=dev)
    g.manual_seed(sd)
    return g


def set_sw(ms=False, srf=False, pinv=False):
    chem4c4b.MULTI_START, chem4c4b.SLOT_ROUND_FIX, chem4c4b.POISSON_INV = ms, srf, pinv


def test_a():
    ref4c4 = _load("chem4c4", os.path.join(S4, "chem4c4.py"))      # lot-s4 originals, loaded by path (read-only)
    ref = _load("chem4c5", os.path.join(S4, "chem4c5.py"))
    mk = _load("make_ref4c3", os.path.join(S4, "make_ref4c3.py"))
    u3709 = [float(v) for v in mk.box2_points()[3709]] + [0.5] * 12
    u12538 = json.load(open(os.path.join(S4, "confirm4c4.json")))["points"]["engine12538"]["u_model"]
    pts = {"p3709": (u3709, {}), "engine12538": (u12538, {}), "p3709_noamph": (u3709, {"amph": 0.0})}
    set_sw()
    out, ok = {"chem4c5_md5": os.popen(f"md5sum {os.path.join(S4, 'chem4c5.py')}").read().split()[0],
               "chem4c5b_md5": os.popen(f"md5sum {os.path.join(HERE, 'chem4c5b.py')}").read().split()[0],
               "chem4c4b_md5": os.popen(f"md5sum {os.path.join(HERE, 'chem4c4b.py')}").read().split()[0]}, True
    for name, (u, ov) in pts.items():
        for seed in (101, 102, 103):
            U = torch.tensor([u] * 4, dtype=torch.float32, device=dev)
            ov2 = dict(ov, return_state=True)
            a = ref.simulate(U, gen(seed), FID, record=True, override=dict(ov2))
            G().check()
            b = chemB.simulate(U, gen(seed), FID, record=True, override=dict(ov2))
            G().check()
            bad = [k for k, v in a.items() if isinstance(v, torch.Tensor) and not (k in b and isinstance(b[k], torch.Tensor)
                   and v.shape == b[k].shape and torch.equal(torch.nan_to_num(v.cpu(), nan=-7.0), torch.nan_to_num(b[k].cpu(), nan=-7.0)))]
            st = all(torch.equal(a["state"][k].cpu(), b["state"][k].cpu()) for k in ("L", "seq", "chem"))
            la, lb = a.get("L3_details"), b.get("L3_details")
            det_eq = (repr(la) == repr(lb))
            n_t = sum(isinstance(v, torch.Tensor) for v in a.values())
            out[f"{name}_s{seed}"] = {"tensors": n_t, "mismatch": bad, "state_equal": bool(st), "L3_details_equal": det_eq}
            ok &= (not bad) and st and det_eq
            print("a", name, seed, out[f"{name}_s{seed}"], flush=True)
    out["pass"] = bool(ok)
    return out


def _cmp(on, off, per):
    a, b = on.double() / per, off.double() / per
    ma, mb = float(a.mean()), float(b.mean())
    se = math.sqrt(float(a.var()) / a.numel() + float(b.var()) / b.numel())
    return {"mean_on": ma, "mean_off": mb, "ratio": ma / mb if mb else float("nan"), "z": (ma - mb) / se if se else float("nan")}


R4 = 1e4
FIX_LOW = {"T": 35.0, "T_m8": 25.0, "k_c": 1.0, "k_t": R4, "thermo": 1.0, "depth": 0.0, "mineral_type": 0.0,
           "x_nuc": 1.0, "init_factor": 0.0, "k_h_fixed": 0.0, "k_UV": 0.0, "tau": 1e12, "eps": 1e-3, "p0": 0.0,
           "avail_time": 1e8, "P": 100.0, "feed": 10.0, "L_t": 12.0, "pH": 7.5, "clamp_mstar": 1000.0, "amph": 0.0,
           "t_end": 0.02 / R4, "init_templates": {"n": 32, "L": 12, "seq": 0}}


def _pool(ms, srf, key, seeds, nb=32):
    xs = []
    for k in range(nb):
        set_sw(ms=ms, srf=srf)
        u = torch.full((2048, chemB.DIM), 0.5, dtype=torch.float32, device=dev)
        o = chem4c4b.simulate(u, gen(seeds + 10 * k), (64, 64, 1), override=dict(FIX_LOW))
        G().check()
        xs.append(o[key])
    set_sw()
    return torch.cat(xs)


def test_b():
    res = {}
    per = 64 * 32
    on, off = _pool(True, False, "copy_tries", 9003), _pool(False, False, "copy_tries", 9004)
    d = _cmp(on, off, per)
    d["tries_on"], d["tries_off"] = float(on.sum()), float(off.sum())
    d["pass"] = abs(d["ratio"] - 1) <= 0.03 and abs(d["z"]) < 3
    res["2a_4c3rule_lowrate"] = d
    print("b 2a", d, flush=True)
    on, off = _pool(True, True, "copies", 9013), _pool(False, True, "copies", 9014)
    d = _cmp(on, off, per)
    d["pass"] = abs(d["ratio"] - 1) <= 0.03 and abs(d["z"]) < 3
    res["2d_slot_round_fix"] = d
    print("b 2d", d, flush=True)
    import validate_p31b as VB
    VB.chem = chem4c4b                         # 2e on chem4c4b._poisson_inv (same seeds and criteria as P3-1b 2e)
    e = VB.test_e()
    res["2e"] = e
    print("b 2e pass", e["pass"], flush=True)
    res["pass"] = res["2a_4c3rule_lowrate"]["pass"] and res["2d_slot_round_fix"]["pass"] and e["pass"]
    return res


def test_c(B=256):
    pf = _load("preflight_mem", os.path.join(S4, "preflight_mem.py"))
    S = _load("scan4c4mem_sets", os.path.join(S4, "scan4c4mem_sets.py"))
    pts = S.first_sets()
    pts["engine12538"] = json.load(open(os.path.join(S4, "confirm4c4.json")))["points"]["engine12538"]["u_model"]
    names = list(pts)
    ov = {"window_open": True, "phi0": 0.3, "s_shear": 10.0, "f_rehyd": 10.0}
    cap = chemB.XPU_FRACTION * torch.xpu.get_device_properties(0).total_memory
    gib = lambda x: None if x is None else round(x / 2 ** 30, 3)
    set_sw(ms=True)
    chem4c4b.MS_MELT_FLOOR = False           # Review 9 (R9-i3): the A2c settings, the heaviest multi-start variant
    rows = []
    for b in range(6):
        u = torch.tensor([pts[names[k % len(names)]] for k in range(B)], dtype=torch.float32, device=dev)
        torch.xpu.empty_cache(); torch.xpu.reset_peak_memory_stats()
        t = time.time()
        o = chemB.simulate(u, gen(800 + b), "screen", record=True, override=dict(ov))
        chemB.paired_control(u, 800 + b, "screen", {k: v for k, v in o.items() if k != "lineage"}, override=dict(ov))
        chemB.paired_control_l2(u, 800 + b, "screen", {k: v for k, v in o.items() if k != "lineage"}, override=dict(ov))
        chemB.paired_inherit(u, 800 + b, "screen", o, override=dict(ov), mask=[d is not None for d in o["L3_details"]])
        rows.append({"batch": b, "seconds": round(time.time() - t, 1), "s_per_run": round((time.time() - t) / B, 3),
                     "peak_xpu_alloc_GiB": gib(torch.xpu.max_memory_allocated()),
                     "peak_xpu_reserved_GiB": gib(torch.xpu.max_memory_reserved()),
                     "cgroup_peak_GiB": gib(pf.cgroup_file("memory.peak")), "cgroup_current_GiB": gib(pf.cgroup_file("memory.current")),
                     "rss_GiB": gib(pf.rss())})
        del o, u
        torch.xpu.empty_cache()
        G().check()
        print("c batch", rows[-1], flush=True)
    set_sw()
    chem4c4b.MS_MELT_FLOOR = True
    a = max(r["peak_xpu_reserved_GiB"] for r in rows)
    h = max((r["cgroup_peak_GiB"] or r["rss_GiB"]) for r in rows)
    growth = max((r["cgroup_current_GiB"] or r["rss_GiB"]) for r in rows[1:]) / max((rows[0]["cgroup_current_GiB"] or rows[0]["rss_GiB"]), 1e-9) - 1
    return {"B": B, "batches": rows, "worst_reserved_GiB": a, "cap_GiB": gib(cap), "worst_host_GiB": h, "host_growth": round(growth, 3),
            "pass": bool(a <= 0.8 * cap / 2 ** 30 and h <= 0.8 * 24 and growth <= 0.10)}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "validate_4c5b.json")
    tests = sys.argv[2] if len(sys.argv) > 2 else "abc"
    res = ckpt.load_json(out) or {}
    for t in tests:
        t0 = time.time()
        if t == "c":
            r = test_c(256)
            if not r["pass"]:
                res["test_c_256"] = r
                ckpt.atomic_json(res, out)
                r = test_c(128)
        else:
            r = globals()[f"test_{t}"]()
        r["seconds"] = round(time.time() - t0, 1)
        res[f"test_{t}"] = r
        ckpt.atomic_json(res, out)
        print(f"TEST {t}", "PASS" if r.get("pass") else "FAIL", flush=True)
        if not r.get("pass"):
            break
    res["thermal"] = G().summary()
    ckpt.atomic_json(res, out)


if __name__ == "__main__":
    main()
