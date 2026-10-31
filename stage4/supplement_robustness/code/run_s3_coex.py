"""P3-1b section 3 (Review 4, R4.3): the 4c-5 point study rerun with chem4c5b per variant (base, A0, A2, A2a, A3, A4).
Copy of lot-s4/run4c5_5a.py (not modified); changes: chem4c5b, variant switches on chem4c4b, a variant override passed
to the main run and to both controls, one JSONL file per variant, the round-cut / MULTI_START records stored per run.
    python run_s3_coex.py out_dir --variants base,A0,A2,A2a,A3,A4 [--batch 256]
--- run4c5_5a description follows ---
Stage 4c-5 step 5a: point study (spec docs/stage4c5_spec.md §4, commit 05017cd). 16 seeds per point,
chem4c5 (chem4c4 v2 dynamics + lineage record + L3 / L4 detectors), screen fidelity, with the paired beta = 0
and random-copy controls (same constants and seed), at:
  * engine12538 (4c-4 joint point, confirm4c4.json);
  * the 4 membrane-compatible 9a points (p1558, p3709, C_engine10653, rep1_1926), each with its first 8
    window-compatible membrane sets (scan4c4mem_sets.sets: the 9a draws; per-set outcomes of the scan were
    not stored, so the joint rate of each set is reported here);
  * the confirmed joint hits of the 9c map (confirm4c4v2.json).
Gate as in the 5b runner: E1 & L2_long (L2_strict from paired_control_l2, copy length >= L_H) & E3_strict;
L3 strict / L4 with the paired random-copy inheritance condition (chem4c5.paired_inherit). Controls: the
controls are scored WITHOUT the L2/E3 gate wherever evaluable (spec a9227f6, choice 19: chem4c5.score_control,
over the main run's joint window, else the whole record), with their own gate status (E1, L2_candidate,
E3_strict) reported next to the scores.
Composition (runs with L3 loose): enclosed chemistry-0 chains >= L_H per vesicle (count, length
percentiles), dominant keys (rising and at the last checkpoint, with sequences), membrane at division
(parent's membrane at the checkpoint before the division interval), generation time (window hours / mean
generations gained by the vesicles alive at the last checkpoint), maximum generation, rho(g, Pi), rho(w, z).
Batches of 256 runs (16 point-sets x 16 seeds); one JSONL line per finished batch (resume: done batches
are skipped).
    python run4c5_5a.py [out_dir] [--batch 256]
"""
import json, math, os, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
S4 = os.path.join(HERE, "..", "..", "lot-s4")
sys.path.insert(0, HERE)
import chem4c5b as chem  # noqa: E402
import chem4c4b  # noqa: E402
import ckpt, thermal  # noqa: E402
sys.path.insert(1, S4)
import scan4c4mem_sets as S9  # noqa: E402  (lot-s4, read only)

VARIANTS = {"base": ({}, {}), "A0": ({"POISSON_INV": True}, {}), "A2": ({"MULTI_START": True}, {}),
            "A2a": ({"SLOT_ROUND_FIX": True}, {}), "A3": ({}, {"f_rehyd": 1.0}), "A4": ({"DUPLEX_HYDRO_FACTOR": 0.1}, {}),
            "A2b": ({"SLOT_CUT_RANDOM": True}, {}), "A2c": ({"MULTI_START": True, "MS_MELT_FLOOR": False}, {})}
SW0 = {"MULTI_START": False, "SLOT_ROUND_FIX": False, "POISSON_INV": False, "DUPLEX_HYDRO_FACTOR": 1.0,
       "SLOT_CUT_RANDOM": False, "MS_MELT_FLOOR": True, "MS_ORDER_SLOT": False}
PERTURB_VARIANTS = ("base", "A2b")       # Review 9: L3 detector perturbation on the records in memory
PSETS = [("ident85", {"L3_IDENT": 0.85}), ("ident95", {"L3_IDENT": 0.95}), ("frac40", {"L3_FRAC": 0.4}),
         ("frac60", {"L3_FRAC": 0.6}), ("up2", {"L3_UP": 2}), ("up4", {"L3_UP": 4})]
VOV = {}


def set_variant(v):
    global VOV
    for k, x in SW0.items():
        setattr(chem4c4b, k, x)
    sw, ov = VARIANTS[v]
    for k, x in sw.items():
        setattr(chem4c4b, k, x)
    VOV = dict(ov)
    chem.PERTURB_SETS = PSETS if v in PERTURB_VARIANTS else None

dev = torch.device("xpu")
FID, NSEED, NSETS, SEED0 = "screen", 16, 8, 95000


def points():
    P = {"engine12538": json.load(open(os.path.join(S4, "confirm4c4.json")))["points"]["engine12538"]["u_model"]}
    for name in ("p1558", "p3709", "C_engine10653", "rep1_1926"):
        for j, s in enumerate(S9.sets(name)[:NSETS]):
            P[f"{name}_set{j}"] = [float(v) for v in s]
    c2 = json.load(open(os.path.join(S4, "confirm4c4v2.json")))
    for k, v in c2["points"].items():
        if v.get("both", 0) > 0:
            P[f"9c_{k}"] = v["u_model"]
    return P


def _f(x):
    x = float(x)
    return None if not math.isfinite(x) else x


def composition(o, b):
    lin = o["lineage"]; r = lin.r; last = chem.N_CHECK - 1
    k0 = int(o["joint_window_k0"][b])
    fc = lin.V.final_chains.get(b) or {}
    al = r["alive"][b, :, last]
    n_ch, Ls = [], []
    if fc:
        for s in range(lin.N):
            for v in np.nonzero(al[s])[0]:
                m = (fc["site"] == s) & (fc["slot"] == v)
                n_ch.append(int(m.sum())); Ls.extend(fc["L"][m].tolist())
    pct = lambda a: [float(x) for x in np.percentile(a, [10, 50, 90])] if len(a) else None
    # membrane at division: parent's membrane at the checkpoint before the division interval
    ev = lin.div
    sel = (ev[:, 0] // lin.N == b) & (ev[:, 4] >= 1)
    a_div = []
    for row, par, d1, d2, k in ev[sel]:
        s = int(row % lin.N)
        pv = np.nonzero(r["alive"][b, s, k - 1] & (r["vid"][b, s, k - 1] == par))[0]
        if len(pv):
            a_div.append(float(r["a"][b, s, k - 1, pv[0]]))
    g0 = r["gen"][b, :, k0][r["alive"][b, :, k0]]
    g1 = r["gen"][b, :, last][al]
    t_end = float(o["t_end"][b])
    dgen = float(g1.mean() - g0.mean()) if len(g0) and len(g1) else float("nan")
    hours = t_end * (last - k0) / chem.N_CHECK
    top = [(int(k), lin.V.key_table.get(int(k))) for k in np.unique(r["top_k"][b, :, last, :, 0][al]) if k != 0][:5]
    d = o["L3_details"][b] or {}
    return {"chains_ge_LH_per_vesicle_pct": pct(n_ch), "chain_len_pct": pct(Ls), "membrane_at_division_pct": pct(a_div),
            "generation_time_h": _f(hours / dgen) if dgen and dgen > 0 else None, "max_generation": int(g1.max()) if len(g1) else None,
            "dominant_keys_last": [(k, None if v is None else [int(v[0]), int(v[1])]) for k, v in top],
            "key_risers": d.get("key_risers"), "rho_g_Pi": _f(o["rho_g_Pi"][b]), "rho_w_z": _f(o["rho_w_z"][b]),
            "inherit_local": d.get("inherit_local")}


def main_runs(u, sd):
    g = torch.Generator(device=dev); g.manual_seed(sd)
    o = chem.simulate(u, g, FID, record=True, override=dict(VOV))
    s2, s_main, s_ctl = chem.paired_control_l2(u, sd, FID, {k: v for k, v in o.items() if k != "lineage"}, override=dict(VOV))
    gate = (o["E1"] & s2 & (o["copy_mean_len"] >= chem.L_H) & o["E3_strict"]).numpy()
    chem.PERTURB_OUT = None
    st, l4, ism, isc = chem.paired_inherit(u, sd, FID, o, override=dict(VOV),
                                           mask=[bool(gate[k] and o["L3_strict"][k]) for k in range(len(gate))])
    o["_perturb"] = chem.PERTURB_OUT
    return o, gate, s2, st, l4, ism, isc


def control_runs(u, sd, cov, k0s):
    chem.PERTURB_OUT = None
    d = chem.score_control(u, sd, FID, dict(VOV, **cov), k0s)
    d["_perturb"] = chem.PERTURB_OUT
    d["gate"] = d["E1"] & d["L2_candidate"] & d["copy_len_ok"] & d["E3_strict"]
    return d


def main():
    out_dir = sys.argv[1]
    batch = int(sys.argv[sys.argv.index("--batch") + 1]) if "--batch" in sys.argv else 256
    variants = sys.argv[sys.argv.index("--variants") + 1].split(",") if "--variants" in sys.argv else list(VARIANTS)
    os.makedirs(out_dir, exist_ok=True)
    guard = thermal.Guard(job="p31b s3 coex")
    for v in variants:
        run_variant(os.path.join(out_dir, v), v, batch, guard)


def run_variant(out_dir, variant, batch, guard):
    os.makedirs(out_dir, exist_ok=True)
    set_variant(variant)
    P = points()
    names = list(P)
    runs = [(n, s) for n in names for s in range(NSEED)]
    ckpt.atomic_json({"spec": "stage4c5_spec.md a9227f6", "fidelity": FID, "seeds": NSEED, "points": P, "batch": batch,
                      "runner": "run_s3_coex.py", "variant": variant, "switches": VARIANTS[variant][0],
                      "override": VARIANTS[variant][1], "chem4c4b_md5": os.popen(f"md5sum {os.path.join(HERE, 'chem4c4b.py')}").read().split()[0],
                      "deviation": os.environ.get("LOT_DEVIATION")}, os.path.join(out_dir, "config.json"))
    path = os.path.join(out_dir, "runs.jsonl")
    done = set()
    if os.path.exists(path):
        for line in open(path):
            done.add(json.loads(line)["batch"])
    for bi in range((len(runs) + batch - 1) // batch):
        if bi in done:
            continue
        t = time.time()
        part = runs[bi * batch:(bi + 1) * batch]
        u = torch.tensor([P[n] for n, _ in part], dtype=torch.float32, device=dev)
        sd = SEED0 + bi
        o, gate, s2, st, l4, ism, isc = main_runs(u, sd)
        guard.check()
        k0s = chem.control_k0(o)
        ctl = {c: control_runs(u, sd, cov, k0s) for c, cov in (("beta0", {"beta": 0.0}), ("randcopy", {"copy_random": True}))}
        guard.check()
        recs = []
        for k, (n, s) in enumerate(part):
            g = bool(gate[k])
            rec = {"point": n, "seed_index": s, "E1": bool(o["E1"][k]), "L2_strict": bool(s2[k]), "E3_strict": bool(o["E3_strict"][k]),
                   "gate": g, "L3_evaluable": bool(g and o["L3_evaluable"][k]), "L3_loose": bool(g and o["L3_loose"][k]),
                   "L3_strict": bool(g and st[k]), "L4": bool(g and l4[k]),
                   "inherit_share_main": _f(ism[k]), "inherit_share_randcopy": _f(isc[k]),
                   "rho_g_Pi": _f(o["rho_g_Pi"][k]), "p_g_Pi": _f(o["p_g_Pi"][k]), "rho_w_z": _f(o["rho_w_z"][k]),
                   "n_pairs": _f(o["n_pairs"][k]), "gen_span": _f(o["gen_span"][k]), "lineage_overflow": bool(o["lineage_overflow"][k]),
                   "copies": _f(o["copies"][k]), "copy_tries": _f(o["copy_tries"][k]),
                   "round_cut_site_steps": _f(o["round_cut_site_steps"][k]), "round_cut_copies_lost": _f(o["round_cut_copies_lost"][k]),
                   "copy_conc_gini": _f(o["copy_conc_gini"][k]), "copy_conc_top10": _f(o["copy_conc_top10"][k]),
                   "copy_conc_bins": _f(o["copy_conc_bins"][k]),
                   "L3_loose_flag": bool(o["L3_loose"][k]),
                   "perturb": ({nm: [_f(a[k]) for a in arr] for nm, arr in o["_perturb"].items()} if o.get("_perturb") else None),
                   **{kk: _f(o[kk][k]) for kk in o if kk.startswith("ms_")}}
            for c, d in ctl.items():
                for key in ("L3_evaluable", "L3_loose", "L3_strict", "L4", "gate", "E1", "L2_candidate", "E3_strict"):
                    rec[f"ctl_{c}_{key}"] = bool(d[key][k])
                for key in ("rho_g_Pi", "p_g_Pi", "inherit_share", "inherit_share_rc"):
                    rec[f"ctl_{c}_{key}"] = _f(d[key][k])
                rec[f"ctl_{c}_perturb"] = ({nm: [_f(a[k]) for a in arr] for nm, arr in d["_perturb"].items()} if d.get("_perturb") else None)
            if rec["L3_loose"]:
                rec["composition"] = composition(o, k)
            recs.append(rec)
        with open(path, "a") as fh:
            fh.write(json.dumps({"batch": bi, "seed": sd, "n": len(part), "seconds": round(time.time() - t, 1), "runs": recs}) + "\n")
            fh.flush(); os.fsync(fh.fileno())
        print(variant, "batch", bi, "runs", len(part), "sec", round(time.time() - t, 1),
              "gate", int(gate.sum()), "L3_loose", sum(r["L3_loose"] for r in recs), "L3_strict", sum(r["L3_strict"] for r in recs),
              "L4", sum(r["L4"] for r in recs), flush=True)
        del o, ctl, u
        torch.xpu.empty_cache()
    # summary per point
    R = [r for line in open(path) for r in json.loads(line)["runs"]]
    summ = {}
    for n in names:
        rr = [r for r in R if r["point"] == n]
        m = len(rr)
        d = {"runs": m, **{k: sum(r[k] for r in rr) for k in ("E1", "L2_strict", "E3_strict", "gate", "L3_evaluable", "L3_loose", "L3_strict", "L4")}}
        for c in ("beta0", "randcopy"):
            d[f"ctl_{c}"] = {k: sum(r[f"ctl_{c}_{k}"] for r in rr) for k in ("gate", "L3_evaluable", "L3_loose", "L3_strict", "L4")}
        vals = lambda key: [r[key] for r in rr if r["gate"] and r[key] is not None]
        d["rho_g_Pi_median"] = _f(np.median(vals("rho_g_Pi"))) if vals("rho_g_Pi") else None
        d["rho_w_z_median"] = _f(np.median(vals("rho_w_z"))) if vals("rho_w_z") else None
        d["inherit_share_main_median"] = _f(np.median(vals("inherit_share_main"))) if vals("inherit_share_main") else None
        summ[n] = d
    ckpt.atomic_json({"variant": variant, "points": summ, "thermal": guard.summary()}, os.path.join(out_dir, "result.json"))
    set_variant("base")
    print(variant, "summary written", flush=True)


if __name__ == "__main__":
    main()
