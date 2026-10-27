"""Stage 4c-3 map (copy of run4c2_map.py, md5 d45d1f57, adapted): joint map of E1, E2, L1 and
L2 (= E1 & L2_strict; linked to E1 only) over all 42 dimensions of chem4c3. After each batch the
runs with L2_loose are rerun with copy_random and the same batch seed (chem4c3.paired_control_l2),
as the runs with L1_loose are rerun with p0 = 0. Detector defaults are those adopted in 4c-2
(8 RAF sites, refresh every step). The 4c-2 description follows.
Stage 4c-2 map (copy of run4c1_map.py adapted): joint map of E1, E2 (= E1 & E2_L8) and
L1 (= E1 & L1_strict; NOT nested in E2) over all 37 dimensions of chem4c2, COMPASS
engine4_soft (FAST) only (no Sobol arm), plus ground truth, importance, partial dependence,
slices and composition summaries. L1_strict needs a paired run without catalysis: after each
batch, the runs with L1_loose are rerun with p0 = 0 and the same batch seed (second pass,
chem4c2.paired_control). Crash-safe: every finished batch is one JSON line (fsync),
every COMPASS round's selection is checkpointed before it is simulated, `--resume` skips
finished work, and seeds come from (config, phase, rep, round, batch), so a resumed run gives
the same simulations as an uninterrupted one. A thermal guard (thermal.py) runs between batches.

  python run4c1_map.py --out-dir map4c1 --n-test 2000 --gt-seeds 32 --budget 32000 --reps 1
  python run4c1_map.py --out-dir map4c1 ... --resume
"""
import argparse, json, math, os, sys, time
import numpy as np
import torch
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "lot-run", "stage1b")))
import chem4c3 as chem  # noqa: E402
import ckpt  # noqa: E402
import thermal  # noqa: E402
import bench1b  # noqa: E402
import surrogate  # noqa: E402

bench1b.FAST = True
PHASES = {"gt": 1, "sobol": 2, "engine": 3, "test": 9}
SLICES = (("T", "depth"), ("pH", "Mg"), ("p0", "kappa"), ("k_t", "eps"), ("T", "T_m8"))
# --focus: the map's unit cube is mapped affinely onto these u sub-intervals of chem4c2.decode
# (p0 1e-2..1, kappa 1e4..1e6, L_cat 3-5); all other dimensions keep their full ranges.
FOCUS_RANGES = {"p0": (4 / 6, 1.0), "kappa": (0.6, 1.0), "L_cat": (0.0, 0.4999)}
FOCUS = None


def to_model_u(x):
    """Map coordinates -> chem4c2 u (identity unless --focus)."""
    if not FOCUS:
        return x
    x = x.copy()
    for nm, (lo, hi) in FOCUS.items():
        d = chem.NAMES.index(nm)
        x[:, d] = lo + (hi - lo) * x[:, d]
    return x


def seed_of(phase, rep, rnd, batch):
    return ((PHASES[phase] * 1000 + rep) * 10000 + rnd) * 100000 + batch


# ---------- JSON-lines store ----------

def load_lines(path):
    """Valid lines of a JSONL file; a torn last line (power cut mid-write) is dropped and the
    file rewritten atomically without it."""
    if not os.path.exists(path):
        return []
    good, bad = [], False
    with open(path) as f:
        for ln in f:
            try:
                good.append(json.loads(ln))
            except json.JSONDecodeError:
                bad = True
    if bad:
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            for d in good:
                f.write(json.dumps(d) + "\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    return good


def append_line(path, d):
    with open(path, "a") as f:
        f.write(json.dumps(d) + "\n")
        f.flush()
        os.fsync(f.fileno())


def run_summary(o, i, strict, ab_cat, ab_ctl, strict2, s_main, s_ctl):
    h = o["length_hist"][i]
    n = float(h[2:].sum())
    mass = o["chain_mass_by_chem"][i]
    nan2none = lambda x: None if math.isnan(float(x)) else float(x)
    return {
        "L1_strict": bool(strict[i]), "L1_loose": bool(o["L1_loose"][i]), "RAF_final": bool(o["RAF_final"][i]),
        "RAF_size": float(o["RAF_size"][i]), "RAF_components": float(o["RAF_components"][i]),
        "RAF_span": float(o["RAF_span_turnovers"][i]), "RAF_ab_cat": float(ab_cat[i]), "RAF_ab_ctl": float(ab_ctl[i]),
        "cat_share_run": float(o["polymer_catalysed_share"][i]),
        "cat_share_RAF_site": nan2none(o["RAF_site_catalysed_share"][i]),
        "RAF_product_mass": float(o["RAF_product_mass"][i]),
        "RAF_catalyst_L": o["RAF_catalyst_L"][i], "RAF_nuc_share": nan2none(o["RAF_nuc_share"][i]),
        "L1_minus_E2_time": nan2none(o["L1_minus_E2_time"][i]),
        "E1": bool(o["E1"][i]), "E2": bool(o["E2"][i]), "E2_L4": bool(o["E2_L4"][i]),
        "E2_L16": bool(o["E2_L16"][i]), "liquid": bool(o["liquid"][i]),
        "reached_end": bool(o["reached_end"][i]), "stationary": bool(o["stationary"][i]),
        "conc": float(o["conc_mstar"][i]), "frac8": float(o["frac_sites_L8"][i]),
        "mean_len": float((h[2:] * torch.arange(2, chem.L_MAX + 1)).sum() / n) if n else 0.0,
        "n_chains": n, "nuc_share_chains": float(mass[0] / mass.sum()) if mass.sum() > 0 else None,
        "distinct": [float(v) for v in o["distinct_seq_by_class"][i]],
        "entropy": [round(float(v), 3) for v in o["entropy_bits_by_class"][i]],
        "hist": [int(v) for v in h[:33]] + [int(h[33:].sum())], "overflow": float(o["overflow_events"][i]),
        "t_first_E2": None if math.isnan(float(o["t_first_E2"][i])) else float(o["t_first_E2"][i]),
        "L2_strict": bool(strict2[i]), "L2_loose": bool(o["L2_loose"][i]),
        "L2_span": float(o["L2_span_turnovers"][i]), "L2_share": nan2none(s_main[i]), "L2_share_ctl": nan2none(s_ctl[i]),
        "L2_share_last5": float(o["L2_share_last5"][i]), "L2_top_len": float(o["L2_top_lineage_len"][i]),
        "t_first_L2": nan2none(o["t_first_L2"][i]), "L2_order_vs_L1": nan2none(o["L2_order_vs_L1"][i]),
        "copies": float(o["copies"][i]), "copy_len": float(o["copy_mean_len"][i]),
        "copy_err": float(o["copy_error_rate"][i]), "copies_q4": float(o["copies_q4"][i]),
        "paired_share": float(o["paired_template_share"][i]),
    }


def run_points(path, phase, rep, rnd, X, spp, batch, dev, fid, guard, done_keys):
    """Simulate points X [K,D], each spp times, in batches; one JSON line per batch."""
    total = len(X) * spp
    nb = (total + batch - 1) // batch
    for b in range(nb):
        if (rnd, b) in done_keys:
            continue
        idx = np.arange(b * batch, min(total, (b + 1) * batch))
        pts = idx // spp
        gen = torch.Generator(device=dev)
        gen.manual_seed(seed_of(phase, rep, rnd, b))
        t = time.time()
        u = torch.tensor(to_model_u(X[pts]), dtype=torch.float32, device=dev)
        o = chem.simulate(u, gen, fid, record=True)
        strict, ab_cat, ab_ctl = chem.paired_control(u, seed_of(phase, rep, rnd, b), fid, o)
        strict2, s_main, s_ctl = chem.paired_control_l2(u, seed_of(phase, rep, rnd, b), fid, o)
        recs = [run_summary(o, i, strict, ab_cat, ab_ctl, strict2, s_main, s_ctl) for i in range(len(idx))]
        append_line(path, {"round": rnd, "batch": b, "i0": int(idx[0]), "n": len(idx),
                           "seconds": time.time() - t, "runs": recs})
        done_keys.add((rnd, b))
        del o, strict, ab_cat, ab_ctl, strict2, s_main, s_ctl, u
        torch.xpu.empty_cache()
        guard.check()


def gather(lines, rnd=None):
    ls = sorted((d for d in lines if rnd is None or d["round"] == rnd), key=lambda d: (d["round"], d["i0"]))
    return [r for d in ls for r in d["runs"]]


LEVELS = ("E1", "E2", "L1", "L2")


def outcomes(recs):
    e1 = np.array([r["E1"] for r in recs], np.float32)
    e2 = np.array([r["E1"] and r["E2"] for r in recs], np.float32)
    l1 = np.array([r["E1"] and r["L1_strict"] for r in recs], np.float32)
    l2 = np.array([r["E1"] and r["L2_strict"] for r in recs], np.float32)
    return e1, e2, l1, l2


# ---------- ground truth ----------

def ground_truth(a, dev, fid, guard):
    path = os.path.join(a.out_dir, "gt.jsonl")
    X = np.random.default_rng(12345).random((a.n_test, chem.DIM)).astype(np.float32)
    lines = load_lines(path)
    run_points(path, "gt", 0, 0, X, a.gt_seeds, a.batch, dev, fid, guard, {(d["round"], d["batch"]) for d in lines})
    recs = gather(load_lines(path))
    ps = [o.reshape(a.n_test, a.gt_seeds).mean(1) for o in outcomes(recs)]
    return X, ps, recs


# ---------- methods ----------

def m_sobol(a, rep, dev, fid, guard):
    path = os.path.join(a.out_dir, f"sobol_rep{rep}.jsonl")
    X = qmc.Sobol(d=chem.DIM, scramble=True, seed=1000 + rep).random(a.budget).astype(np.float32)
    lines = load_lines(path)
    run_points(path, "sobol", rep, 0, X, 1, a.batch, dev, fid, guard, {(d["round"], d["batch"]) for d in lines})
    return X, gather(load_lines(path)), {}


def m_engine(a, rep, dev, fid, guard, explore=0.3, n_init=2048, ens=5, pool_bits=19):
    """engine4_soft on the joint map: score = sqrt(sd_E1^2 + sd_E2^2 + sd_L1^2 + sd_L2^2) of bootstrap ensembles,
    70% of each round sampled with probability ~ score^2, 30% sequential Sobol."""
    path = os.path.join(a.out_dir, f"engine_rep{rep}.jsonl")
    cpath = os.path.join(a.out_dir, f"engine_rep{rep}_rounds.json")
    pool = qmc.Sobol(d=chem.DIM, scramble=True, seed=1000 + rep).random(1 << pool_bits).astype(np.float32)
    state = ckpt.load_json(cpath) or {"rounds": [], "next_seq": 0}
    lines = load_lines(path)
    done = {(d["round"], d["batch"]) for d in lines}
    used = np.zeros(len(pool), bool)
    for sel in state["rounds"]:
        used[sel] = True
    pool_t = None
    rnd = 0
    while True:
        total = sum(len(s) for s in state["rounds"][:rnd])
        if total >= a.budget:
            break
        if rnd >= len(state["rounds"]):
            nxt = state["next_seq"]
            def take_seq(m):
                nonlocal nxt
                idx = []
                while len(idx) < m and nxt < len(pool):
                    if not used[nxt]:
                        idx.append(nxt)
                    nxt += 1
                return idx
            if rnd == 0:
                sel = take_seq(min(n_init, a.budget))
            else:
                batch = min(max(1024, total // 3), a.budget - total)
                seq_idx = take_seq(int(batch * explore))
                used[seq_idx] = True
                prev = [r for k in range(rnd) for r in gather(lines, k)]
                xs = np.concatenate([pool[state["rounds"][k]] for k in range(rnd)])
                if pool_t is None:
                    pool_t = torch.tensor(pool, device=dev)
                x = torch.tensor(xs, device=dev)
                sd2 = np.zeros(len(pool))
                for j, s in enumerate(outcomes(prev)):
                    st = torch.tensor(s, device=dev)
                    preds = surrogate.predict_multi(surrogate.fit_es_multi(
                        x, st, torch.ones_like(st), [len(x)] * ens, seed=rep * 100 + rnd * 10 + j, hidden=128,
                        max_steps=3000, bootstrap=True), pool_t)
                    sd2 += preds.std(0).cpu().numpy().astype(np.float64) ** 2
                score = np.sqrt(sd2)
                score[used] = -np.inf
                m = batch - len(seq_idx)
                rng = np.random.default_rng(3000 + rep * 1000 + rnd)
                w = np.where(np.isfinite(score) & (score > 0), score, 0.0) ** 2
                k = min(m, int((w > 0).sum()))
                act = rng.choice(len(pool), size=k, replace=False, p=w / w.sum()) if k > 0 else np.zeros(0, int)
                if k < m:
                    rest = np.argsort(-np.where(w > 0, -np.inf, score))[:m - k]
                    act = np.concatenate([act, rest])
                sel = list(map(int, seq_idx)) + list(map(int, act))
            used[sel] = True
            state["rounds"].append([int(i) for i in sel])
            state["next_seq"] = nxt
            ckpt.atomic_json(state, cpath)          # selection is saved before it is simulated
        sel = state["rounds"][rnd]
        run_points(path, "engine", rep, rnd, pool[sel], 1, a.batch, dev, fid, guard, done)
        lines = load_lines(path)
        rnd += 1
    X = np.concatenate([pool[s] for s in state["rounds"]])[:a.budget]
    recs = [r for k in range(len(state["rounds"])) for r in gather(lines, k)][:a.budget]
    return X, recs, {"rounds": len(state["rounds"])}


# ---------- analysis ----------

def accuracy_curves(X, recs, test_x, gt, dev, a, rep):
    out = {}
    for j, name in enumerate(LEVELS):
        s = outcomes(recs)[j]
        tp = gt[j]
        hard = ((tp >= 0.1) & (tp <= 0.4)) | ((tp >= 0.6) & (tp <= 0.9))
        out[name] = bench1b.evaluate(X, s, np.ones(len(s), np.float32), test_x, tp, hard, dev, a.budget, rep)
        out[name + "_hard_n"] = int(hard.sum())
    return out


def final_maps(X, recs, test_x, gt, dev):
    """Importance, partial dependence and slices from a map fitted on the given points."""
    res = {}
    rng = np.random.default_rng(5)
    bg = torch.tensor(rng.random((2000, chem.DIM)).astype(np.float32), device=dev)
    tq = torch.tensor(test_x, device=dev)
    xt = torch.tensor(X, device=dev)
    for j, name in enumerate(LEVELS):
        s = torch.tensor(outcomes(recs)[j], device=dev)
        net = surrogate.fit_es(xt, s, torch.ones_like(s), seed=99, hidden=256)
        base_p = surrogate.predict(net, tq).cpu().numpy()
        lab = gt[j] > 0.5
        acc0 = float(((base_p > 0.5) == lab).mean())
        imp = {}
        for d, nm in enumerate(chem.NAMES):
            q = tq.clone()
            q[:, d] = q[torch.randperm(len(q), generator=torch.Generator().manual_seed(d)).to(dev), d]
            p = surrogate.predict(net, q).cpu().numpy()
            imp[nm] = round(acc0 - float(((p > 0.5) == lab).mean()), 4)
        grid = np.linspace(0, 1, 13)
        pd = {}
        for d, nm in enumerate(chem.NAMES):
            vals = []
            for g in grid:
                q = bg.clone()
                q[:, d] = float(g)
                vals.append(round(float(surrogate.predict(net, q).mean()), 4))
            pd[nm] = vals
        slices = {}
        for (dx, dy) in SLICES:
            ix, iy = chem.NAMES.index(dx), chem.NAMES.index(dy)
            gx = np.linspace(0, 1, 21)
            P = np.zeros((21, 21))
            for yi, gy in enumerate(gx):
                for xi, gxv in enumerate(gx):
                    q = bg[:500].clone()
                    q[:, ix] = float(gxv)
                    q[:, iy] = float(gy)
                    P[yi, xi] = float(surrogate.predict(net, q).mean())
            slices[f"{dx}_x_{dy}"] = {"x": dx, "y": dy, "u": gx.round(3).tolist(), "p": P.round(3).tolist()}
        res[name] = {"final_map_accuracy": acc0, "importance": imp, "pd_u": grid.round(3).tolist(),
                     "partial_dependence": pd, "slices": slices}
    return res


def composition_summary(recs, X, spp):
    """Composition of E2-positive runs (E1 & E2), from the ground-truth runs."""
    pos = [(i, r) for i, r in enumerate(recs) if r["E1"] and r["E2"]]
    if not pos:
        return {"n_positive_runs": 0}
    arr = lambda key: np.array([r[key] for _, r in pos if r[key] is not None], float)
    xs = np.array([X[i // spp] for i, _ in pos])
    nuc = arr("nuc_share_chains")
    ml = arr("mean_len")
    hist = np.array([r["hist"] for _, r in pos], float)
    histn = hist / np.maximum(hist[:, 2:].sum(1, keepdims=True), 1)
    q = lambda v: [round(float(x), 4) for x in np.percentile(v, [10, 25, 50, 75, 90])] if len(v) else None
    out = {"n_positive_runs": len(pos), "nuc_share_in_chains_pct10_25_50_75_90": q(nuc),
           "share_runs_mostly_nucleotide(>0.9)": float((nuc > 0.9).mean()) if len(nuc) else None,
           "share_runs_mostly_peptide(<0.1)": float((nuc < 0.1).mean()) if len(nuc) else None,
           "mean_len_pct10_25_50_75_90": q(ml),
           "mean_length_hist": histn.mean(0).round(5).tolist(),
           "distinct_by_class_median": np.median(np.array([r["distinct"] for _, r in pos]), 0).tolist(),
           "entropy_bits_by_class_median": np.median(np.array([r["entropy"] for _, r in pos]), 0).tolist(),
           "length_classes": [[2, 3], [4, 7], [8, 15], [16, 32]]}
    by = {}
    for nm, edges in (("T", [0, 1 / 3, 2 / 3, 1.0001]), ("x_nuc", [0, 1 / 3, 2 / 3, 1.0001]),
                      ("mineral_type", [0, 1 / 3, 2 / 3, 1.0001])):
        d = chem.NAMES.index(nm)
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (xs[:, d] >= lo) & (xs[:, d] < hi)
            sub = [r for (_, r), k in zip(pos, m) if k]
            nv = [r["nuc_share_chains"] for r in sub if r["nuc_share_chains"] is not None]
            rows.append({"u_range": [round(lo, 3), round(min(hi, 1), 3)], "n": len(sub),
                         "nuc_share_median": float(np.median(nv)) if nv else None,
                         "mean_len_median": float(np.median([r["mean_len"] for r in sub])) if sub else None})
        by[nm] = rows
    out["by_condition"] = by
    return out


def l1_summary(recs, X, spp):
    """L1 (strict) positive GT runs: RAF size, catalyst lengths, catalysed share, product mass,
    nucleotide share, overlap with E2, timing relative to E2."""
    pos = [(i, r) for i, r in enumerate(recs) if r["E1"] and r["L1_strict"]]
    loose = [r for r in recs if r["E1"] and r["L1_loose"]]
    if not pos:
        return {"n_positive_runs": 0, "n_loose_runs": len(loose)}
    arr = lambda key: np.array([r[key] for _, r in pos if r[key] is not None], float)
    q = lambda v: [round(float(x), 4) for x in np.percentile(v, [10, 25, 50, 75, 90])] if len(v) else None
    catL = [x for _, r in pos for x in r["RAF_catalyst_L"]]
    xs = np.array([X[i // spp] for i, _ in pos])
    out = {"n_positive_runs": len(pos), "n_loose_runs": len(loose),
           "share_with_E2": float(np.mean([r["E2"] for _, r in pos])),
           "RAF_size_pct": q(arr("RAF_size")), "RAF_span_turnovers_pct": q(arr("RAF_span")),
           "catalyst_length_pct": q(np.array(catL, float)) if catL else None,
           "cat_share_RAF_site_pct": q(arr("cat_share_RAF_site")), "RAF_product_mass_pct": q(arr("RAF_product_mass")),
           "RAF_nuc_share_pct": q(arr("RAF_nuc_share")), "mean_len_pct": q(arr("mean_len")),
           "abundance_cat_over_control_pct": q(arr("RAF_ab_cat") / np.maximum(arr("RAF_ab_ctl"), 1.0)),
           "L1_minus_E2_time_pct": q(arr("L1_minus_E2_time"))}
    by = {}
    for nm in ("T", "pH", "x_nuc", "mineral_type", "p0"):
        d = chem.NAMES.index(nm)
        by[nm] = [int(((xs[:, d] >= lo) & (xs[:, d] < hi)).sum()) for lo, hi in ((0, 1 / 3), (1 / 3, 2 / 3), (2 / 3, 1.0001))]
    out["count_by_condition_thirds"] = by
    return out


def l2_summary(recs, X, spp):
    """L2 (strict) positive GT runs: lineage persistence, information share vs control, copying,
    overlap with E2 and L1, order of first L1 and first L2."""
    pos = [(i, r) for i, r in enumerate(recs) if r["E1"] and r["L2_strict"]]
    loose = [r for r in recs if r["E1"] and r["L2_loose"]]
    if not pos:
        return {"n_positive_runs": 0, "n_loose_runs": len(loose)}
    arr = lambda key: np.array([r[key] for _, r in pos if r[key] is not None], float)
    q = lambda v: [round(float(x), 4) for x in np.percentile(v, [10, 25, 50, 75, 90])] if len(v) else None
    xs = np.array([X[i // spp] for i, _ in pos])
    out = {"n_positive_runs": len(pos), "n_loose_runs": len(loose),
           "share_with_E2": float(np.mean([r["E2"] for _, r in pos])),
           "share_with_L1": float(np.mean([r["L1_strict"] for _, r in pos])),
           "L2_span_pct": q(arr("L2_span")), "L2_share_pct": q(arr("L2_share")),
           "L2_share_minus_ctl_pct": q(arr("L2_share") - np.nan_to_num(arr("L2_share_ctl"))),
           "top_lineage_len_pct": q(arr("L2_top_len")), "copies_pct": q(arr("copies")),
           "copy_len_pct": q(arr("copy_len")), "copy_err_pct": q(arr("copy_err")),
           "paired_share_pct": q(arr("paired_share")), "order_vs_L1_pct": q(arr("L2_order_vs_L1"))}
    by = {}
    for nm in ("T", "depth", "k_t", "eps", "T_m8", "L_t"):
        d = chem.NAMES.index(nm)
        by[nm] = [int(((xs[:, d] >= lo) & (xs[:, d] < hi)).sum()) for lo, hi in ((0, 1 / 3), (1 / 3, 2 / 3), (2 / 3, 1.0001))]
    out["count_by_condition_thirds"] = by
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "map4c3"))
    ap.add_argument("--n-test", type=int, default=1000)
    ap.add_argument("--gt-seeds", type=int, default=16)
    ap.add_argument("--budget", type=int, default=32000)
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--batch", type=int, default=512)
    ap.add_argument("--fidelity", default="screen")
    ap.add_argument("--methods", default="engine")
    ap.add_argument("--skip-analysis", action="store_true")
    ap.add_argument("--focus", action="store_true", help="restrict p0, kappa, L_cat (FOCUS_RANGES)")
    ap.add_argument("--detector", default="", help="JSON file with adopted detector settings "
                    "(ns_raf, refresh, cat_after_e1), e.g. sens4c2/result.json")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    global FOCUS
    det = {"ns_raf": chem.NS_RAF, "refresh": chem.REFRESH, "cat_after_e1": chem.CAT_AFTER_E1}
    if a.detector:
        det.update(ckpt.load_json(a.detector)["adopt"])
    chem.NS_RAF, chem.REFRESH, chem.CAT_AFTER_E1 = int(det["ns_raf"]), int(det["refresh"]), bool(det["cat_after_e1"])
    if a.focus:
        FOCUS = dict(FOCUS_RANGES)
    os.makedirs(a.out_dir, exist_ok=True)
    cfg = {k: getattr(a, k) for k in ("n_test", "gt_seeds", "budget", "reps", "batch", "fidelity", "methods")}
    if a.focus:
        cfg["focus"] = True
        cfg["focus_u_ranges"] = {k: list(v) for k, v in FOCUS_RANGES.items()}
        cfg["focus_physical"] = {"p0": [1e-2, 1.0], "kappa": [1e4, 1e6], "L_cat": [3, 5]}
    if a.detector:
        cfg["detector"] = a.detector
    cfg["detector_settings"] = det
    cfg_path = os.path.join(a.out_dir, "config.json")
    old = ckpt.load_json(cfg_path)
    if old is not None and not a.resume:
        sys.exit(f"{a.out_dir} already has a run; use --resume or a new --out-dir")
    if old is not None and "detector_settings" not in old:        # folder started before the knob existed
        cfg.pop("detector_settings")
    if old is not None and old != cfg:
        sys.exit(f"settings differ from the saved run: {old} vs {cfg}")
    ckpt.atomic_json(cfg, cfg_path)
    dev = torch.device(a.device)
    guard = thermal.Guard(job=f"4c3 {os.path.basename(a.out_dir)}")
    t0 = time.time()
    test_x, ps, gt_recs = ground_truth(a, dev, a.fidelity, guard)
    p1, p2, p3, p4 = ps
    loose = np.array([r["E1"] and r["L1_loose"] for r in gt_recs], np.float32).reshape(a.n_test, a.gt_seeds).mean(1)
    loose2 = np.array([r["E1"] and r["L2_loose"] for r in gt_recs], np.float32).reshape(a.n_test, a.gt_seeds).mean(1)
    print(f"ground truth {a.n_test} x {a.gt_seeds}: E1 life {np.mean(p1 > 0.5):.3f}, E2 life {np.mean(p2 > 0.5):.3f}, "
          f"L1 life {np.mean(p3 > 0.5):.3f}, L2 life {np.mean(p4 > 0.5):.3f} ({time.time() - t0:.0f}s)", flush=True)
    res = {"config": cfg, "dims": chem.NAMES, "inactive": chem.INACTIVE,
           "life_fraction": {"E1": float(np.mean(p1 > 0.5)), "E2": float(np.mean(p2 > 0.5)),
                             "L1_strict": float(np.mean(p3 > 0.5)), "L1_loose": float(np.mean(loose > 0.5)),
                             "L1_and_E2": float(np.mean((p3 > 0.5) & (p2 > 0.5))),
                             "L1_without_E2": float(np.mean((p3 > 0.5) & (p2 <= 0.5))),
                             "L2_strict": float(np.mean(p4 > 0.5)), "L2_loose": float(np.mean(loose2 > 0.5)),
                             "L2_and_E2": float(np.mean((p4 > 0.5) & (p2 > 0.5))),
                             "L2_and_L1": float(np.mean((p4 > 0.5) & (p3 > 0.5)))},
           "gt_run_rates": {"L1_strict": float(np.mean([r["E1"] and r["L1_strict"] for r in gt_recs])),
                            "L1_loose": float(np.mean([r["E1"] and r["L1_loose"] for r in gt_recs])),
                            "L1_and_E2": float(np.mean([r["E1"] and r["L1_strict"] and r["E2"] for r in gt_recs])),
                            "L2_strict": float(np.mean([r["E1"] and r["L2_strict"] for r in gt_recs])),
                            "L2_loose": float(np.mean([r["E1"] and r["L2_loose"] for r in gt_recs]))},
           "gt_flags": {k: float(np.mean([r[k] for r in gt_recs])) for k in ("liquid", "reached_end", "stationary")},
           "gt_nesting_violation_E2raw_not_E1": float(np.mean([r["E2"] and not r["E1"] for r in gt_recs])),
           "gt_overflow_run_share": float(np.mean([r["overflow"] > 0 for r in gt_recs])),
           "gt_ms_per_run": 1000 * sum(d["seconds"] for d in load_lines(os.path.join(a.out_dir, "gt.jsonl")))
           / (a.n_test * a.gt_seeds),
           "methods": {}}
    final_data = None
    for name in a.methods.split(","):
        res["methods"][name] = []
        for rep in range(a.reps):
            t = time.time()
            X, recs, extra = (m_sobol if name == "sobol" else m_engine)(a, rep, dev, a.fidelity, guard)
            tm = time.time() - t
            curves = accuracy_curves(X, recs, test_x, ps, dev, a, rep)
            res["methods"][name].append({"curves": curves, "seconds_this_session": tm, **extra})
            print(f"{name} rep {rep}: E1 {curves['E1'][-1]} E2 {curves['E2'][-1]} L1 {curves['L1'][-1]} "
                  f"L2 {curves['L2'][-1]} ({tm:.0f}s)",
                  flush=True)
            if name == "engine" and rep == 0:
                final_data = (X, recs)
            guard.check(force=True)
    if not a.skip_analysis:
        if final_data is None:
            final_data = (X, recs)
        res["maps"] = final_maps(*final_data, test_x, ps, dev)
        res["composition_E2_positive"] = composition_summary(gt_recs, test_x, a.gt_seeds)
        res["composition_L1_positive"] = l1_summary(gt_recs, test_x, a.gt_seeds)
        res["composition_L2_positive"] = l2_summary(gt_recs, test_x, a.gt_seeds)
    res["thermal"] = guard.summary()
    res["total_seconds_this_session"] = time.time() - t0
    ckpt.atomic_json(res, os.path.join(a.out_dir, "result.json"))
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
