"""Stage-3d: condition maps for L3 and L4 on the spatial model (spatial.py), with
ground truth for L1..L4 from the same simulations; ground-truth life fractions also
restricted to lk < 0.001 (space alone, almost no leak); lk_crit for L2 and L4 at the
reference from direct simulation sweeps; S_max by sigma for L3/L4 at lk = 0.

(stage-3c description follows)
Stage-3c: condition map for life level L4 (Darwinian evolution) on evolve.py,
adapted from stage-2e maps.py. Ground truth for L1..L4 comes from the same
simulations; COMPASS and Sobol maps are drawn for L4 only (L1-L3 maps are in 2e).
Extra read-offs: minimum eta for P(L4) > 0.5 at the reference, and the share of
L3 conditions that also reach L4.

(original 2e description follows)
Stage-2e: condition maps for life levels L1, L2, L3 on the protocell model.

For each level, COMPASS (engine4_soft: bench1b.m_engine3, soft selection, --fast
map-maker) and a Sobol reference each spend BUDGET full-fidelity simulations
targeting that level. One ground-truth set (N_TEST conditions x GT_SEEDS runs,
simulated once with return_stats) scores all three levels.

Analysis of each level's final COMPASS map (fit_es_multi, hidden 256, all points):
  (a) life fraction over the space (map vs ground truth)
  (b) permutation importance: drop in test accuracy when one condition is shuffled
  (c) partial dependence P(level) vs each condition (averaged over random
      backgrounds) and 2D slices S x sigma, mu x S, km x food at the 2b
      reference values for the other conditions
  (d) nesting L3 within L2 within L1 (ground truth and maps)
  (e) numbers read off the maps for the known 2a/2b behaviour
Writes maps.json (all numbers and grids); plot_maps.py draws the PNG figures.
"""
import argparse, json, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "stage2"))
sys.path.insert(2, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "stage1b")))
import bench1b  # noqa: E402
import spatial as pc  # noqa: E402
import surrogate  # noqa: E402

bench1b.model = pc
bench1b.FAST = True
EVAL_FIT = None
REF = dict(S=64, sigma=1.5, mu=0.01, K=2, d=0.005, food=1.0, km=2.0, dm=0.005, eta=0.1, lk=0.0, z1=0.5, z2=0.5)
LEVELS = (1, 2, 3, 4)
MAP_LEVELS = (3, 4)
# physical value <-> u for each condition (inverse of protocell.decode)
TO_U = {
    "S": lambda v: math.log(v / 8) / math.log(64),
    "sigma": lambda v: (v - 1) / 2,
    "mu": lambda v: math.log10(v / 1e-3) / math.log10(300),
    "K": lambda v: (v - 0.5) / 3,
    "d": lambda v: v / 0.02,
    "food": lambda v: math.log(v / 0.25) / math.log(4),
    "km": lambda v: (math.log10(v) + 1) / math.log10(300),
    "dm": lambda v: v / 0.02,
    "eta": lambda v: v / 0.2,
    "lk": lambda v: v / 0.01,
    "z1": lambda v: v,
    "z2": lambda v: v,
}
FROM_U = {
    "S": lambda u: 8 * 64 ** u, "sigma": lambda u: 1 + 2 * u, "mu": lambda u: 1e-3 * 300 ** u,
    "K": lambda u: np.minimum(1 + np.floor(u * 3), 3), "d": lambda u: 0.02 * u, "food": lambda u: 0.25 * 4 ** u,
    "km": lambda u: 0.1 * 300 ** u, "dm": lambda u: 0.02 * u, "eta": lambda u: 0.2 * u, "lk": lambda u: 0.01 * u, "z1": lambda u: u, "z2": lambda u: u,
}
SLICES = [("lk", "S"), ("lk", "sigma"), ("eta", "food"), ("S", "sigma")]


def ref_u(**over):
    v = dict(REF, **over)
    return np.array([TO_U[k](v[k]) for k in pc.NAMES], np.float32)


def predict(net, u, dev):
    return surrogate.predict_multi(net, torch.tensor(np.asarray(u, np.float32), device=dev)).cpu().numpy()[0]


def threshold_along(net, name, grid_vals, dev, rising, **over):
    """First value along grid where P crosses 0.5 (upwards if rising else downwards)."""
    u = np.stack([ref_u(**dict(over, **{name: v})) for v in grid_vals])
    p = predict(net, u, dev)
    for i in range(1, len(p)):
        if (rising and p[i - 1] < 0.5 <= p[i]) or (not rising and p[i - 1] >= 0.5 > p[i]):
            return float(grid_vals[i]), p.tolist()
    return None, p.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--budget", type=int, default=32000)
    ap.add_argument("--n-test", type=int, default=4000)
    ap.add_argument("--gt-seeds", type=int, default=64)
    ap.add_argument("--out", default="maps_sp.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    t0 = time.time()

    # ---- one ground truth for all three levels ----
    g = torch.Generator(device=dev).manual_seed(12345)
    test_x = np.random.default_rng(12345).random((a.n_test, pc.DIM)).astype(np.float32)
    xr = torch.tensor(np.repeat(test_x, a.gt_seeds, 0), device=dev)
    stats = {f"l{L}": [] for L in LEVELS}
    for i in range(0, len(xr), 1 << 15):
        _, st = pc.simulate(xr[i:i + (1 << 15)], generator=g, return_stats=True)
        for k in stats:
            stats[k].append(st[k].float().cpu())
    gt = {L: torch.cat(stats[f"l{L}"]).view(a.n_test, a.gt_seeds).mean(1).numpy() for L in LEVELS}
    print(f"ground truth {a.n_test} x {a.gt_seeds} in {time.time()-t0:.0f}s: "
          + ", ".join(f"L{L} fraction {np.mean(gt[L] > 0.5):.3f}" for L in LEVELS), flush=True)
    res = {"budget": a.budget, "n_test": a.n_test, "gt_seeds": a.gt_seeds, "ref": REF, "names": pc.NAMES,
           "levels": {}}
    nest_gt = {"L2_not_L1": float(np.mean((gt[2] > 0.5) & (gt[1] <= 0.5))),
               "L3_not_L2": float(np.mean((gt[3] > 0.5) & (gt[2] <= 0.5))),
               "L4_not_L3": float(np.mean((gt[4] > 0.5) & (gt[3] <= 0.5)))}
    res["nesting_ground_truth"] = nest_gt
    res["life_fraction_gt"] = {f"L{L}": float(np.mean(gt[L] > 0.5)) for L in LEVELS}
    res["share_of_L3_reaching_L4"] = float(np.mean(gt[4][gt[3] > 0.5] > 0.5))
    low_lk = test_x[:, pc.NAMES.index("lk")] < 0.1          # lk < 0.001
    res["life_fraction_gt_lk_below_0.001"] = {f"L{L}": float(np.mean(gt[L][low_lk] > 0.5)) for L in LEVELS}
    res["n_lk_below_0.001"] = int(low_lk.sum())
    maps = {}

    for L in MAP_LEVELS:
        pc.LEVEL = L
        p_true = gt[L]
        hard = ((p_true >= 0.1) & (p_true <= 0.4)) | ((p_true >= 0.6) & (p_true <= 0.9))
        lvl = {"life_fraction_gt": float(np.mean(p_true > 0.5)), "hard": int(hard.sum())}
        for name, fn in (("engine4_soft", lambda sim: bench1b.m_engine3(sim, a.budget, 0, dev, soft=True)),
                         ("sobol_1seed", lambda sim: bench1b.m_points(sim, a.budget, 0, 1, True))):
            t = time.time()
            x, s, n, extra = fn(bench1b.Sim(dev, 2000 + 10 * L, 1 << 15))
            curve = bench1b.evaluate(x, s, n, test_x, p_true, hard, dev, a.budget, 0)
            lvl[name] = {"curve": curve, "seconds": time.time() - t, **extra}
            print(f"L{L} {name}: {curve[-1]} ({time.time()-t:.0f}s)", flush=True)
            if name == "engine4_soft":
                pts = (x, s, n)
        # final COMPASS map on all points
        x, s, n = pts
        net = surrogate.fit_es_multi(torch.tensor(x, device=dev), torch.tensor(s, device=dev),
                                     torch.tensor(n, device=dev), [len(x)], seed=99, hidden=256, max_steps=10000)
        maps[L] = net
        p_map = predict(net, test_x, dev)
        lab = p_true > 0.5
        acc = float(np.mean((p_map > 0.5) == lab))
        lvl["final_map"] = {"accuracy": acc, "hard_accuracy": float(np.mean(((p_map > 0.5) == lab)[hard])),
                            "life_fraction_map": float(np.mean(p_map > 0.5))}
        # (b) permutation importance
        rng = np.random.default_rng(L)
        imp = {}
        for j, nm in enumerate(pc.NAMES):
            drops = []
            for _ in range(3):
                xp = test_x.copy()
                xp[:, j] = xp[rng.permutation(len(xp)), j]
                drops.append(acc - float(np.mean((predict(net, xp, dev) > 0.5) == lab)))
            imp[nm] = float(np.mean(drops))
        lvl["importance"] = imp
        # (c) partial dependence over random backgrounds
        bg = np.random.default_rng(100 + L).random((1000, pc.DIM)).astype(np.float32)
        grid = np.linspace(0, 1, 25, dtype=np.float32)
        pd_ = {}
        for j, nm in enumerate(pc.NAMES):
            xs = np.repeat(bg[None], len(grid), 0)
            xs[:, :, j] = grid[:, None]
            p = predict(net, xs.reshape(-1, pc.DIM), dev).reshape(len(grid), -1).mean(1)
            pd_[nm] = {"u": grid.tolist(), "value": [float(FROM_U[nm](q)) for q in grid], "p": p.tolist()}
        lvl["partial_dependence"] = pd_
        # 2D slices at the reference values
        g2 = np.linspace(0, 1, 41, dtype=np.float32)
        sl = {}
        for ax_x, ax_y in SLICES:
            jx, jy = pc.NAMES.index(ax_x), pc.NAMES.index(ax_y)
            u = np.repeat(ref_u()[None], len(g2) ** 2, 0)
            uu, vv = np.meshgrid(g2, g2)
            u[:, jx], u[:, jy] = uu.ravel(), vv.ravel()
            p = predict(net, u, dev).reshape(len(g2), len(g2))
            sl[f"{ax_x}_x_{ax_y}"] = {"x": ax_x, "y": ax_y, "x_values": [float(FROM_U[ax_x](q)) for q in g2],
                                      "y_values": [float(FROM_U[ax_y](q)) for q in g2], "p": p.tolist()}
        lvl["slices"] = sl
        res["levels"][f"L{L}"] = lvl
        print(f"L{L}: final map accuracy {acc:.3f}; importance "
              + ", ".join(f"{k} {v:.3f}" for k, v in sorted(imp.items(), key=lambda kv: -kv[1])), flush=True)

    # (e) read-offs: S_max by sigma at lk = 0 from the L3 and L4 maps; lk_crit by direct sweeps
    Sg = np.exp(np.linspace(math.log(8), math.log(512), 61))
    known = {}
    for L in MAP_LEVELS:
        known[f"L{L}"] = {"S_max_by_sigma_lk0": {f"{sg:.1f}": threshold_along(maps[L], "S", Sg, dev, rising=False, sigma=sg)[0]
                                                  for sg in (1.0, 1.2, 1.5, 2.0, 2.5, 3.0)}}
    lks = np.linspace(0, 0.01, 21)
    g = torch.Generator(device=dev).manual_seed(4242)
    sweep = {}
    for lv in lks:
        u = torch.tensor(np.repeat(ref_u(lk=float(lv))[None], 400, 0), device=dev)
        _, st = pc.simulate(u, generator=g, return_stats=True)
        sweep[f"{lv:.4f}"] = {f"L{L}": float(st[f"l{L}"].float().mean()) for L in LEVELS}
    for L in (2, 4):
        ps = [sweep[f"{lv:.4f}"][f"L{L}"] for lv in lks]
        crit = next((float(lks[i]) for i in range(1, len(ps)) if ps[i - 1] >= 0.5 > ps[i]), None)
        known.setdefault(f"L{L}", {})["lk_crit_direct"] = crit
    res["lk_sweep_reference_400pops"] = sweep
    res["known_behaviour"] = known
    res["seconds"] = time.time() - t0
    json.dump(res, open(a.out, "w"))
    print(f"\nlife fraction (ground truth): {res['life_fraction_gt']}")
    print(f"life fraction with lk < 0.001 (n={res['n_lk_below_0.001']}): {res['life_fraction_gt_lk_below_0.001']}")
    print(f"nesting violations (ground truth): {nest_gt}")
    print(f"share of L3 conditions that reach L4: {res['share_of_L3_reaching_L4']:.3f}")
    for L in MAP_LEVELS:
        v = res["levels"][f"L{L}"]
        e, s_ = v["engine4_soft"]["curve"][-1], v["sobol_1seed"]["curve"][-1]
        print(f"L{L} map: COMPASS {a.budget // 1000}k acc {e[1]:.3f} / hard {e[2]:.3f}; Sobol {s_[1]:.3f} / {s_[2]:.3f}; "
              f"hard set {v['hard']}; final map acc {v['final_map']['accuracy']:.3f}")
    print("lk sweep at the reference (400 pops; binomial sd <= 0.025):")
    for k, v in sweep.items():
        print(f"  lk={k}: " + " ".join(f"{q:.2f}" for q in v.values()))
    print("read-offs:", json.dumps(known, indent=1))
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
