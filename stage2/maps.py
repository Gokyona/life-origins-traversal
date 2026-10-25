"""Stage-2e: condition maps for life levels L1, L2, L3 on the protocell model.

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
sys.path.insert(1, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "stage1b")))
import bench1b  # noqa: E402
import protocell as pc  # noqa: E402
import surrogate  # noqa: E402

bench1b.model = pc
bench1b.FAST = True
EVAL_FIT = None
REF = dict(S=64, sigma=1.5, mu=0.01, K=2, d=0.005, food=1.0, km=2.0, dm=0.005, z1=0.5, z2=0.5)
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
    "z1": lambda v: v,
    "z2": lambda v: v,
}
FROM_U = {
    "S": lambda u: 8 * 64 ** u, "sigma": lambda u: 1 + 2 * u, "mu": lambda u: 1e-3 * 300 ** u,
    "K": lambda u: np.minimum(1 + np.floor(u * 3), 3), "d": lambda u: 0.02 * u, "food": lambda u: 0.25 * 4 ** u,
    "km": lambda u: 0.1 * 300 ** u, "dm": lambda u: 0.02 * u, "z1": lambda u: u, "z2": lambda u: u,
}
SLICES = [("S", "sigma"), ("mu", "S"), ("km", "food")]


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
    ap.add_argument("--out", default="maps.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    t0 = time.time()

    # ---- one ground truth for all three levels ----
    g = torch.Generator(device=dev).manual_seed(12345)
    test_x = np.random.default_rng(12345).random((a.n_test, pc.DIM)).astype(np.float32)
    xr = torch.tensor(np.repeat(test_x, a.gt_seeds, 0), device=dev)
    stats = {k: [] for k in ("l1", "l2", "l3")}
    for i in range(0, len(xr), 1 << 15):
        _, st = pc.simulate(xr[i:i + (1 << 15)], generator=g, return_stats=True)
        for k in stats:
            stats[k].append(st[k].float().cpu())
    gt = {L: torch.cat(stats[f"l{L}"]).view(a.n_test, a.gt_seeds).mean(1).numpy() for L in (1, 2, 3)}
    print(f"ground truth {a.n_test} x {a.gt_seeds} in {time.time()-t0:.0f}s: "
          + ", ".join(f"L{L} fraction {np.mean(gt[L] > 0.5):.3f}" for L in (1, 2, 3)), flush=True)
    res = {"budget": a.budget, "n_test": a.n_test, "gt_seeds": a.gt_seeds, "ref": REF, "names": pc.NAMES,
           "levels": {}}
    nest_gt = {"L2_not_L1": float(np.mean((gt[2] > 0.5) & (gt[1] <= 0.5))),
               "L3_not_L2": float(np.mean((gt[3] > 0.5) & (gt[2] <= 0.5)))}
    res["nesting_ground_truth"] = nest_gt
    maps = {}

    for L in (1, 2, 3):
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

    # (d) nesting in the maps
    p1, p2, p3 = (predict(maps[L], test_x, dev) > 0.5 for L in (1, 2, 3))
    res["nesting_maps"] = {"L2_not_L1": float(np.mean(p2 & ~p1)), "L3_not_L2": float(np.mean(p3 & ~p2))}

    # (e) known behaviour read off the maps (reference values elsewhere)
    Sg = np.exp(np.linspace(math.log(8), math.log(512), 61))
    known = {}
    for L in (2, 3):
        net = maps[L]
        k = {}
        k["S_max_by_sigma"] = {f"{sg:.1f}": threshold_along(net, "S", Sg, dev, rising=False, sigma=sg)[0]
                               for sg in (1.0, 1.2, 1.5, 2.0, 2.5, 3.0)}
        k["mu_critical"] = threshold_along(net, "mu", 1e-3 * 300 ** np.linspace(0, 1, 61), dev, rising=False)[0]
        kmg = 0.1 * 300 ** np.linspace(0, 1, 61)
        k["km_window"] = [threshold_along(net, "km", kmg, dev, rising=True)[0],
                          threshold_along(net, "km", kmg, dev, rising=False)[0]]
        k["K3_S_min"] = threshold_along(net, "S", Sg, dev, rising=True, K=3, sigma=1.0, mu=0.001, km=3.0)[0]
        k["food_min"] = threshold_along(net, "food", 0.25 * 4 ** np.linspace(0, 1, 61), dev, rising=True)[0]
        known[f"L{L}"] = k
    res["known_behaviour"] = known
    res["seconds"] = time.time() - t0
    json.dump(res, open(a.out, "w"))

    print(f"\nnesting violations (share of test points): ground truth {nest_gt}, maps {res['nesting_maps']}")
    print("\n| level | life fraction (truth / map) | COMPASS acc all / hard | Sobol acc all / hard | hard set |")
    print("|---|---|---|---|---|")
    for L in (1, 2, 3):
        v = res["levels"][f"L{L}"]
        e, s_ = v["engine4_soft"]["curve"][-1], v["sobol_1seed"]["curve"][-1]
        print(f"| L{L} | {v['life_fraction_gt']:.3f} / {v['final_map']['life_fraction_map']:.3f} | "
              f"{e[1]:.3f} / {e[2]:.3f} | {s_[1]:.3f} / {s_[2]:.3f} | {v['hard']} |")
    print("\nknown behaviour read off the maps:", json.dumps(known, indent=1))
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
