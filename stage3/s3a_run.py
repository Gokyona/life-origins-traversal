"""Stage-3a run: life regions of the protocell model with each mechanism knocked out.

For the baseline and every knockout in knockouts.KNOCKOUTS:
  * one ground-truth set (the same N_TEST conditions for all; GT_SEEDS runs each,
    simulated once with return_stats) gives P(L1), P(L2), P(L3) per condition;
  * life fraction per level, and the share of conditions gained / lost relative
    to the baseline (labels p > 0.5);
  * a COMPASS map per level (engine4_soft, BUDGET runs, --fast map-maker) scored
    against that ground truth (skipped where the level is undefined: L3 for
    no_compartments and no_coupling).
Writes s3a.json and prints the tables.
"""
import argparse, json, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "stage2"))
sys.path.insert(2, os.environ.get("LOT_STAGE1B", os.path.join(HERE, "..", "stage1b")))
import bench1b  # noqa: E402
import knockouts as ko  # noqa: E402

bench1b.model = ko
bench1b.FAST = True
UNDEFINED = {("no_compartments", 3), ("no_coupling", 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--budget", type=int, default=16000)
    ap.add_argument("--n-test", type=int, default=4000)
    ap.add_argument("--gt-seeds", type=int, default=64)
    ap.add_argument("--knockouts", default=",".join(str(k) for k in ko.KNOCKOUTS))
    ap.add_argument("--out", default="s3a.json")
    a = ap.parse_args()
    dev = torch.device(a.device)
    kos = [None if k == "None" else k for k in a.knockouts.split(",")]
    test_x = np.random.default_rng(12345).random((a.n_test, ko.DIM)).astype(np.float32)
    xr = torch.tensor(np.repeat(test_x, a.gt_seeds, 0), device=dev)
    res = {"budget": a.budget, "n_test": a.n_test, "gt_seeds": a.gt_seeds, "knockouts": {}}
    labels = {}
    t_all = time.time()
    for k in kos:
        name = str(k)
        ko.KNOCKOUT = k
        t = time.time()
        g = torch.Generator(device=dev).manual_seed(12345)
        acc = {f"l{L}": [] for L in (1, 2, 3)}
        for i in range(0, len(xr), 1 << 15):
            _, st = ko.simulate(xr[i:i + (1 << 15)], generator=g, return_stats=True)
            for key in acc:
                acc[key].append(st[key].float().cpu())
        p = {L: torch.cat(acc[f"l{L}"]).view(a.n_test, a.gt_seeds).mean(1).numpy() for L in (1, 2, 3)}
        labels[name] = {L: p[L] > 0.5 for L in (1, 2, 3)}
        entry = {"gt_seconds": time.time() - t, "levels": {}}
        print(f"\n== {name}: ground truth in {time.time()-t:.0f}s", flush=True)
        for L in (1, 2, 3):
            lab = labels[name][L]
            e = {"life_fraction": float(lab.mean())}
            if "None" in labels and name != "None":
                base = labels["None"][L]
                e["gained"] = float(np.mean(lab & ~base))
                e["lost"] = float(np.mean(~lab & base))
            if (k, L) not in UNDEFINED:
                ko.LEVEL = L
                hard = ((p[L] >= 0.1) & (p[L] <= 0.4)) | ((p[L] >= 0.6) & (p[L] <= 0.9))
                t2 = time.time()
                x, s, n, extra = bench1b.m_engine3(bench1b.Sim(dev, 2000 + L, 1 << 15), a.budget, 0, dev, soft=True)
                curve = bench1b.evaluate(x, s, n, test_x, p[L], hard, dev, a.budget, 0)
                e.update({"compass_curve": curve, "compass_seconds": time.time() - t2, "hard": int(hard.sum()), **extra})
                print(f"  L{L}: life {e['life_fraction']:.3f} gained {e.get('gained', 0):.3f} lost {e.get('lost', 0):.3f} | "
                      f"COMPASS {a.budget // 1000}k acc {curve[-1][1]:.3f} hard {curve[-1][2]:.3f} (n={hard.sum()}) "
                      f"({time.time()-t2:.0f}s)", flush=True)
            else:
                print(f"  L{L}: undefined for this knockout", flush=True)
            entry["levels"][f"L{L}"] = e
        res["knockouts"][name] = entry
        json.dump(res, open(a.out, "w"))

    print("\nLife fraction per level (ground truth, p > 0.5); gained / lost vs baseline in brackets")
    print("| knockout | L1 | L2 | L3 |")
    print("|---|---|---|---|")
    for name, entry in res["knockouts"].items():
        cells = []
        for L in (1, 2, 3):
            e = entry["levels"][f"L{L}"]
            if ("no_compartments" == name or "no_coupling" == name) and L == 3:
                cells.append("undefined")
            elif "gained" in e:
                cells.append(f"{e['life_fraction']:.3f} (+{e['gained']:.3f} / -{e['lost']:.3f})")
            else:
                cells.append(f"{e['life_fraction']:.3f}")
        print(f"| {name} | " + " | ".join(cells) + " |")
    print("\nCOMPASS map accuracy at the final budget (all / hard)")
    print("| knockout | L1 | L2 | L3 |")
    print("|---|---|---|---|")
    for name, entry in res["knockouts"].items():
        cells = []
        for L in (1, 2, 3):
            e = entry["levels"][f"L{L}"]
            cells.append(f"{e['compass_curve'][-1][1]:.3f} / {e['compass_curve'][-1][2]:.3f}" if "compass_curve" in e else "-")
        print(f"| {name} | " + " | ".join(cells) + " |")
    print(f"\ntotal {time.time()-t_all:.0f}s")


if __name__ == "__main__":
    main()
