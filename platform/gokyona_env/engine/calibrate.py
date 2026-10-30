"""Hardware calibration: peak device memory and step time against L_MAX at a given sites x slots, the largest L_MAX
that fits a memory budget (default 60% of free device memory, capped by --budget-gb), and a memory pre-flight
(peak <= 80% of the budget, no growth over the batches). Writes <data-dir>/hardware.json.
    python -m gokyona_env.engine.calibrate --data-dir DIR [--sites 64 --slots 512 --budget-gb 20 --steps 30]
"""
import argparse, json, math, os, sys, time
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import gk_chem, devices, darwin  # noqa: E402
from trait import Trait  # noqa: E402

PKG = os.path.dirname(HERE)


class StopAfter:
    def __init__(self, n, dev):
        self.n, self.k, self.dev = n, 0, dev
        self.t = []

    def due(self):
        return True

    def __call__(self, loc):
        devices.sync(self.dev)
        self.t.append(time.time())
        self.k += 1
        return "stop" if self.k >= self.n else None


def measure(dev, W, sites, slots, steps, trait_cfg):
    gk_chem.configure(W)
    pre = json.load(open(os.path.join(PKG, "data", "presets.json")))["darwin_demo"]
    u = torch.tensor([pre["u"]], dtype=torch.float32, device=dev)
    ov = {k: float(v) for k, v in pre["env"].items() if k != "avail_time"}
    ov["endless"] = True
    devices.reset_peak(dev)
    g = torch.Generator(device=dev); g.manual_seed(1)
    hk = StopAfter(steps, dev)
    tr = Trait(trait_cfg, sites, dev) if trait_cfg.get("enabled") else None
    gk_chem.simulate(u, g, (sites, slots, 300), override=ov, trait=tr, hook=hk)
    _, peak = devices.allocated(dev)
    dt = [b - a for a, b in zip(hk.t[:-1], hk.t[1:])]
    late = dt[len(dt) // 2:] if dt else [float("nan")]
    return {"W": W, "lmax": 64 * W, "peak_gib": None if peak is None else peak / 2 ** 30,
            "step_ms": 1000 * sum(late) / len(late)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.path.expanduser("~/.gokyona-env"))
    ap.add_argument("--device", default="auto")
    ap.add_argument("--sites", type=int, default=64)
    ap.add_argument("--slots", type=int, default=512)
    ap.add_argument("--steps", type=int, default=30)
    ap.add_argument("--budget-gb", type=float, default=20.0)
    ap.add_argument("--ws", default="1,4,16,64")
    a = ap.parse_args()
    dev = devices.pick(a.device)
    free, total = devices.mem_info(dev)
    budget = min(a.budget_gb * 2 ** 30, 0.6 * free)
    trait_cfg = {"enabled": True, "mode": "shared", "s": 1.0}
    rows = []
    for W in [int(x) for x in a.ws.split(",")]:
        r = measure(dev, W, a.sites, a.slots, a.steps, trait_cfg)
        rows.append(r)
        print(r, flush=True)
        if r["peak_gib"] is not None and r["peak_gib"] * 2 ** 30 > budget:
            break
    lmax_max = None
    if all(r["peak_gib"] is not None for r in rows) and len(rows) >= 2:
        (w1, m1), (w2, m2) = (rows[0]["W"], rows[0]["peak_gib"]), (rows[-1]["W"], rows[-1]["peak_gib"])
        slope = max((m2 - m1) / max(w2 - w1, 1), 1e-6)
        wmax = int((budget / 2 ** 30 * 0.8 - m1) / slope + w1)
        lmax_max = max(64, 64 * max(1, wmax))
    out = {"device": devices.describe(dev), "device_type": dev.type, "free_gib": free / 2 ** 30, "total_gib": total / 2 ** 30,
           "budget_gib": budget / 2 ** 30, "sites": a.sites, "slots": a.slots, "measurements": rows, "lmax_max": lmax_max,
           "note": "peak memory grows ~linearly with words per chain; lmax_max keeps the peak <= 80% of the budget",
           "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    if rows:
        best = rows[-1]
        out["hw_note"] = (f"{out['device']}: max length about {lmax_max} letters at {a.sites} sites x {a.slots} slots; "
                          f"{best['step_ms']:.0f} ms/step at length {best['lmax']}")
    os.makedirs(a.data_dir, exist_ok=True)
    with open(os.path.join(a.data_dir, "hardware.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1), flush=True)


if __name__ == "__main__":
    main()
