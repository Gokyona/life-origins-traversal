"""Check of the 3b boxes: P(life) inside each 95% box restricted to the corner S > 256 and sigma > 2.

Regenerates the confirmation points of s3b_box.py (same seeds) and simulates them again
(same generator seeds, so the same runs), then reports the corner subset.
"""
import json, math, os, sys
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "stage2"))
import knockouts as ko  # noqa: E402

dev = torch.device("xpu")
r = json.load(open(os.path.join(HERE, "s3b.json")))
S_u = math.log(256 / 8) / math.log(64)
sig_u = (2 - 1) / 2
for L in (1, 2, 3):
    box = r["levels"][f"L{L}"]["target_95"]
    lo, hi = np.zeros(ko.DIM), np.ones(ko.DIM)
    for st in box["steps"]:
        j = ko.NAMES.index(st["condition"])
        lo[j], hi[j] = st["u_range"]
    U = lo + np.random.default_rng(600 + L).random((2000, ko.DIM)) * (hi - lo)
    g = torch.Generator(device=dev).manual_seed(500 + L)
    ur = torch.tensor(np.repeat(U, 16, 0), dtype=torch.float32, device=dev)
    y = torch.cat([ko.simulate(ur[i:i + (1 << 15)], generator=g, level=L)
                   for i in range(0, len(ur), 1 << 15)]).float().view(2000, 16).mean(1).cpu().numpy()
    m = (U[:, 0] > S_u) & (U[:, 1] > sig_u)
    print(f"L{L}: whole box share P>0.5 {np.mean(y > 0.5):.3f}; corner S>256 & sigma>2: n={m.sum()}"
          + (f", share P>0.5 {np.mean(y[m] > 0.5):.3f}, mean P {y[m].mean():.3f}" if m.any() else " (box excludes it)"))
