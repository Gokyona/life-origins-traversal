"""How well does the cheap protocell simulator (fidelity 0) agree with the reference (fidelity 1)?

On --n random conditions, estimate P(L) at both fidelities from --seeds runs each and
report: cost per simulation, label agreement (p > 0.5) with a confusion matrix,
correlation of the probabilities, and agreement split by how clear-cut the
high-fidelity probability is.
"""
import argparse, os, sys, time
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import protocell as pc  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--n", type=int, default=4000)
    ap.add_argument("--seeds", type=int, default=16)
    ap.add_argument("--level", type=int, default=3)
    a = ap.parse_args()
    dev = torch.device(a.device)
    sync = (lambda: torch.xpu.synchronize()) if dev.type == "xpu" else (lambda: None)
    g = torch.Generator(device=dev).manual_seed(11)
    x = torch.rand(a.n, pc.DIM, device=dev, generator=g)
    xr = x.repeat_interleave(a.seeds, 0)
    p, cost = {}, {}
    for f in (0, 1):
        pc.simulate(xr[:256], generator=g, fidelity=f)
        sync()
        t = time.time()
        y = torch.cat([pc.simulate(xr[i:i + 16384], generator=g, level=a.level, fidelity=f)
                       for i in range(0, len(xr), 16384)])
        sync()
        cost[f] = (time.time() - t) / len(xr)
        p[f] = y.float().view(a.n, a.seeds).mean(1).cpu().numpy()
        print(f"fidelity {f} {pc.FIDELITY[f]}: {1000 * cost[f]:.3f} ms per simulation; "
              f"P(L{a.level}) over the space {np.mean(p[f] > 0.5):.3f}", flush=True)
    lo, hi = p[0] > 0.5, p[1] > 0.5
    print(f"\ncost ratio high/low: {cost[1] / cost[0]:.1f}x")
    print(f"label agreement (p>0.5): {np.mean(lo == hi):.3f}")
    print(f"confusion (rows high-fi, cols low-fi): no/no {np.sum(~hi & ~lo)}, no/yes {np.sum(~hi & lo)}, "
          f"yes/no {np.sum(hi & ~lo)}, yes/yes {np.sum(hi & lo)}")
    print(f"correlation of probabilities: {np.corrcoef(p[0], p[1])[0, 1]:.3f}")
    clear = (p[1] <= 0.1) | (p[1] >= 0.9)
    print(f"agreement where high-fi is clear-cut (p<=0.1 or >=0.9, {clear.mean():.0%} of points): "
          f"{np.mean((lo == hi)[clear]):.3f}; elsewhere: {np.mean((lo == hi)[~clear]):.3f}")
    # which way does low-fi err? per condition, mean low-fi minus high-fi probability by high-fi bin
    bins = [0, 0.1, 0.4, 0.6, 0.9, 1.01]
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (p[1] >= b0) & (p[1] < b1)
        if m.any():
            print(f"  high-fi p in [{b0:.1f},{min(b1, 1):.1f}): n={m.sum():4d}, mean low-fi p {p[0][m].mean():.2f}, "
                  f"mean high-fi p {p[1][m].mean():.2f}")


if __name__ == "__main__":
    main()
