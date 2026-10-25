"""Sanity sweeps for scm.py against the stochastic corrector's known behaviour.

  (a) parasite advantage sigma > 1: survival falls as the division size S grows
  (b) K >= 2 essential genes, no parasite advantage: very small S loses gene
      types at division (assortment load), K = 1 does not
  (c) a critical mutation rate mu above which populations collapse
P(L1), P(L2) per point over --pops independent populations; plus cost per simulation.
"""
import argparse, math, os, sys, time
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scm  # noqa: E402

def to_u(S=64, sigma=1.5, mu=0.01, K=2, d=0.005, food=1.0, z=0.5):
    return [math.log(S / 8) / math.log(64), (sigma - 1) / 2, math.log10(mu / 1e-3) / math.log10(300),
            (K - 0.5) / 3, d / 0.02, math.log(food / 0.25) / math.log(4), z, z]


def sweep(title, key, values, pops, dev, g, **fixed):
    print(f"\n{title}")
    print(f"{key:>8} | P(L1)  P(L2) | functional  gens")
    for v in values:
        kw = dict(fixed, **{key: v})
        u = torch.tensor([to_u(**kw)] * pops, device=dev, dtype=torch.float32)
        _, st = scm.simulate(u, generator=g, return_stats=True)
        print(f"{v:>8} | {st['l1'].float().mean():.2f}   {st['l2'].float().mean():.2f}  |"
              f"    {st['frac_alive'].mean():.2f}    {st['gens'].mean():5.1f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--pops", type=int, default=200)
    a = ap.parse_args()
    dev = torch.device(a.device)
    g = torch.Generator(device=dev).manual_seed(0)

    u = torch.rand(4096, scm.DIM, device=dev, generator=g)
    scm.simulate(u[:64], generator=g)
    torch.xpu.synchronize() if dev.type == "xpu" else None
    t = time.time()
    y = scm.simulate(u, generator=g)
    torch.xpu.synchronize() if dev.type == "xpu" else None
    dt = time.time() - t
    print(f"cost: {len(u)} random populations in {dt:.1f}s = {1000 * dt / len(u):.2f} ms each; "
          f"P(L{scm.LEVEL}) over the space {y.float().mean():.3f}")
    _, st = scm.simulate(u, generator=g, return_stats=True)
    print(f"over the space: P(L1) {st['l1'].float().mean():.3f}, P(L2) {st['l2'].float().mean():.3f}")

    Ss = [8, 16, 32, 64, 128, 256, 512]
    sweep("(a) S sweep, sigma=2, mu=0.01, K=2, d=0.005, food=1", "S", Ss, a.pops, dev, g,
          sigma=2.0, mu=0.01, K=2)
    sweep("(a') S sweep, sigma=1.2", "S", Ss, a.pops, dev, g, sigma=1.2, mu=0.01, K=2)
    for K in (1, 2, 3):
        sweep(f"(b) S sweep, K={K}, sigma=1, mu=0.001", "S", Ss, a.pops, dev, g, sigma=1.0, mu=0.001, K=K)
    sweep("(d) food sweep, S=64, sigma=1.5, mu=0.01, K=2", "food", [0.25, 0.35, 0.5, 0.7, 1.0], a.pops, dev, g,
          S=64, sigma=1.5, mu=0.01, K=2)
    sweep("(c) mu sweep, S=64, sigma=1.5, K=2", "mu", [0.001, 0.003, 0.01, 0.03, 0.1, 0.3], a.pops, dev, g,
          S=64, sigma=1.5, K=2)


if __name__ == "__main__":
    main()
