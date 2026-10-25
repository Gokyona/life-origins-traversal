"""Short direct sweeps for each knockout (stage 3a): reference setting and the S line at sigma = 2.

Per point over --pops populations: P(L1), P(L2), P(L3), functional share, generations,
content/capacity ratio; plus cost per simulation for each knockout.
"""
import argparse, math, os, sys, time
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stage2"))
import knockouts as ko  # noqa: E402

REF = dict(S=64, sigma=1.5, mu=0.01, K=2, d=0.005, food=1.0, km=2.0, dm=0.005)


def to_u(S, sigma, mu, K, d, food, km, dm, z=0.5):
    return [math.log(S / 8) / math.log(64), (sigma - 1) / 2, math.log10(mu / 1e-3) / math.log10(300),
            (K - 0.5) / 3, d / 0.02, math.log(food / 0.25) / math.log(4), (math.log10(km) + 1) / math.log10(300),
            dm / 0.02, z, z]


def row(label, kw, k, pops, dev, g):
    u = torch.tensor([to_u(**kw)] * pops, device=dev, dtype=torch.float32)
    _, st = ko.simulate(u, generator=g, return_stats=True, knockout=k)
    m = lambda key: st[key].float().mean().item()
    print(f"  {label:>10} | {m('l1'):.2f}  {m('l2'):.2f}  {m('l3'):.2f}  | {m('frac_alive'):.2f}  {m('gens'):6.1f}  "
          f"{m('ratio'):.2f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--pops", type=int, default=200)
    a = ap.parse_args()
    dev = torch.device(a.device)
    sync = (lambda: torch.xpu.synchronize()) if dev.type == "xpu" else (lambda: None)
    g = torch.Generator(device=dev).manual_seed(0)
    for k in ko.KNOCKOUTS:
        u = torch.rand(2048, ko.DIM, device=dev, generator=g)
        ko.simulate(u[:64], generator=g, knockout=k)
        sync()
        t = time.time()
        ko.simulate(u, generator=g, knockout=k)
        sync()
        print(f"\n== knockout {k}: {1000 * (time.time() - t) / len(u):.2f} ms per simulation")
        print(f"  {'setting':>10} | P(L1) P(L2) P(L3) | funct  gens  ratio")
        row("reference", REF, k, a.pops, dev, g)
        for S in (8, 16, 64, 128, 256, 512):
            row(f"s2 S={S}", dict(REF, sigma=2.0, S=S), k, a.pops, dev, g)


if __name__ == "__main__":
    main()
