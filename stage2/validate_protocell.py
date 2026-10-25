"""Sanity sweeps for protocell.py (stage 2b).

  (a) reference setting: membrane-triggered division works and L3 holds
  (b) km sweep: failure at low km (crowding, slow division, bursting), at high km
      (under-filled vesicles, homeostasis lost), success in between
  (c) S sweep at sigma = 2 and K sweep at small S: the 2a trends survive coupling
  (d) mu sweep: a critical mutation rate is still present
Per point over --pops populations: P(L1), P(L2), P(L3), functional share,
generations per slot, content/capacity ratio, lysis per division (second half).
"""
import argparse, math, os, sys, time
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import protocell as pc  # noqa: E402

REF = dict(S=64, sigma=1.5, mu=0.01, K=2, d=0.005, food=1.0, km=2.0, dm=0.005)


def to_u(S, sigma, mu, K, d, food, km, dm, z=0.5):
    return [math.log(S / 8) / math.log(64), (sigma - 1) / 2, math.log10(mu / 1e-3) / math.log10(300),
            (K - 0.5) / 3, d / 0.02, math.log(food / 0.25) / math.log(4), (math.log10(km) + 1) / math.log10(300),
            dm / 0.02, z, z]


def sweep(title, key, values, pops, dev, g, **over):
    print(f"\n{title}")
    print(f"{key:>7} | P(L1) P(L2) P(L3) | funct  gens  ratio  lysis/div")
    for v in values:
        kw = dict(REF, **over, **{key: v})
        u = torch.tensor([to_u(**kw)] * pops, device=dev, dtype=torch.float32)
        _, st = pc.simulate(u, generator=g, return_stats=True)
        m = lambda k: st[k].float().mean().item()
        print(f"{v:>7} | {m('l1'):.2f}  {m('l2'):.2f}  {m('l3'):.2f}  | {m('frac_alive'):.2f}  "
              f"{m('gens'):5.1f}  {m('ratio'):.2f}   {m('lysis_frac'):.3f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--pops", type=int, default=200)
    a = ap.parse_args()
    dev = torch.device(a.device)
    g = torch.Generator(device=dev).manual_seed(0)
    sync = (lambda: torch.xpu.synchronize()) if dev.type == "xpu" else (lambda: None)

    u = torch.rand(4096, pc.DIM, device=dev, generator=g)
    pc.simulate(u[:64], generator=g)
    sync()
    t = time.time()
    _, st = pc.simulate(u, generator=g, return_stats=True)
    sync()
    dt = time.time() - t
    print(f"cost: {len(u)} random populations in {dt:.1f}s = {1000 * dt / len(u):.2f} ms each")
    print("over the space: " + ", ".join(f"P({k.upper()}) {st[k].float().mean():.3f}" for k in ("l1", "l2", "l3")))

    sweep("(a)+(b) km sweep at the reference setting " + str(REF), "km",
          [0.1, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0], a.pops, dev, g)
    Ss = [8, 16, 32, 64, 128, 256, 512]
    sweep("(c) S sweep, sigma=2", "S", Ss, a.pops, dev, g, sigma=2.0)
    for K in (1, 2, 3):
        sweep(f"(c) S sweep, K={K}, sigma=1, mu=0.001, km={K}", "S", Ss, a.pops, dev, g,
              sigma=1.0, mu=0.001, K=K, km=float(K))
    sweep("(d) mu sweep", "mu", [0.001, 0.003, 0.01, 0.03, 0.1, 0.3], a.pops, dev, g)


if __name__ == "__main__":
    main()
