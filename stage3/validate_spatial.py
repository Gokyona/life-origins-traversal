"""Validation of spatial.py (stage 3d) against the well-mixed evolve.py.

At the reference (2b reference, eta = 0.1, lk = 0) and along the key lines from
2a/2b/3c, P(L1..L4) for the spatial and the well-mixed model side by side; the
equal_split and no_selection knockouts in both; a leak sweep; cost per simulation
and P(L1..L4) over the whole space.
"""
import argparse, os, sys, time
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "stage2"))
import evolve as ev  # noqa: E402
import spatial as sp  # noqa: E402
from validate_knockouts import to_u as to_u10, REF  # noqa: E402


def u_for(model, eta=0.1, lk=0.0, **kw):
    u = to_u10(**dict(REF, **kw))
    if model is sp:
        return u[:8] + [eta / sp.ETA_MAX, lk / sp.LK_MAX] + u[8:]
    return u[:8] + [eta / ev.ETA_MAX] + u[8:]


def levels(model, pops, dev, g, knockout=None, **kw):
    u = torch.tensor([u_for(model, **kw)] * pops, dtype=torch.float32, device=dev)
    _, st = model.simulate(u, generator=g, return_stats=True, knockout=knockout)
    return [st[f"l{L}"].float().mean().item() for L in (1, 2, 3, 4)], st["mean_q"].mean().item()


def row(label, pops, dev, g, knockout=None, models=(sp, ev), **kw):
    cells = []
    for m in models:
        p, q = levels(m, pops, dev, g, knockout=knockout, **kw)
        cells.append(" ".join(f"{x:.2f}" for x in p) + f" q {q:.2f}")
    print(f"  {label:>28} | " + " | ".join(cells), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--pops", type=int, default=200)
    a = ap.parse_args()
    dev = torch.device(a.device)
    sync = (lambda: torch.xpu.synchronize()) if dev.type == "xpu" else (lambda: None)
    g = torch.Generator(device=dev).manual_seed(0)
    print(f"{'':>30} | spatial P(L1 L2 L3 L4), mean q | well-mixed P(L1 L2 L3 L4), mean q")
    print("reference")
    row("reference", a.pops, dev, g)
    print("S line, sigma = 2")
    for S in (16, 64, 128, 256, 512):
        row(f"S={S}", a.pops, dev, g, S=S, sigma=2.0)
    print("mu sweep")
    for mu in (0.003, 0.01, 0.03, 0.1):
        row(f"mu={mu}", a.pops, dev, g, mu=mu)
    print("km sweep")
    for km in (0.1, 0.3, 2.0, 10.0):
        row(f"km={km}", a.pops, dev, g, km=km)
    print("knockouts")
    for k in ("equal_split", "no_selection"):
        row(f"{k} ref", a.pops, dev, g, knockout=k)
        row(f"{k} S=256 sigma=2", a.pops, dev, g, knockout=k, S=256, sigma=2.0)
    row("baseline S=256 sigma=2", a.pops, dev, g, S=256, sigma=2.0)
    row("K=3 S=8 sigma=1 mu=.001", a.pops, dev, g, K=3, S=8, sigma=1.0, mu=0.001, km=3.0)
    print("leak (spatial only)")
    for lk in (0.0, 0.001, 0.003, 0.01):
        row(f"lk={lk}", a.pops, dev, g, lk=lk, models=(sp,))
        row(f"lk={lk} S=256 sigma=2", a.pops, dev, g, lk=lk, S=256, sigma=2.0, models=(sp,))

    print("\nover the whole space (4096 random conditions each)")
    for m in (sp, ev):
        x = torch.rand(4096, m.DIM, device=dev, generator=g)
        m.simulate(x[:64], generator=g)
        sync()
        t = time.time()
        _, st = m.simulate(x, generator=g, return_stats=True)
        sync()
        print(f"  {m.__name__:>8}: {1000 * (time.time() - t) / len(x):.2f} ms per simulation; "
              + ", ".join(f"P(L{L}) {st[f'l{L}'].float().mean():.3f}" for L in (1, 2, 3, 4)), flush=True)


if __name__ == "__main__":
    main()
