"""Validation of evolve.py (stage 3c, level L4).

  * reference setting (2b reference, eta = 0.1): adaptive rise of mean q against the
    paired neutral marker q', with the time course of both
  * double-neutral control (q has no effect): spread of dq = mean(q) - mean(q') at the
    reference and over random conditions -> false-positive rate of L4 at L4_DELTA
  * eta sweep (eta = 0 must give no L4), knockouts (no_selection, equal_split), and
    a setting where L3 fails (no L4 there by construction)
  * cost per simulation and P(L1..L4) over the whole space
"""
import argparse, math, os, sys, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "stage2"))
import evolve as ev  # noqa: E402
from validate_knockouts import to_u as to_u10, REF  # noqa: E402


def to_u(eta=0.1, **kw):
    u = to_u10(**dict(REF, **kw))
    return u[:8] + [eta / ev.ETA_MAX] + u[8:]


def run(label, pops, dev, g, eta=0.1, knockout=None, neutral_q=False, **kw):
    u = torch.tensor([to_u(eta, **kw)] * pops, dtype=torch.float32, device=dev)
    _, st = ev.simulate(u, generator=g, return_stats=True, knockout=knockout, neutral_q=neutral_q)
    dq = st["dq"].cpu().numpy()
    m = lambda k: st[k].float().mean().item()
    print(f"  {label:>34} | L3 {m('l3'):.2f} L4 {m('l4'):.2f} | mean q {m('mean_q'):.3f} q' {m('mean_q_neutral'):.3f} "
          f"dq {dq.mean():+.3f} ± {dq.std():.3f} | z {st['z'].mean().item():+.1f} | gens {m('gens'):.0f}", flush=True)
    return dq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="xpu")
    ap.add_argument("--pops", type=int, default=200)
    a = ap.parse_args()
    dev = torch.device(a.device)
    sync = (lambda: torch.xpu.synchronize()) if dev.type == "xpu" else (lambda: None)
    g = torch.Generator(device=dev).manual_seed(0)
    print(f"L4_DELTA = {ev.L4_DELTA}, L4_Z = {ev.L4_Z}, M_NEUTRAL = {ev.M_NEUTRAL}, Q_MAX = {ev.Q_MAX}, h(q) = 2q/(1+q)")

    print("\n(1) reference, time course of mean q and mean q' (functional cells, mean over populations)")
    u = torch.tensor([to_u(0.1)] * a.pops, dtype=torch.float32, device=dev)
    _, st = ev.simulate(u, generator=g, return_stats=True, trace_every=100)
    for step, mq, mqn in st["trace"]:
        print(f"  step {step:4d}: q {mq:.3f}  q' {mqn:.3f}")
    print(f"  end: q {st['mean_q'].mean():.3f}  q' {st['mean_q_neutral'].mean():.3f}, P(L3) {st['l3'].float().mean():.2f}, "
          f"P(L4) {st['l4'].float().mean():.2f}")

    print("\n(2) eta sweep at the reference")
    for eta in (0.0, 0.02, 0.05, 0.1, 0.2):
        run(f"eta={eta}", a.pops, dev, g, eta=eta)

    print("\n(3) controls at the reference (eta=0.1)")
    dq_neutral = run("double-neutral (q has no effect)", a.pops, dev, g, neutral_q=True)
    run("knockout no_selection", a.pops, dev, g, knockout="no_selection")
    run("knockout equal_split", a.pops, dev, g, knockout="equal_split")
    run("L3 fails: S=512 sigma=2", a.pops, dev, g, S=512, sigma=2.0)
    run("L3 fails: km=0.1", a.pops, dev, g, km=0.1)

    print("\n(4) over the whole space (4096 random conditions)")
    x = torch.rand(4096, ev.DIM, device=dev, generator=g)
    ev.simulate(x[:64], generator=g)
    sync()
    t = time.time()
    _, st = ev.simulate(x, generator=g, return_stats=True)
    sync()
    print(f"  cost {1000 * (time.time() - t) / len(x):.2f} ms per simulation; "
          + ", ".join(f"P(L{L}) {st[f'l{L}'].float().mean():.3f}" for L in (1, 2, 3, 4)))
    _, stn = ev.simulate(x, generator=g, return_stats=True, neutral_q=True)
    dqn = stn["dq"][stn["l3"]].cpu().numpy()
    print(f"  double-neutral over the space: P(L4) {stn['l4'].float().mean():.3f} "
          f"(= false-positive rate among L3 populations {stn['l4'][stn['l3']].float().mean():.3f})")
    dqs = st["dq"][st["l3"]].cpu().numpy()
    print(f"  double-neutral over the space (L3 populations, n={len(dqn)}): dq {dqn.mean():+.3f} ± {dqn.std():.3f}, "
          f"99th pct {np.percentile(dqn, 99):.3f}")
    print(f"  with selection (L3 populations, n={len(dqs)}): dq {dqs.mean():+.3f} ± {dqs.std():.3f}, "
          f"P(L4 | L3) {st['l4'][st['l3']].float().mean():.3f}")
    print(f"  reference double-neutral: dq sd {dq_neutral.std():.3f}, 99th pct {np.percentile(dq_neutral, 99):.3f}")


if __name__ == "__main__":
    main()
