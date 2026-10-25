"""Shared map-maker: a small neural network fitted to (condition, successes, runs).

Every method's tested points are turned into a life/no-life map with the same
fit, so methods differ only in which conditions they chose to simulate.
An ensemble of bootstrap fits gives the engine an estimate of *epistemic*
uncertainty (where more data would change the map), as opposed to outcomes
that are simply random.
"""
import torch
import torch.nn as nn


def _net(dim, hidden):
    return nn.Sequential(nn.Linear(dim, hidden), nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU(),
                         nn.Linear(hidden, 1))


def fit(x, s, n, seed=0, steps=600, batch=4096, hidden=128, lr=3e-3, bootstrap=False):
    """x [N,D], s/n successes and runs per condition (tensors on the same device)."""
    g = torch.Generator(device="cpu").manual_seed(seed)
    torch.manual_seed(seed)
    net = _net(x.shape[1], hidden).to(x.device)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
    target = s / n
    N = x.shape[0]
    pool = torch.randint(0, N, (N,), generator=g) if bootstrap else torch.arange(N)
    for _ in range(steps):
        idx = pool[torch.randint(0, N, (min(batch, N),), generator=g)].to(x.device)
        logit = net(x[idx]).squeeze(1)
        loss = (nn.functional.binary_cross_entropy_with_logits(logit, target[idx], reduction="none") * n[idx]).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
        sched.step()
    net.eval()
    return net


@torch.no_grad()
def predict(net, q, chunk=1 << 16):
    return torch.cat([torch.sigmoid(net(q[i:i + chunk]).squeeze(1)) for i in range(0, q.shape[0], chunk)])
