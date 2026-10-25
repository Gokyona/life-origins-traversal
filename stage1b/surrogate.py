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


def fit_es(x, s, n, seed=0, max_steps=10000, batch=4096, hidden=256, lr=3e-3, bootstrap=False,
           val_frac=0.1, eval_every=200, patience=5):
    """Fit with early stopping on a held-out split, so training length adapts to data size.

    A fixed number of steps underfits large datasets and overfits small ones
    (stage 1c diagnosis); here the network trains until held-out loss stops improving.
    """
    g = torch.Generator(device="cpu").manual_seed(seed)
    torch.manual_seed(seed)
    N = x.shape[0]
    perm = torch.randperm(N, generator=g)
    n_val = max(1, int(N * val_frac))
    val, tr = perm[:n_val].to(x.device), perm[n_val:]
    if bootstrap:
        tr = tr[torch.randint(0, len(tr), (len(tr),), generator=g)]
    net = _net(x.shape[1], hidden).to(x.device)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    target = s / n
    best, best_state, bad = float("inf"), None, 0
    for step in range(1, max_steps + 1):
        idx = tr[torch.randint(0, len(tr), (min(batch, len(tr)),), generator=g)].to(x.device)
        loss = (nn.functional.binary_cross_entropy_with_logits(net(x[idx]).squeeze(1), target[idx],
                                                               reduction="none") * n[idx]).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
        if step % eval_every == 0:
            with torch.no_grad():
                vl = (nn.functional.binary_cross_entropy_with_logits(net(x[val]).squeeze(1), target[val],
                                                                     reduction="none") * n[val]).mean().item()
            if vl < best - 1e-5:
                best, bad = vl, 0
                best_state = {k: v.clone() for k, v in net.state_dict().items()}
            else:
                bad += 1
                if bad >= patience:
                    break
    if best_state is not None:
        net.load_state_dict(best_state)
    net.eval()
    return net


class Ensemble(nn.Module):
    """K independent 2-hidden-layer networks evaluated as one batched matrix product."""

    def __init__(self, k, dim, hidden, device):
        super().__init__()

        def lin(i, o):
            b = 1.0 / i ** 0.5
            return (nn.Parameter(torch.empty(k, i, o, device=device).uniform_(-b, b)),
                    nn.Parameter(torch.empty(k, 1, o, device=device).uniform_(-b, b)))

        self.w1, self.b1 = lin(dim, hidden)
        self.w2, self.b2 = lin(hidden, hidden)
        self.w3, self.b3 = lin(hidden, 1)

    def forward(self, x):
        """x [K,B,D] (one batch per member) or [B,D] (same batch for all) -> logits [K,B]."""
        h = nn.functional.silu(torch.matmul(x, self.w1) + self.b1)
        h = nn.functional.silu(torch.bmm(h, self.w2) + self.b2)
        return (torch.bmm(h, self.w3) + self.b3).squeeze(2)


@torch.no_grad()
def predict_multi(ens, q, chunk=1 << 15):
    """Probabilities [K,Q] for every member on query points q [Q,D]."""
    return torch.cat([torch.sigmoid(ens(q[i:i + chunk])) for i in range(0, q.shape[0], chunk)], 1)


def fit_es_multi(x, s, n, sizes, seed=0, max_steps=10000, batch=4096, hidden=256, lr=3e-3,
                 bootstrap=False, val_frac=0.1, eval_every=200, patience=5):
    """Train len(sizes) networks at once; member i sees only the first sizes[i] points.

    Same early-stopping rule as fit_es, applied per member (each keeps its own best
    weights); training ends when every member has stopped improving. Used both for
    the engine's bootstrap ensemble (equal sizes, bootstrap=True) and for scoring all
    checkpoints of one run together (growing sizes). Adam is elementwise, so one
    optimizer over the stacked weights equals one optimizer per member.
    """
    dev = x.device
    g = torch.Generator(device="cpu").manual_seed(seed)
    torch.manual_seed(seed)
    K = len(sizes)
    tr_sets, val_sets = [], []
    for k in sizes:
        perm = torch.randperm(k, generator=g)
        n_val = max(1, int(k * val_frac))
        tr = perm[n_val:]
        if bootstrap:
            tr = tr[torch.randint(0, len(tr), (len(tr),), generator=g)]
        tr_sets.append(tr)
        val_sets.append(perm[:n_val])
    pad = lambda sets: torch.stack([torch.cat([t, t[:1].repeat(max(map(len, sets)) - len(t))]) for t in sets]).to(dev)
    tr_pad, val_pad = pad(tr_sets), pad(val_sets)
    tr_len = torch.tensor([len(t) for t in tr_sets], device=dev, dtype=torch.float32)
    val_len = torch.tensor([len(t) for t in val_sets], device=dev)
    val_mask = (torch.arange(val_pad.shape[1], device=dev)[None, :] < val_len[:, None]).float()
    target = s / n
    ens = Ensemble(K, x.shape[1], hidden, dev)
    opt = torch.optim.Adam(ens.parameters(), lr=lr)
    gd = torch.Generator(device=dev).manual_seed(seed)
    B = min(batch, int(tr_len.max().item()))
    best = torch.full((K,), float("inf"), device=dev)
    best_state = [p.detach().clone() for p in ens.parameters()]
    bad = torch.zeros(K, device=dev)
    for step in range(1, max_steps + 1):
        pos = (torch.rand(K, B, device=dev, generator=gd) * tr_len[:, None]).long()
        idx = torch.gather(tr_pad, 1, pos)
        loss = (nn.functional.binary_cross_entropy_with_logits(ens(x[idx]), target[idx], reduction="none")
                * n[idx]).mean(1).sum()
        opt.zero_grad()
        loss.backward()
        opt.step()
        if step % eval_every == 0:
            with torch.no_grad():
                vl = nn.functional.binary_cross_entropy_with_logits(ens(x[val_pad]), target[val_pad],
                                                                    reduction="none") * n[val_pad]
                vl = (vl * val_mask).sum(1) / val_len
                better = vl < best - 1e-5
                best = torch.where(better, vl, best)
                bad = torch.where(better, torch.zeros_like(bad), bad + 1)
                for b, p in zip(best_state, ens.parameters()):
                    b[better] = p[better]
                if bool((bad >= patience).all()):
                    break
    with torch.no_grad():
        for b, p in zip(best_state, ens.parameters()):
            p.copy_(b)
    ens.eval()
    return ens
