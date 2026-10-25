"""Ordered condition-space testing engine (rules 1-6 of the paper).

Rule 1  canonical key + permanent archive: a condition is never simulated twice
Rule 2  Sobol low-discrepancy order for global coverage
Rule 3  monotone inference with spot checks
Rule 5  refine only where the surrogate is unsure (the life / no-life boundary)
Rule 6  sequential testing (SPRT-style) over enumerated seeds
(Rule 4, the multi-fidelity funnel, needs a second cheaper model and is left
for stage 2.)
"""
import numpy as np
import torch
from scipy.stats import beta, qmc

KEY_BITS = 12


def canonical_key(x):
    return tuple(np.floor(np.asarray(x) * (1 << KEY_BITS)).astype(np.int64).clip(0, (1 << KEY_BITS) - 1))


def knn_predict(train_x, train_p, train_w, query, k=10, chunk=4096):
    """Distance-weighted kNN estimate of success probability."""
    out = []
    for i in range(0, query.shape[0], chunk):
        q = query[i:i + chunk]
        d = torch.cdist(q, train_x)
        kk = min(k, train_x.shape[0])
        dist, idx = d.topk(kk, largest=False)
        w = train_w[idx] / (dist + 1e-3)
        out.append((w * train_p[idx]).sum(1) / w.sum(1))
    return torch.cat(out)


class Archive:
    def __init__(self):
        self.index = {}
        self.x, self.s, self.n, self.inferred = [], [], [], []

    def __contains__(self, x):
        return canonical_key(x) in self.index

    def add(self, x, s, n, inferred=False):
        key = canonical_key(x)
        if key in self.index:  # rule 1: never record a condition twice
            return False
        self.index[key] = len(self.x)
        self.x.append(np.asarray(x, dtype=np.float32))
        self.s.append(s)
        self.n.append(n)
        self.inferred.append(inferred)
        return True

    def tensors(self, device):
        x = torch.tensor(np.stack(self.x), device=device)
        s = torch.tensor(self.s, dtype=torch.float32, device=device)
        n = torch.tensor(self.n, dtype=torch.float32, device=device)
        return x, (s + 0.5) / (n + 1.0), n


class Runner:
    """Counts every stochastic run; each run gets a fresh, never-reused seed slot."""

    def __init__(self, simulate, device, seed):
        self.simulate = simulate
        self.device = device
        self.runs = 0
        self.gen = torch.Generator(device=device)
        self.gen.manual_seed(seed)

    def __call__(self, x_np):
        x = torch.tensor(x_np, dtype=torch.float32, device=self.device)
        self.runs += x.shape[0]
        return self.simulate(x, generator=self.gen).cpu().numpy()


def sequential_test(runner, xs, round_seeds=4, max_seeds=24, conf=0.95):
    """Rule 6: add seeds in rounds until P(p>0.5) is confidently high or low."""
    m = len(xs)
    s = np.zeros(m)
    n = np.zeros(m)
    active = np.arange(m)
    while active.size:
        batch = np.repeat(xs[active], round_seeds, axis=0)
        y = runner(batch).reshape(active.size, round_seeds)
        s[active] += y.sum(1)
        n[active] += round_seeds
        post = beta.sf(0.5, 1 + s[active], 1 + n[active] - s[active])
        done = (post > conf) | (post < 1 - conf) | (n[active] >= max_seeds)
        active = active[~done]
    return s, n


class OrderedEngine:
    def __init__(self, simulate, dim, monotone, device, seed=0, use_sprt=True,
                 use_active=True, use_monotone=True, fixed_seeds=8,
                 n_init=256, batch=256, explore_frac=0.2, delta=0.08,
                 spot_check=0.1, pool_size=1 << 18):
        self.dim = dim
        self.mono = np.asarray(monotone)
        self.device = device
        self.run = Runner(simulate, device, seed)
        self.rng = np.random.default_rng(seed)
        self.sobol = qmc.Sobol(d=dim, scramble=True, seed=seed)
        self.archive = Archive()
        self.use_sprt, self.use_active, self.use_monotone = use_sprt, use_active, use_monotone
        self.fixed_seeds = fixed_seeds
        self.n_init, self.batch = n_init, batch
        self.explore_frac, self.delta, self.spot_check = explore_frac, delta, spot_check
        self.pool = self.sobol.random(pool_size).astype(np.float32)
        self.pool_used = np.zeros(pool_size, dtype=bool)
        self.next_seq = 0  # next pool index in plain Sobol order
        self.inferred_count = 0
        self.spot_checks = 0
        self.violations = 0
        self.history = []

    # ---- rule 3 -------------------------------------------------------
    def _infer(self, xs):
        """Return +1/-1/0 per candidate from dominance by confident results."""
        if not self.use_monotone or not self.archive.x:
            return np.zeros(len(xs), dtype=int)
        ax = np.stack(self.archive.x)
        s = np.asarray(self.archive.s, dtype=float)
        n = np.asarray(self.archive.n, dtype=float)
        sim = ~np.asarray(self.archive.inferred)
        if self.use_sprt:
            post = beta.sf(0.5, 1 + s, 1 + n - s)
        else:  # single-seed results: judge confidence from the smoothed map
            tx, tp, tw = self.archive.tensors(self.device)
            post = knn_predict(tx, tp, tw, tx).cpu().numpy()
        decl = self.mono != 0
        und = ~decl
        sign = self.mono[decl]
        out = np.zeros(len(xs), dtype=int)
        for label, mask in ((1, sim & (post > 0.95)), (-1, sim & (post < 0.05))):
            ref = ax[mask]
            if not len(ref):
                continue
            xd = torch.tensor(xs[:, decl] * sign)
            rd = torch.tensor(ref[:, decl] * sign)
            if label == 1:  # candidate at least as favourable as a success
                ok = (xd[:, None, :] >= rd[None, :, :]).all(-1)
            else:           # candidate at most as favourable as a failure
                ok = (xd[:, None, :] <= rd[None, :, :]).all(-1)
            if und.any():
                near = (torch.tensor(xs[:, und])[:, None, :] - torch.tensor(ref[:, und])[None, :, :]).abs() <= self.delta
                ok &= near.all(-1)
            hit = ok.any(1).numpy() & (out == 0)
            out[hit] = label
        return out

    # ---- evaluate a set of candidate conditions ------------------------
    def _test(self, xs):
        xs = np.array([x for x in xs if x not in self.archive], dtype=np.float32)
        if not len(xs):
            return
        inf = self._infer(xs)
        check = (inf != 0) & (self.rng.random(len(xs)) < self.spot_check)
        for x, lab in zip(xs[(inf != 0) & ~check], inf[(inf != 0) & ~check]):
            self.archive.add(x, 4 if lab > 0 else 0, 4, inferred=True)
            self.inferred_count += 1
        sim_mask = (inf == 0) | check
        xs_sim = xs[sim_mask]
        if not len(xs_sim):
            return
        if self.use_sprt:
            s, n = sequential_test(self.run, xs_sim)
        else:
            y = self.run(np.repeat(xs_sim, self.fixed_seeds, axis=0)).reshape(-1, self.fixed_seeds)
            s, n = y.sum(1).astype(float), np.full(len(xs_sim), float(self.fixed_seeds))
        lab_sim = inf[sim_mask]
        for x, si, ni, lab in zip(xs_sim, s, n, lab_sim):
            self.archive.add(x, si, ni)
            if lab != 0:
                self.spot_checks += 1
                post = beta.sf(0.5, 1 + si, 1 + ni - si)
                if not self.use_sprt:
                    post = si / ni
                if (lab > 0 and post < 0.05) or (lab < 0 and post > 0.95):
                    self.violations += 1
        if self.violations >= 3 and self.use_monotone:
            self.use_monotone = False  # assumption refuted: stop relying on it

    def _take_sequential(self, m):
        idx = []
        while len(idx) < m and self.next_seq < len(self.pool):
            if not self.pool_used[self.next_seq]:
                idx.append(self.next_seq)
            self.next_seq += 1
        return np.array(idx, dtype=int)

    def _take_active(self, m):
        free = np.flatnonzero(~self.pool_used)
        tx, tp, tw = self.archive.tensors(self.device)
        q = torch.tensor(self.pool[free], device=self.device)
        p = knn_predict(tx, tp, tw, q).cpu().numpy()
        score = p * (1 - p) + 1e-3 * self.rng.random(len(free))
        return free[np.argsort(-score)[:m]]

    def step(self):
        if not self.archive.x:
            idx = self._take_sequential(self.n_init)
        elif self.use_active:
            # batches grow with the archive so the number of rounds stays ~log(budget)
            batch = max(self.batch, len(self.archive.x) // 4)
            n_exp = int(batch * self.explore_frac)
            idx = self._take_sequential(n_exp)
            self.pool_used[idx] = True
            idx = np.concatenate([idx, self._take_active(batch - n_exp)])
        else:
            idx = self._take_sequential(self.batch)
        self.pool_used[idx] = True
        self._test(self.pool[idx])

    def predict(self, query):
        tx, tp, tw = self.archive.tensors(self.device)
        return knn_predict(tx, tp, tw, query)
