"""Proxy for clades from a COPY of a checkpoint (CPU only): the 24-letter functional domains of nucleotide-like chains
(>= 24 letters), strand-canonical (min of domain and its complement, since copies are complements), per half; for the
200 most common domains, the share of the half's chains within Hamming distance <= 1, 2, 3 (a "ball"); reports the
largest ball per radius and the share of the single most common domain.
    python tests/domain_clusters.py <checkpoint copy> <sites>"""
import sys
import numpy as np
import torch
ck = torch.load(sys.argv[1], map_location="cpu", weights_only=False)
N = int(sys.argv[2])
L, seq, chem = ck["L"], ck["seq"], ck["chem"]
R, S = L.shape
m = (L >= 24) & (chem == 0)
rr, cc = m.nonzero(as_tuple=True)
dom = (seq[rr, cc, 0] & ((1 << 24) - 1)).numpy().astype(np.int64)
comp = (~dom) & ((1 << 24) - 1)
canon = np.minimum(dom, comp)
half = ((rr % N) >= N // 2).numpy()
POP = np.array([bin(i).count("1") for i in range(1 << 12)], np.int64)
pc = lambda x: POP[x & 0xFFF] + POP[(x >> 12) & 0xFFF]
out = {"sim_years": float(ck["t"][0]) / 8766}
for name, hm in (("real", ~half), ("control", half)):
    d = canon[hm]
    n = len(d)
    u, cnt = np.unique(d, return_counts=True)
    top = u[np.argsort(-cnt)[:200]]
    best = {}
    for r in (0, 1, 2, 3):
        shares = [(pc(d ^ t) <= r).mean() for t in top]
        best[r] = float(max(shares))
    out[name] = {"chains_ge24": int(n), "distinct_domains": int(len(u)), "top1_domain_share": float(cnt.max() / n),
                 "max_ball_share_r1": best[1], "max_ball_share_r2": best[2], "max_ball_share_r3": best[3]}
print(out)
