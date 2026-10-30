"""Why are low-z templates copied more? (coordinator request, 29 Oct 2026 12:31). CPU only. Short run (Darwin demo,
16 sites x 512 slots, 64 vesicle slots, L_MAX 1024, seed 1). At every step the hook tags every chain with its slot
(write-only for the engine; copies inherit their template's tag) and records, for a sample of scored real-half
templates (chemistry 0, >= 12 letters, single-stranded, enclosed): z, length, slot position in the site, vesicle mean z,
vesicle age (steps since its id appeared), chains in the vesicle, activated monomers of the site, chains in the site.
Outcome: copies of that template made in the next step. Engine fact used: monomers are a per-site pool and started
copies are accepted in slot order while the cumulative need fits (gk_chem: within_c <= Ms).
    python tests/copy_bias.py out.json [wall seconds] [skip years] [engine options as JSON]"""
import json, os, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gokyona_env", "engine"))
import gk_chem  # noqa: E402
from trait import Trait  # noqa: E402

torch.set_num_threads(4)
OUT = sys.argv[1]
WALL = float(sys.argv[2]) if len(sys.argv) > 2 else 180
SKIP_Y = float(sys.argv[3]) if len(sys.argv) > 3 else 8
EXTRA = json.loads(sys.argv[4]) if len(sys.argv) > 4 else {}      # engine options, e.g. {"ration": "random"}
PRE = json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
YEAR = 8766.0
NC, SC = 16, 512
gk_chem.configure(16, 64)
NV = gk_chem.NV
u = torch.tensor([PRE["u"]], dtype=torch.float32)
ov = {k: float(v) for k, v in PRE["env"].items() if k != "avail_time"}
ov["mineral_type"] = float(round(ov.get("mineral_type", 0)))
ov["endless"] = True
ov.update(EXTRA)
trC = Trait(dict(PRE["trait"], enabled=True), NC, torch.device("cpu"))
trS = Trait(dict(PRE["trait"], enabled=True), NC, torch.device("cpu"))
FEAT = ["z", "length", "slot_pos", "ves_zbar", "ves_age_steps", "ves_chains", "site_monomers", "site_chains"]
rng = np.random.default_rng(3)


class Hook:
    def __init__(self):
        self.k, self.t0 = 0, time.time()
        self.prev = None
        self.first_seen = {}
        self.X, self.Y = [], []

    def due(self):
        return True

    def __call__(self, loc):
        self.k += 1
        L, seq, chem, lab, pair, ft, VS, Ms = loc["L"], loc["seq"], loc["chem"], loc["lab"], loc["pair"], loc["ftag"], loc["VS"], loc["Ms"]
        R, S = L.shape
        t = float(loc["t"][0])
        # outcome for the previous step's sample: copies = chains in another slot carrying the template's tag, paired with it
        if self.prev is not None:
            gidx, feats = self.prev
            tag = ft.reshape(-1)
            occ = (L.reshape(-1) > 0)
            orig = tag - 1
            pos = torch.arange(R * S)
            prc = pair.reshape(-1)
            isc = occ & (tag > 0) & (orig != pos) & ((orig // S) == (pos // S)) & (prc == (orig % S))
            cnt = torch.zeros(R * S).index_add_(0, orig[isc], torch.ones(int(isc.sum())))
            self.X.append(feats); self.Y.append(cnt[gidx].numpy())
            self.prev = None
        # vesicle ages
        al = VS.alive.nonzero(as_tuple=True)
        for r_, v_ in zip(al[0].tolist(), VS.vid[al].tolist()):
            self.first_seen.setdefault((r_, v_), self.k)
        if t / YEAR >= SKIP_Y:
            real = ((torch.arange(R) % NC) < NC // 2)[:, None].expand(R, S)
            sc = (chem == 0) & (L >= 12) & real & (lab > 0) & (pair < 0)
            rr, cc = sc.nonzero(as_tuple=True)
            if rr.numel():
                z = trS.zscore(L[rr, cc], seq[rr, cc])
                # vesicle stats over all scored chains (paired or not)
                sc2 = (chem == 0) & (L >= 12) & (lab > 0)
                r2, c2 = sc2.nonzero(as_tuple=True)
                z2 = trS.zscore(L[r2, c2], seq[r2, c2])
                comp2 = r2 * NV + lab[r2, c2].long() - 1
                n_s = torch.zeros(R * NV).index_add_(0, comp2, torch.ones_like(z2))
                s_s = torch.zeros(R * NV).index_add_(0, comp2, z2)
                occ = (L > 0) & (lab > 0)
                r3, c3 = occ.nonzero(as_tuple=True)
                n_all = torch.zeros(R * NV).index_add_(0, r3 * NV + lab[r3, c3].long() - 1, torch.ones(r3.numel()))
                comp = rr * NV + lab[rr, cc].long() - 1
                vid = VS.vid[rr, lab[rr, cc].long() - 1]
                age = torch.tensor([self.k - self.first_seen.get((a, b), self.k) for a, b in zip(rr.tolist(), vid.tolist())], dtype=torch.float32)
                site_ch = (L > 0).sum(1).float()
                f = torch.stack([z, L[rr, cc].float(), cc.float() / S, s_s[comp] / n_s[comp].clamp(min=1), age, n_all[comp],
                                 Ms[rr, 0].float(), site_ch[rr]], 1).numpy()
                pick = rng.choice(len(f), size=min(len(f), 800), replace=False)
                self.prev = ((rr * S + cc)[pick], f[pick])
        ft.copy_(torch.where(L > 0, torch.arange(1, R * S + 1).view(R, S), torch.zeros_like(ft)))
        if self.k % 100 == 0:
            print("step", self.k, "years", round(t / YEAR, 1), round(time.time() - self.t0), "s", flush=True)
        if time.time() - self.t0 > WALL:
            return "stop"


hk = Hook()
g = torch.Generator(); g.manual_seed(1)
gk_chem.simulate(u, g, (NC, SC, 300), override=ov, trait=trC, hook=hk)
X = np.concatenate(hk.X); Y = np.concatenate(hk.Y)
yb = (Y > 0).astype(float)
res = {"size": f"{NC} sites x {SC} slots, Darwin demo, seed 1", "options": EXTRA, "templates": int(len(Y)), "copy_prob_per_step": float(yb.mean()),
       "copies_per_template_per_step": float(Y.mean()), "features": FEAT}
Xs = (X - X.mean(0)) / X.std(0).clip(1e-9)
res["univariate_corr_with_copied"] = {f: float(np.corrcoef(X[:, i], yb)[0, 1]) for i, f in enumerate(FEAT)}
A = np.c_[np.ones(len(Xs)), Xs]
beta = np.linalg.lstsq(A, yb, rcond=None)[0]
res["ols_standardised_beta_copied"] = {f: float(b) for f, b in zip(FEAT, beta[1:])}
pred = A @ beta
res["ols_r2"] = float(1 - ((yb - pred) ** 2).sum() / ((yb - yb.mean()) ** 2).sum())


def partial(i):
    others = [j for j in range(len(FEAT)) if j != i]
    B = np.c_[np.ones(len(Xs)), Xs[:, others]]
    rx = Xs[:, i] - B @ np.linalg.lstsq(B, Xs[:, i], rcond=None)[0]
    ry = yb - B @ np.linalg.lstsq(B, yb, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


res["partial_corr_with_copied"] = {f: partial(i) for i, f in enumerate(FEAT)}
res["feature_corr_with_z"] = {f: float(np.corrcoef(X[:, 0], X[:, i])[0, 1]) for i, f in enumerate(FEAT)}
# copy probability by template-length and slot-position quintiles, and mean z of copied vs not copied
q = lambda a: np.quantile(a, [0.2, 0.4, 0.6, 0.8])
for name, i in (("length", 1), ("slot_pos", 2), ("ves_zbar", 3), ("z", 0)):
    b = np.digitize(X[:, i], q(X[:, i]))
    res[f"copy_prob_by_{name}_quintile"] = [{"mean_" + name: float(X[b == k, i].mean()), "copy_prob": float(yb[b == k].mean()),
                                             "z_mean": float(X[b == k, 0].mean())} for k in range(5)]
res["z_mean_copied"] = float(X[yb > 0, 0].mean()); res["z_mean_not_copied"] = float(X[yb == 0, 0].mean())
json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps(res, indent=1))
