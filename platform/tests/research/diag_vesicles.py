"""Why membrane mode does not change reproduction (coordinator, 29 Oct 08:33): per vesicle (id) log membrane uptake,
divisions (as parent), z_mean of enclosed scored chains, enclosed chain count; per step the division probability of
eligible vesicles (1 - exp(-k_div dt)) and the share of vesicles that are eligible (a >= 2 a0).
    python tests/diag_vesicles.py [steps] [device]"""
import json, os, sys
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gokyona_env", "engine"))
import gk_chem, devices
from trait import Trait

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 150
dev = devices.pick(sys.argv[2] if len(sys.argv) > 2 else "cpu")
gk_chem.configure(4)
pre = json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
N, S = 16, 256
tr = Trait({"enabled": True, "mode": "membrane", "s": 4, "landscape": "domain", "pair_score": True}, N, dev)
log = {"uptake": {}, "div": {}, "z": {}, "nch": {}, "steps": {}}
per_step = {"p_div": [], "elig_share": [], "alive": [], "trig_size_limited": 0, "trig_prob_limited": 0}
V = gk_chem.Vesicles
orig_step, orig_div = V.step, V._divide


def step(self, t, dt, ph, L, seq, chem, pair, dlen, lab, unpair, Ms):
    va0, vid0, al0 = self.va.clone(), self.vid.clone(), self.alive.clone()
    self._divd = torch.zeros_like(self.alive)
    dtr = dt.to(self.dev, torch.float32)[self.rr]
    pdiv = (1 - torch.exp(-self.kdiv * dtr))
    elig = al0 & (va0 >= 2 * gk_chem.A0_VES)
    per_step["p_div"].append(float(pdiv.mean())); per_step["elig_share"].append(float(elig.float().sum() / al0.float().sum().clamp(min=1)))
    per_step["alive"].append(float(al0.sum()) / self.R)
    orig_step(self, t, dt, ph, L, seq, chem, pair, dlen, lab, unpair, Ms)
    zm, has = tr.vesicle_z(L, chem, lab, self.R) if tr._all is not None else (torch.zeros_like(va0), torch.zeros_like(al0))
    nch = torch.zeros(self.R * 8, device=self.dev).index_add_(0, ((torch.arange(self.R, device=self.dev)[:, None] * 8 + lab.long() - 1)[lab > 0]),
                                                             torch.ones(int((lab > 0).sum()), device=self.dev)).view(self.R, 8)
    gain = torch.where(self._divd, self.va * 2 - va0, self.va - va0)        # a division halves the parent
    rr, cc = al0.nonzero(as_tuple=True)
    for r, c in zip(rr.tolist(), cc.tolist()):
        k = (r, int(vid0[r, c]))
        log["uptake"][k] = log["uptake"].get(k, 0.0) + float(gain[r, c])
        log["div"][k] = log["div"].get(k, 0) + int(self._divd[r, c])
        log["z"].setdefault(k, []).append(float(zm[r, c]) if bool(has[r, c]) else np.nan)
        log["nch"].setdefault(k, []).append(float(nch[r, c]))
        log["steps"][k] = log["steps"].get(k, 0) + 1


def _divide(self, t, div, L, lab, pair):
    self._divd = self._divd | div
    orig_div(self, t, div, L, lab, pair)


V.step, V._divide = step, _divide


class Stop:
    k = 0

    def due(self):
        return True

    def __call__(self, loc):
        Stop.k += 1
        return "stop" if Stop.k >= steps else None


u = torch.tensor([pre["u"]], dtype=torch.float32, device=dev)
ov = {k: float(v) for k, v in pre["env"].items() if k != "avail_time"}
ov["endless"] = True
g = torch.Generator(device=dev); g.manual_seed(1)
gk_chem.simulate(u, g, (N, S, 300), override=ov, trait=tr, hook=Stop())
ks = [k for k in log["steps"] if log["steps"][k] >= 3]
up = np.array([log["uptake"][k] / log["steps"][k] for k in ks]); dv = np.array([log["div"][k] / log["steps"][k] for k in ks])
z = np.array([np.nanmean(log["z"][k]) if np.isfinite(log["z"][k]).any() else np.nan for k in ks]); nc = np.array([np.mean(log["nch"][k]) for k in ks])
m = np.isfinite(z)
sp = lambda a, b: float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1]) if len(a) > 3 and np.std(a) > 0 and np.std(b) > 0 else None
print(json.dumps({"vesicles_tracked": len(ks), "mean_steps_per_vesicle": float(np.mean([log["steps"][k] for k in ks])),
                  "division_prob_per_step_if_eligible_mean": float(np.mean(per_step["p_div"])),
                  "eligible_share_of_alive_mean": float(np.mean(per_step["elig_share"])),
                  "alive_vesicles_per_site_mean": float(np.mean(per_step["alive"])),
                  "divisions_per_step_mean": float(dv.mean()), "uptake_per_step_mean": float(up.mean()),
                  "spearman_div_vs_uptake": sp(dv, up), "spearman_div_vs_z": sp(dv[m], z[m]), "spearman_uptake_vs_z": sp(up[m], z[m]),
                  "enclosed_chains_per_vesicle_mean": float(nc.mean()), "share_vesicles_with_scored_chains": float(m.mean())}, indent=1))
