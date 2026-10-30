"""Physics-only coupling check (coordinator 29 Oct 23:12, item a): does a vesicle's internal osmotic load Pi (its
enclosed letters per enclosed volume, relative to the outside: the engine's step-6 rule, membrane uptake share
proportional to a (1 + beta Pi)) predict its membrane growth and its divisions, beyond membrane area and age? CPU only.
Darwin demo (as darwin_long1: retention trait on), fair random rationing, external division, 16 sites x 256 slots,
64 vesicle slots, L_MAX 1024, seed 1. Every 5 steps after year 10 each alive vesicle's Pi, membrane a, age and
cumulative divisions are recorded; over the next interval: divisions (same vesicle line at the slot) and the relative
membrane change of non-dividing vesicles. Partial correlations (controlling for log a and age), real vs control half;
mother-daughter Pi correlation (heritability of the load).
    python tests/osmotic_coupling.py out.json [wall seconds]"""
import json, os, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gokyona_env", "engine"))
import gk_chem  # noqa: E402
from trait import Trait  # noqa: E402

torch.set_num_threads(4)
OUT = sys.argv[1]
WALL = float(sys.argv[2]) if len(sys.argv) > 2 else 150
PRE = json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
N, S = 16, 256
gk_chem.configure(16, 64)
NV = gk_chem.NV
A0 = gk_chem.A0_VES


class Hook:
    def __init__(self):
        self.t0, self.k, self.prev, self.first = time.time(), 0, None, {}
        self.rec = []            # (half, pi, log_a, age, ndiv, dlog_a or nan)
        self.md = []

    def due(self):
        return True

    def __call__(self, loc):
        self.k += 1
        VS = loc["VS"]; L, lab = loc["L"], loc["lab"]
        al = VS.alive.cpu()
        for r_, j_ in al.nonzero().tolist():
            self.first.setdefault((r_, int(VS.vid[r_, j_])), self.k)
        if float(loc["t"][0]) / 8766 < 10 or self.k % 5:
            return "stop" if time.time() - self.t0 > WALL else None
        R = L.shape[0]
        ench = VS._letters_in(L, lab).cpu()
        outl = torch.where((lab == 0) & (L > 0), L, torch.zeros_like(L)).sum(1).float().cpu()
        phi = (VS.phi0[:, None] * (VS.va / A0) ** 1.5).clamp(max=1.0).cpu() * al
        outvol = (1 - phi.sum(1)).clamp(min=0.05)
        Pi = ench / phi.clamp(min=1e-30) / (outl / outvol + 1.0)[:, None]
        va = VS.va.cpu(); vid = VS.vid.cpu(); vpar = VS.vpar.cpu(); vdiv = VS.vdiv.cpu()
        cur = {}
        for r_, j_ in al.nonzero().tolist():
            key = (r_, j_)
            cur[key] = (float(Pi[r_, j_]), float(va[r_, j_]), int(vid[r_, j_]), int(vpar[r_, j_]), float(vdiv[r_, j_]),
                        self.k - self.first.get((r_, int(vid[r_, j_])), self.k))
        if self.prev is not None:
            for key, (pi0, a0, id0, par0, d0, age0) in self.prev.items():
                c = cur.get(key)
                if c is None:
                    continue
                pi1, a1, id1, par1, d1, _ = c
                if not (id1 == id0 or par1 == id0):
                    continue                                # slot taken by another vesicle line
                nd = d1 - d0
                dla = float(np.log(max(a1, 1e-12)) - np.log(max(a0, 1e-12))) if nd == 0 else float("nan")
                half = "control" if (key[0] % N) >= N // 2 else "real"
                if a0 > 0 and pi0 > 0:
                    self.rec.append((half, pi0, float(np.log(a0)), age0, nd, dla))
            # mother-daughter Pi: daughters born in the interval (new id, parent id alive at the last record)
            pid = {(k_[0], v[2]): v[0] for k_, v in self.prev.items()}
            for key, (pi1, a1, id1, par1, d1, age1) in cur.items():
                if (key[0], id1) not in pid and par1 >= 0 and (key[0], par1) in pid and pi1 > 0:
                    self.md.append((pid[(key[0], par1)], pi1))
        self.prev = cur
        return "stop" if time.time() - self.t0 > WALL else None


def partial(y, x, Z):
    Zc = np.c_[np.ones(len(y)), Z]
    ry = y - Zc @ np.linalg.lstsq(Zc, y, rcond=None)[0]
    rx = x - Zc @ np.linalg.lstsq(Zc, x, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1]) if rx.std() > 0 and ry.std() > 0 else None


u = torch.tensor([PRE["u"]], dtype=torch.float32)
ov = {k: float(v) for k, v in PRE["env"].items() if k != "avail_time"}
ov.update({"mineral_type": float(round(ov.get("mineral_type", 0))), "endless": True, "ration": "random", "vesicle_stats": True})
tr = Trait(dict(PRE["trait"], enabled=True), N, torch.device("cpu"))
hk = Hook()
g = torch.Generator(); g.manual_seed(1)
gk_chem.simulate(u, g, (N, S, 300), override=ov, trait=tr, hook=hk)
p_ = gk_chem.decode(u)
res = {"beta": float(p_["beta"][0]), "phi0": float(p_["phi0"][0]), "records": len(hk.rec), "halves": {}}
for half in ("real", "control"):
    R_ = [r for r in hk.rec if r[0] == half]
    if len(R_) < 50:
        continue
    pi = np.log(np.array([r[1] for r in R_])); la = np.array([r[2] for r in R_]); age = np.array([r[3] for r in R_], float)
    nd = np.array([r[4] for r in R_], float); dla = np.array([r[5] for r in R_])
    Z = np.c_[la, age]
    m = ~np.isnan(dla)
    res["halves"][half] = {
        "n": len(R_), "log_pi_sd": float(pi.std()), "division_share": float((nd > 0).mean()),
        "corr_logpi_divisions": float(np.corrcoef(pi, nd)[0, 1]) if nd.std() > 0 else None,
        "partial_logpi_divisions_given_loga_age": partial(nd, pi, Z) if nd.std() > 0 else None,
        "corr_loga_divisions": float(np.corrcoef(la, nd)[0, 1]) if nd.std() > 0 else None,
        "partial_logpi_membrane_growth_given_loga_age": partial(dla[m], pi[m], Z[m]) if m.sum() > 50 else None,
        "n_nondividing": int(m.sum())}
md = np.array(hk.md) if hk.md else np.zeros((0, 2))
res["mother_daughter_logpi_corr"] = float(np.corrcoef(np.log(md[:, 0]), np.log(md[:, 1]))[0, 1]) if len(md) > 10 else None
res["mother_daughter_n"] = int(len(md))
json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps(res, indent=1))
