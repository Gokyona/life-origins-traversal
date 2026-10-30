"""Sources of fresh (unselected) domains per step in the Darwin demo (coordinator request, 29 Oct 2026 12:26). CPU only.
Short run (Darwin demo, 16 sites x 512 slots, 64 vesicle slots, L_MAX 1024, seed 1). At every step the hook gives each
chain a tag = its slot + 1 (the tag is write-only for the engine; fragments and copies inherit their parent's tag) and
stores each scored chain's z (real landscape, pair score). At the next step every scored chain (chemistry 0, >= 12
letters) of the real half is classified:
  same slot, same tag : old chain (a left piece shortened below 24 letters is counted as "left piece truncated")
  other slot, paired with the parent's slot : copy (complement strand; partial if shorter than 24 letters)
  other slot, otherwise                    : right fragment of a bond break (its domain = letters that were inside the
                                             parent), or a copy whose duplex dissolved in the same step
  tag 0                                    : grown de novo
For each class: chains per step, their mean z, the parent's z, and the change of the real-half chain mean z per step
they cause by entering (n x (z_class - z_mean) / N); compared with the selection differential from washout.
    python tests/fresh_domains.py out.json [wall seconds] [skip years]"""
import json, os, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gokyona_env", "engine"))
import gk_chem  # noqa: E402
from trait import Trait  # noqa: E402

torch.set_num_threads(4)
OUT = sys.argv[1]
WALL = float(sys.argv[2]) if len(sys.argv) > 2 else 200
SKIP_Y = float(sys.argv[3]) if len(sys.argv) > 3 else 8
PRE = json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
YEAR = 8766.0
NC, SC = 16, 512
gk_chem.configure(16, 64)
NV = gk_chem.NV
u = torch.tensor([PRE["u"]], dtype=torch.float32)
ov = {k: float(v) for k, v in PRE["env"].items() if k != "avail_time"}
ov["mineral_type"] = float(round(ov.get("mineral_type", 0)))
ov["endless"] = True
trC = Trait(dict(PRE["trait"], enabled=True), NC, torch.device("cpu"))
trS = Trait(dict(PRE["trait"], enabled=True), NC, torch.device("cpu"))
CLS = ("old", "left_truncated", "copy_full", "copy_partial", "fragment", "de_novo")


class Hook:
    def __init__(self):
        self.k, self.t0 = 0, time.time()
        self.prev_z = None
        self.acc = {c: {"n": 0, "z": 0.0, "zpar": 0.0, "npar": 0, "contrib": 0.0, "z_lt0": 0} for c in CLS}
        self.steps_used, self.N_sum, self.zmean_sum, self.sel = 0, 0, 0.0, []
        self.copy_err = {"n": 0, "dz_ne0": 0}

    def due(self):
        return True

    def __call__(self, loc):
        self.k += 1
        L, seq, chem, lab, pair, ft, VS = loc["L"], loc["seq"], loc["chem"], loc["lab"], loc["pair"], loc["ftag"], loc["VS"]
        R, S = L.shape
        t = float(loc["t"][0])
        real = ((torch.arange(R) % NC) < NC // 2)[:, None].expand(R, S)
        sc = (chem == 0) & (L >= 12) & real
        rr, cc = sc.nonzero(as_tuple=True)
        z = trS.zscore(L[rr, cc], seq[rr, cc]) if rr.numel() else torch.zeros(0)
        zfull = torch.full((R * S,), float("nan"))
        zfull[rr * S + cc] = z
        if self.prev_z is not None and t / YEAR >= SKIP_Y and rr.numel():
            g = rr * S + cc
            tag = ft[rr, cc]
            orig = tag - 1
            zm = float(z.mean()); N = int(z.numel())
            same = orig == g
            newc = (tag > 0) & ~same
            pr_col = pair[rr, cc]
            is_copy = newc & ((orig // S) == rr) & (pr_col == (orig % S))
            frag = newc & ~is_copy
            Lnow = L[rr, cc]
            trunc = same & (self.prev_L[g] >= 24) & (Lnow < 24)
            cls = {"old": same & ~trunc, "left_truncated": trunc, "copy_full": is_copy & (Lnow >= 24),
                   "copy_partial": is_copy & (Lnow < 24), "fragment": frag, "de_novo": tag == 0}
            zpar_all = self.prev_z[orig.clamp(min=0)]
            for c, m in cls.items():
                if c == "old" or not bool(m.any()):
                    continue
                a = self.acc[c]
                n = int(m.sum())
                a["n"] += n; a["z"] += float(z[m].sum()); a["z_lt0"] += int((z[m] < 0).sum())
                zp = zpar_all[m]; ok = ~torch.isnan(zp)
                a["npar"] += int(ok.sum()); a["zpar"] += float(zp[ok].sum())
                a["contrib"] += n * (float(z[m].mean()) - zm) / N
            cf = cls["copy_full"]
            if bool(cf.any()):
                zp = zpar_all[cf]; ok = ~torch.isnan(zp)
                self.copy_err["n"] += int(ok.sum()); self.copy_err["dz_ne0"] += int(((z[cf][ok] - zp[ok]).abs() > 1e-5).sum())
            self.steps_used += 1; self.N_sum += N; self.zmean_sum += zm
            # selection differential from washout this step (chain-weighted and vesicle-weighted), real half
            if VS.wash_mult is not None:
                enc = lab[rr, cc] > 0
                comp = rr[enc] * NV + lab[rr, cc][enc].long() - 1
                n_v = torch.zeros(R * NV).index_add_(0, comp, torch.ones(int(enc.sum())))
                s_v = torch.zeros(R * NV).index_add_(0, comp, z[enc])
                has = n_v > 0
                zb = s_v[has] / n_v[has]
                dt = float(loc["dt"][0])
                pw = 1 - torch.exp(-(VS.w[:, None] * dt).expand(R, NV).reshape(-1)[has] * VS.wash_mult.reshape(-1)[has])
                sv = 1 - pw
                self.sel.append(((float((sv * zb).sum() / sv.sum() - zb.mean())),
                                 float((sv * n_v[has] * zb).sum() / (sv * n_v[has]).sum() - (n_v[has] * zb).sum() / n_v[has].sum()),
                                 float(pw.mean())))
        # re-tag every chain with its own slot; remember z and length
        ft.copy_(torch.where(L > 0, torch.arange(1, R * S + 1).view(R, S), torch.zeros_like(ft)))
        self.prev_z, self.prev_L = zfull, L.reshape(-1).clone()
        if self.k % 100 == 0:
            print("step", self.k, "years", round(t / YEAR, 1), round(time.time() - self.t0), "s", flush=True)
        if time.time() - self.t0 > WALL:
            return "stop"


hk = Hook()
g = torch.Generator(); g.manual_seed(1)
gk_chem.simulate(u, g, (NC, SC, 300), override=ov, trait=trC, hook=hk)
st = max(1, hk.steps_used)
res = {"size": f"{NC} sites x {SC} slots, 64 vesicle slots, L_MAX 1024, Darwin demo, seed 1", "steps_used": hk.steps_used,
       "skip_years": SKIP_Y, "scored_chains_real_mean": hk.N_sum / st, "z_mean_real": hk.zmean_sum / st, "classes": {}}
for c in CLS:
    a = hk.acc[c]
    if c == "old":
        continue
    res["classes"][c] = {"per_step": a["n"] / st, "per_step_share_of_scored": a["n"] / max(1, hk.N_sum),
                         "z_mean": a["z"] / a["n"] if a["n"] else None, "share_z_below_0": a["z_lt0"] / a["n"] if a["n"] else None,
                         "parent_z_mean": a["zpar"] / a["npar"] if a["npar"] else None,
                         "dmean_z_per_step": a["contrib"] / st}
sel = np.array(hk.sel) if hk.sel else np.zeros((1, 3))
res["selection_diff_per_step"] = {"vesicle_weighted": float(sel[:, 0].mean()), "chain_weighted": float(sel[:, 1].mean()), "washout_p": float(sel[:, 2].mean())}
res["full_copies_with_dz"] = {"n": hk.copy_err["n"], "share_z_differs_from_template": hk.copy_err["dz_ne0"] / max(1, hk.copy_err["n"])}
res["note"] = ("dmean_z_per_step = n x (z_class - z_mean) / N_scored per step, i.e. the shift of the real-half chain mean z caused by "
               "these chains entering; compare with the chain-weighted selection differential (both chain-level)")
json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps(res, indent=1))
