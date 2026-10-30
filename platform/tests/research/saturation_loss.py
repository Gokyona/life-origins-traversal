"""Paper-model diagnostics at the 5a points (coordinator request, 29 Oct 08:03): gk_chem with W = 1 (bit-identical to
chem4c4), paper settings (screen fidelity 64 x 512 x 300, the paper's dt, no trait, single-start copying), one run per
point, seed 1. Reports per point: the fraction of template-steps with lambda dt > 3 and > 1, the median lambda dt
(from a 0.1-decade histogram), and the loss budget of lineage chains (chemistry 0, L >= L_H): bond breaks, washout of
free chains, loss inside washed-out vesicles, and fragments lost for want of a slot (all chains).
    python tests/saturation_loss.py out.json [device]"""
import json, os, sys, time
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gokyona_env", "engine"))
import gk_chem, devices  # noqa: E402

out_path = sys.argv[1]
dev = devices.pick(sys.argv[2] if len(sys.argv) > 2 else "auto")
gk_chem.configure(1)
pre = json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))
LOT = os.environ.get("GK_CHEM4C4_DIR")          # paper repository, code/lot-s4/ (confirm4c4v2.json, scan4c4mem_sets.py)
if not LOT:
    sys.exit("set GK_CHEM4C4_DIR to the paper repository folder code/lot-s4/")
c2 = json.load(open(os.path.join(LOT, "confirm4c4v2.json")))["points"]
sys.path.insert(0, LOT)
import scan4c4mem_sets as S9  # noqa: E402
P = {"engine12538": pre["engine12538"]["u"]}
for k, v in S9.first_sets().items():
    P[k + "_set0"] = v
for k, v in c2.items():
    P["9c_" + k] = v["u_model"]
names = list(P)
u = torch.tensor([P[n] for n in names], dtype=torch.float32, device=dev)


class Keep:
    loc = None

    def due(self):
        return True

    def __call__(self, loc):
        Keep.loc = loc


g = torch.Generator(device=dev); g.manual_seed(1)
t0 = time.time()
o = gk_chem.simulate(u, g, "screen", hook=Keep())
loc = Keep.loc
res = {"device": devices.describe(dev), "seconds": round(time.time() - t0, 1), "points": {}}
edges = [10 ** (-6 + 0.1 * (i + 0.5)) for i in range(121)]
for b, n in enumerate(names):
    h = loc["sat_hist"][b]
    tot = float(h.sum())
    cum = torch.cumsum(h, 0)
    med = edges[int(torch.searchsorted(cum, torch.tensor(tot / 2, dtype=torch.float64)))] if tot else None
    lv = float(loc["VS"].loss_ves.view(len(names), -1)[b].sum()) if hasattr(loc["VS"], "loss_ves") else 0.0
    lh, lw, ov = float(loc["loss_hyd"][b]), float(loc["loss_wash"][b]), float(loc["overflow"][b])
    tl = lh + lw + lv
    res["points"][n] = {"template_steps": float(loc["sat_tot"][b]),
                        "frac_lambda_dt_gt3": float(loc["sat3"][b]) / max(float(loc["sat_tot"][b]), 1),
                        "frac_lambda_dt_gt1": float(loc["cp_sat"][b]) / max(float(loc["sat_tot"][b]), 1),
                        "median_lambda_dt": med, "dt_target_h": float(loc["dt_target"][b]), "t_end_h": float(loc["t_end"][b]),
                        "loss_lineage_chains": {"bond_breaks": lh, "washout_free": lw, "vesicle_washout": lv,
                                                "share_breaks": lh / tl if tl else None, "share_washout": lw / tl if tl else None,
                                                "share_vesicle": lv / tl if tl else None},
                        "fragments_lost_no_slot": ov, "copies": float(o["copies"][b])}
    print(n, json.dumps(res["points"][n]), flush=True)
json.dump(res, open(out_path, "w"), indent=1)
