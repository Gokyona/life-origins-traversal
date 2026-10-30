"""Population statistics of a Gokyona run from a COPY of its checkpoint (CPU only; never touches the running job).
    python tests/population_snapshot.py <copied checkpoint.pt> <run config.json> <history.jsonl> <out.json>"""
import json, math, os, pickle, sys
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gokyona_env", "engine"))
import gk_chem, words as wd, gk_cat3 as cat  # noqa: E402
from trait import Trait  # noqa: E402

ck = torch.load(sys.argv[1], map_location="cpu", weights_only=False)
cfg = json.load(open(sys.argv[2]))
hist = [json.loads(l) for l in open(sys.argv[3])]
out_path = sys.argv[4]
gk_chem.configure(int(ck["W"]), int(cfg.get("vesicle_slots", 8)))
L, seq, chem, lab = ck["L"], ck["seq"], ck["chem"], ck["lab"]
R, S = L.shape
N = int(cfg["sites"])
tr = Trait(dict(cfg["trait"], enabled=True), N, torch.device("cpu"))
occ = L > 0
rows = torch.arange(R)[:, None].expand(R, S)
half_ctrl = (rows % N) >= N // 2
q = lambda a: {f"p{p}": float(np.percentile(a, p)) for p in (10, 25, 50, 75, 90)} if len(a) else None
res = {"sim_h": float(ck["t"][0]), "sim_years": float(ck["t"][0]) / 8766.0, "wall_s_at_checkpoint": float(ck.get("wall_s", 0)),
       "steps": int(ck["steps"]), "W": int(ck["W"]), "L_MAX": gk_chem.L_MAX, "sites": N, "slots": S}
Ls = L[occ].numpy()
res["chains"] = {"total": int(occ.sum()), "real_half": int((occ & ~half_ctrl).sum()), "control_half": int((occ & half_ctrl).sum()),
                 "length": {"mean": float(Ls.mean()), "median": float(np.median(Ls)), "p90": float(np.percentile(Ls, 90)),
                            "max": int(Ls.max()), "share_at_LMAX": float((Ls == gk_chem.L_MAX).mean())},
                 "nucleotide_share_letters": float(L[occ & (chem == 0)].sum()) / float(L[occ].sum()),
                 "enclosed_share_chains": float((occ & (lab > 0)).sum()) / float(occ.sum())}
VSd = ck["VS"]
alive = VSd["alive"]
res["vesicles"] = {"per_site": float(alive.sum()) / R, "divisions_so_far": float(VSd["n_div"].sum()) if "n_div" in VSd else None,
                   "generation_max": int(VSd["vgen"][alive].max()), "generation_mean": float(VSd["vgen"][alive].float().mean())}
# lineages (chemistry 0, L >= 6; identity = first 24 letters), per half
long_ = occ & (chem == 0) & (L >= 6)
rr, cc = long_.nonzero(as_tuple=True)
Lk = L[rr, cc].clamp(max=24)
keys = cat.species_key(chem[rr, cc], Lk, seq[rr, cc] & wd.mask(Lk)).numpy()
z = tr.zscore(L[rr, cc], seq[rr, cc]).numpy()
lens = L[rr, cc].numpy()
hc = half_ctrl[rr, cc].numpy()
site = (rr % N).numpy()
lin = {}
for hname, m in (("real", ~hc), ("control", hc)):
    k, zz, ll = keys[m], z[m], lens[m]
    uk, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
    p = cnt / cnt.sum()
    order = np.argsort(-cnt)
    top = []
    for i in order[:10]:
        mm = inv == i
        top.append({"key": "%06x" % (int(uk[i]) & 0xFFFFFF), "chains": int(cnt[i]), "share": float(p[i]),
                    "z_mean": float(zz[mm].mean()), "length_mean": float(ll[mm].mean()), "sites": int(len(set(site[m][mm])))})
    H = float(-(p * np.log(p)).sum())
    lin[hname] = {"chains_ge_6": int(m.sum()), "distinct_lineages": int(len(uk)), "top1_share": float(p[order[0]]),
                  "top5_share": float(p[order[:5]].sum()), "shannon_H": H, "effective_number": float(math.exp(H)),
                  "top10": top}
res["lineages"] = lin
res["trait_z"] = {h: {"mean": float(z[m].mean()), "sd": float(z[m].std()), **q(z[m])} for h, m in (("real", ~hc), ("control", hc))}
# per-site mean z (real sites)
per_site = {}
for s_ in range(N):
    m = site == s_
    if m.any():
        per_site[s_] = float(z[m].mean())
real_sites = [per_site[s_] for s_ in range(N // 2) if s_ in per_site]
ctrl_sites = [per_site[s_] for s_ in range(N // 2, N) if s_ in per_site]
res["per_site_mean_z"] = {"real": {"n": len(real_sites), "mean": float(np.mean(real_sites)), "sd": float(np.std(real_sites)),
                                   "min": float(np.min(real_sites)), "max": float(np.max(real_sites))},
                          "control": {"n": len(ctrl_sites), "mean": float(np.mean(ctrl_sites)), "sd": float(np.std(ctrl_sites)),
                                      "min": float(np.min(ctrl_sites)), "max": float(np.max(ctrl_sites))}}
cpn = float(ck["cp_n"][0])
p = gk_chem.decode(torch.tensor([json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))[cfg["preset"]]["u"]]))
eps = float(p["eps"][0])
res["copying"] = {"copies_total": cpn, "copies_per_chain_present": cpn / max(1, int(occ.sum())),
                  "eps_per_letter": eps, "expected_errors_per_copy_at_mean_length": eps * float(Ls.mean()),
                  "expected_mutations_per_24_letter_domain_per_copy": eps * 24,
                  "note": "errors per copy are not counted in the run; values are expectations from eps"}
# trait over time (history), every ~250 simulated years
ser = []
nxt = 0.0
for h in hist:
    yr = h["sim_h"] / 8766.0
    if yr >= nxt and h.get("score_mean_real") is not None:
        ser.append({"sim_years": round(yr, 1), "z_real": h["score_mean_real"], "z_control": h["score_mean_ctl"],
                    "distinct": h["distinct"], "mean_len": round(h["mean_len"], 1), "vesicles_per_site": round(h["vesicles"], 2)})
        nxt = (math.floor(yr / 250) + 1) * 250
res["trait_over_time_every_250y"] = ser
trk = pickle.loads(ck["tracker"])
res["turnover"] = {"sweeps_recorded": len(trk.sweeps), "dominant_lineage_changes": None,
                   "oldest_lineage_age": None,
                   "note": ("the run records sweeps (new dominant lineage with >= 30% share) only; changes of a dominant "
                            "lineage below that share and the first-seen time of every lineage are not recorded, so "
                            "dominant turnover and lineage ages cannot be reconstructed from the run's files. The top lineage "
                            "holds ~1-2% of chains in each half (see top1_share), i.e. no lineage dominates."),
                   "tracker_nodes": list(trk.nodes.values())}
json.dump(res, open(out_path, "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("lineages", "trait_over_time_every_250y", "turnover")}, indent=1))
for h in ("real", "control"):
    x = lin[h]
    print(h, {k: v for k, v in x.items() if k != "top10"})
    for t in x["top10"]:
        print("   ", t)
print("over time", ser)
