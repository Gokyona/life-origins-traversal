"""First-shift analysis for the Darwin demo (coordinator request, 29 Oct 2026 11:58). CPU only; reads a COPY of the
darwin_long1 checkpoint and the run histories, and makes one short small CPU run with per-step logging.
    python tests/analysis_first_shift.py <copy of darwin_long1 checkpoint.pt> <out.json> [cpu run wall seconds]
Parts:
  A  early trajectory (0-40 years) from the histories of darwin_long1 and darwin_mig1; onset of the trait-response rule
  B  vesicle-level selection at the plateau (darwin_long1, 2009 years): chains per vesicle, between/within-vesicle
     variance of chain z, SD of vesicle mean z, selection differential per step from washout x exp(-s z_mean)
  C  short CPU run (Darwin demo, 16 sites x 512 slots, 64 vesicle slots, L_MAX 1024, seed 1) logged every step:
     mother/daughter and sister z_mean correlations at division, inflow of chains of new ancestry (initiations)
     into vesicles, loss composition of lineage chains, copy limits, saturation, division and slot fill, and the same
     vesicle statistics as B over the early years
  D  mutation load (mean dz of one random domain mutation at z ~ 0.4 and ~ 2) and the landscape ceiling (exhaustive
     over all 2^24 domains)"""
import json, math, os, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
ENG = os.path.join(HERE, "..", "..", "gokyona_env", "engine")
sys.path.insert(0, ENG)
import gk_chem, words as wd, darwin  # noqa: E402
from trait import Trait  # noqa: E402

torch.set_num_threads(4)
CK, OUT = sys.argv[1], sys.argv[2]
WALL = float(sys.argv[3]) if len(sys.argv) > 3 else 420
RUNS = os.path.expanduser("~/.gokyona-env/runs")
PRE = json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
YEAR = 8766.0
res = {"created": time.strftime("%Y-%m-%d %H:%M:%S"), "note": "CPU only; darwin_mig1 not touched; checkpoint read from a copy"}


def stats(a):
    a = np.asarray(a, float)
    return {"n": int(a.size), "mean": float(a.mean()), "sd": float(a.std()), "p10": float(np.percentile(a, 10)),
            "p50": float(np.median(a)), "p90": float(np.percentile(a, 90)), "max": float(a.max())} if a.size else None


# ---------------- A: early trajectory ----------------
def part_a(run):
    h = [json.loads(l) for l in open(os.path.join(RUNS, run, "history.jsonl"))]
    h = [x for x in h if x.get("score_mean_real") is not None]
    zc, rn, rows, onset, lamp_at = [], 0, [], None, None
    for i, x in enumerate(h):
        zc.append(x["score_mean_ctl"])
        d = x["score_mean_real"] - x["score_mean_ctl"]
        sd = float(np.std(zc[-20:], ddof=1)) if len(zc) >= 5 else None
        passing = sd is not None and d > 3 * sd
        if lamp_at is None:
            rn = rn + 1 if passing else 0
            if rn == 5:
                lamp_at, onset = i, i - 4
        yr = x["sim_h"] / YEAR
        if yr <= 40:
            rows.append({"years": round(yr, 2), "z_real": round(x["score_mean_real"], 3), "z_control": round(x["score_mean_ctl"], 3),
                         "diff": round(d, 3), "ctl_sd": None if sd is None else round(sd, 3), "passes_3sd": passing,
                         "lamp": lamp_at is not None and i >= lamp_at, "chains": x["n_chains"], "mean_len": round(x["mean_len"], 1),
                         "max_len": x["max_len"], "vesicles_per_site": round(x["vesicles"], 2), "distinct": x["distinct"],
                         "gen_max": x["gen_max"], "rho_wz_real": x.get("rho_wz_real"), "null99_real": x.get("null99_real")})
    return {"reports_0_40y": rows,
            "rule_onset_years": None if onset is None else round(h[onset]["sim_h"] / YEAR, 2),
            "rule_lamp_years": None if lamp_at is None else round(h[lamp_at]["sim_h"] / YEAR, 2),
            "first_report_years": round(h[0]["sim_h"] / YEAR, 2), "first_report_diff": round(h[0]["score_mean_real"] - h[0]["score_mean_ctl"], 3),
            "note": ("the rule needs 5 control reports for its SD, then 5 passing reports; nucleotide share, enclosed share and "
                     "copies per chain are not in the history (see part C for these over the early years)")}


res["A_early_trajectory"] = {r: part_a(r) for r in ("darwin_long1", "darwin_mig1")}
print("A done", {r: (v["first_report_years"], v["first_report_diff"], v["rule_onset_years"], v["rule_lamp_years"]) for r, v in res["A_early_trajectory"].items()}, flush=True)


# ---------------- vesicle statistics (shared by B and C) ----------------
def ves_stats(L, seq, chem, lab, alive, wmult, w_row, dt_row, N, trait, NV, half="real"):
    """Real-half (or control-half) vesicle statistics with the REAL landscape z of each scored enclosed chain."""
    R = L.shape[0]
    rows = torch.arange(R)[:, None].expand_as(L)
    hm = ((rows % N) >= N // 2) if half == "control" else ((rows % N) < N // 2)
    sc = (chem == 0) & (L >= 12) & (lab > 0) & hm
    rr, cc = sc.nonzero(as_tuple=True)
    z = trait.zscore(L[rr, cc], seq[rr, cc])
    comp = rr * NV + lab[rr, cc].long() - 1
    n = torch.zeros(R * NV).index_add_(0, comp, torch.ones_like(z))
    s1 = torch.zeros(R * NV).index_add_(0, comp, z.float())
    has = n > 0
    zbar = torch.where(has, s1 / n.clamp(min=1), torch.zeros_like(s1))
    within = float(((z - zbar[comp]) ** 2).sum())
    zm = float(z.mean())
    between = float((n * (zbar - zm) ** 2)[has].sum())
    al = alive.reshape(-1).bool() & ((torch.arange(R * NV) // NV % N < N // 2) if half == "real" else (torch.arange(R * NV) // NV % N >= N // 2))
    all_n = torch.zeros(R * NV).index_add_(0, (rows * NV + lab.long() - 1)[(lab > 0) & hm], torch.ones(int(((lab > 0) & hm).sum())))
    # washout probability per vesicle this step: 1 - exp(-w dt x mult)
    kw = (w_row * dt_row)[:, None].expand(R, NV).reshape(-1) * wmult.reshape(-1)
    pw = 1 - torch.exp(-kw)
    sel = has & al
    surv = 1 - pw[sel]
    zb, nn = zbar[sel], n[sel]
    sd_v = float(zb.std()) if zb.numel() > 1 else None
    s_ves = float((surv * zb).sum() / surv.sum() - zb.mean())                       # vesicle-weighted
    s_chain = float((surv * nn * zb).sum() / (surv * nn).sum() - (nn * zb).sum() / nn.sum())  # chain-weighted
    return {"scored_chains": int(z.numel()), "chain_z": stats(z.numpy()),
            "vesicles_alive": int(al.sum()), "vesicles_with_scored": int(sel.sum()),
            "scored_chains_per_vesicle": stats(n[sel].numpy()), "all_chains_per_vesicle": stats(all_n[al].numpy()),
            "var_total": (within + between) / max(1, z.numel()), "var_between": between / max(1, z.numel()),
            "var_within": within / max(1, z.numel()), "between_share": between / max(1e-12, within + between),
            "sd_vesicle_mean_z": sd_v, "washout_p_per_step": stats(pw[sel].numpy()),
            "selection_diff_per_step_vesicle_weighted": s_ves, "selection_diff_per_step_chain_weighted": s_chain}


# ---------------- B: plateau, darwin_long1 final checkpoint ----------------
cfg1 = json.load(open(os.path.join(RUNS, "darwin_long1", "config.json")))
ck = torch.load(CK, map_location="cpu", weights_only=False)
gk_chem.configure(int(ck["W"]), int(cfg1.get("vesicle_slots", 64)))
N1 = int(cfg1["sites"])
tr1 = Trait(dict(cfg1["trait"], enabled=True), N1, torch.device("cpu"))
VSd = ck["VS"]
dt_mean = float(ck["t"][0]) / ck["steps"]
B = ves_stats(ck["L"], ck["seq"], ck["chem"], ck["lab"], VSd["alive"], VSd["wash_mult"], VSd["w"], torch.full((ck["L"].shape[0],), dt_mean), N1, tr1, gk_chem.NV)
B["control_half_real_landscape"] = ves_stats(ck["L"], ck["seq"], ck["chem"], ck["lab"], VSd["alive"], torch.ones_like(VSd["wash_mult"]), VSd["w"],
                                             torch.full((ck["L"].shape[0],), dt_mean), N1, tr1, gk_chem.NV, "control")
B.update({"sim_years": float(ck["t"][0]) / YEAR, "steps": ck["steps"], "dt_mean_h": dt_mean, "tau_h": float(1 / VSd["w"][0]),
          "s": tr1.s, "vesicle_generations_mean": float(VSd["vgen"][VSd["alive"]].float().mean())})
B["steps_per_vesicle_generation"] = ck["steps"] / B["vesicle_generations_mean"]
B["selection_diff_per_vesicle_generation_vesicle_weighted"] = B["selection_diff_per_step_vesicle_weighted"] * B["steps_per_vesicle_generation"]
occ = (ck["L"] > 0) & (ck["chem"] == 0) & (ck["L"] >= 12)
rr, cc = occ.nonzero(as_tuple=True)
zall = tr1.zscore(ck["L"][rr, cc], ck["seq"][rr, cc])
realm = (rr % N1) < N1 // 2
B["best_chain_z_real_half"] = float(zall[realm].max())
B["best_chain_z_all"] = float(zall.max())
B["control_half_note"] = "control half scored with the real landscape (its selection used a re-shuffled table every step), washout mult set to 1"
res["B_plateau_long1_2009y"] = B
print("B done", {k: B[k] for k in ("between_share", "sd_vesicle_mean_z", "selection_diff_per_step_vesicle_weighted", "selection_diff_per_step_chain_weighted", "steps_per_vesicle_generation", "best_chain_z_real_half")},
      B["scored_chains_per_vesicle"], flush=True)


# ---------------- D: mutation load and landscape ceiling ----------------
def domain_z_all(trait):
    """Pair-score z of every 24-letter domain (letters 0..23 = bits 0..23), exhaustively."""
    D, k = trait.D, trait.k
    nw = D - k + 1
    tab = trait.table.cpu()
    sd = float(trait.sd[nw])
    out = torch.empty(1 << D)
    step = 1 << 20
    pos = torch.arange(nw)
    for s0 in range(0, 1 << D, step):
        x = torch.arange(s0, s0 + step, dtype=torch.int64)
        def sc(v):
            idx = (v[:, None] >> pos[None, :]) & ((1 << k) - 1)
            vals = tab[pos[None, :].expand_as(idx), idx] if tab.dim() == 2 else tab[idx]
            return vals.mean(1) / sd
        out[s0:s0 + step] = 0.5 * (sc(x) + sc((~x) & ((1 << D) - 1)))
    return out


zdom = domain_z_all(tr1)
# check against trait.zscore on random chains (length 200)
g = torch.Generator().manual_seed(7)
Lc = torch.full((64,), 200, dtype=torch.int64)
sq = wd.pack(torch.randint(0, 2, (64, 200), generator=g)) if hasattr(wd, "pack") else None
chk = None
if sq is not None:
    zt = tr1.zscore(Lc, sq)
    dom = (sq[:, 0] & ((1 << 24) - 1))
    chk = float((zt - zdom[dom]).abs().max())
zmax = float(zdom.max()); imax = int(zdom.argmax())
D = {"landscape_max_z": zmax, "landscape_argmax_domain": format(imax, "024b")[::-1], "landscape_min_z": float(zdom.min()),
     "share_of_domains_above": {str(t): float((zdom > t).float().mean()) for t in (0.4, 1, 2, 3)},
     "check_vs_trait_zscore_maxabs": chk, "best_chain_present_long1_real": B["best_chain_z_real_half"]}
# mutation load: all 24 single-letter mutations of domains near z 0.4 and near z 2
flip = torch.tensor([1 << i for i in range(24)], dtype=torch.int64)
for target in (0.4, 2.0):
    cand = ((zdom - target).abs() < 0.05).nonzero().squeeze(1)
    pick = cand[torch.randperm(cand.numel(), generator=g)[:5000]]
    dz = zdom[pick[:, None] ^ flip[None, :]] - zdom[pick][:, None]
    D[f"single_mutation_dz_at_z{target}"] = {"n_domains": int(pick.numel()), "mean_dz": float(dz.mean()), "sd_dz": float(dz.std()),
                                             "share_beneficial": float((dz > 0).float().mean()), "share_dz_gt_0p1": float((dz > 0.1).float().mean())}
eps_dom = 0.03107872039553474                    # domain mutations per copy at eps = 0.00129 (population snapshot)
for target in (0.4, 2.0):
    D[f"expected_dz_per_copy_at_z{target}"] = eps_dom * D[f"single_mutation_dz_at_z{target}"]["mean_dz"]
# hill climbing from 200 random domains (greedy best single flip) for comparison with the exhaustive maximum
starts = torch.randint(0, 1 << 24, (200,), generator=g)
x = starts.clone()
for _ in range(40):
    nb = x[:, None] ^ flip[None, :]
    best = zdom[nb].max(1)
    up = best.values > zdom[x]
    x = torch.where(up, nb[torch.arange(len(x)), best.indices], x)
D["hill_climb_200_starts"] = {"local_optima_z": stats(zdom[x].numpy()), "share_reaching_global_max": float((x == imax).float().mean())}
res["D_landscape_and_mutation"] = D
print("D done", json.dumps({k: v for k, v in D.items() if k != "landscape_argmax_domain"})[:900], flush=True)
json.dump(res, open(OUT, "w"), indent=1)


# ---------------- C: short CPU run with per-step logging ----------------
NC, SC = 16, 512
gk_chem.configure(16, 64)
darwin.install_ancestry(gk_chem)
dev = torch.device("cpu")
u = torch.tensor([PRE["u"]], dtype=torch.float32)
ov = {k: float(v) for k, v in PRE["env"].items() if k != "avail_time"}
ov["mineral_type"] = float(round(ov.get("mineral_type", 0)))
ov["endless"] = True
trC = Trait(dict(PRE["trait"], enabled=True), NC, dev)
trS = Trait(dict(PRE["trait"], enabled=True), NC, dev)      # scoring copy (does not touch the run's shuffle RNG)
NV = gk_chem.NV


class Log:
    def __init__(self):
        self.prev = None; self.k = 0; self.t0 = time.time()
        self.md, self.sis = [], []
        self.rows = []; self.acc = None; self.last_year = -1
        self.inflow, self.encl = [], []
        self.cum0 = None

    def due(self):
        return True

    def __call__(self, loc):
        self.k += 1
        VS = loc["VS"]; L, seq, chem, lab = loc["L"], loc["seq"], loc["chem"], loc["lab"]
        R = L.shape[0]
        t = float(loc["t"][0]); dt = float(loc["dt"][0]) if torch.is_tensor(loc["dt"]) else float(loc["dt"])
        real_rows = (torch.arange(R) % NC) < NC // 2
        # vesicle mean z (real landscape) of all alive vesicles in the real half
        sc = (chem == 0) & (L >= 12) & (lab > 0) & real_rows[:, None]
        rr, cc = sc.nonzero(as_tuple=True)
        z = trS.zscore(L[rr, cc], seq[rr, cc]) if rr.numel() else torch.zeros(0)
        comp = rr * NV + lab[rr, cc].long() - 1
        n = torch.zeros(R * NV).index_add_(0, comp, torch.ones_like(z))
        s1 = torch.zeros(R * NV).index_add_(0, comp, z.float())
        zb = (s1 / n.clamp(min=1)).view(R, NV)
        has = (n > 0).view(R, NV) & VS.alive & real_rows[:, None]
        cur = {}
        hr, hc = has.nonzero(as_tuple=True)
        vids, vpars = VS.vid[hr, hc].tolist(), VS.vpar[hr, hc].tolist()
        for r_, v_, p_, zz in zip(hr.tolist(), vids, vpars, zb[hr, hc].tolist()):
            cur[(r_, v_)] = (zz, p_)
        if self.prev is not None:
            fam = {}
            for (r_, v_), (zz, p_) in cur.items():
                if (r_, v_) not in self.prev and p_ >= 0 and (r_, p_) in self.prev:
                    self.md.append((self.prev[(r_, p_)][0], zz))
                    fam.setdefault((r_, p_), []).append(zz)
            for f in fam.values():
                if len(f) == 2:
                    self.sis.append(tuple(f))
        self.prev = cur
        # inflow of chains of outside-born ancestry into vesicles: at every hook enclosed chains get tag 1, outside chains
        # keep theirs; new initiations get 0 (engine), so an enclosed chain with tag 0 now came in from an outside-born
        # lineage since the last step (vesicle formation or capture); ftag is write-only for the engine
        ft = loc["ftag"]
        lin = (chem == 0) & (L >= 12)
        encl = lin & (lab > 0) & real_rows[:, None]
        self.inflow.append(int((encl & (ft == 0)).sum())); self.encl.append(int(encl.sum()))
        # scored chains (chemistry 0, >= 12 letters) -> tag 1; shorter chains keep theirs (0 = grown de novo from an
        # initiation, 1 = fragment or copy of a scored chain); so a scored chain with tag 0 now reached 12 letters by de-novo
        # growth since the last step, not by copying from or breaking off an existing scored chain
        ft.copy_(torch.where(lin, torch.ones_like(ft), torch.where(L > 0, ft, torch.zeros_like(ft))))
        cum = {k: float(loc[k].sum()) for k in ("loss_hyd", "loss_wash", "overflow", "lim_try", "lim_mono", "lim_slot", "sat_tot", "sat3", "cp_n")}
        cum["loss_ves"] = float(VS.loss_ves.sum()) if hasattr(VS, "loss_ves") else 0.0
        cum["n_div"] = float(VS.n_div.sum())
        if self.cum0 is None:
            self.cum0 = dict(cum)
        yr = t / YEAR
        if int(yr) > self.last_year:                         # one row per simulated year
            self.last_year = int(yr)
            row = {"years": round(yr, 2), "step": self.k, "dt_h": dt}
            occ = L > 0
            row["chains"] = int(occ.sum()); row["mean_len"] = float(L[occ].float().mean())
            row["nucleotide_share_letters"] = float(L[occ & (chem == 0)].sum()) / max(1.0, float(L[occ].sum()))
            row["enclosed_share_chains"] = float((occ & (lab > 0)).sum()) / max(1, int(occ.sum()))
            row["chains_ge24_share"] = float((occ & (chem == 0) & (L >= 24)).sum()) / max(1, int(occ.sum()))
            row["vesicles_per_site"] = float(VS.alive.sum()) / R
            row["vesicle_gen_mean"] = float(VS.vgen[VS.alive].float().mean()) if bool(VS.alive.any()) else 0.0
            ctl = (chem == 0) & (L >= 12) & ~real_rows[:, None]
            zr = z.mean().item() if z.numel() else None
            r2, c2 = ctl.nonzero(as_tuple=True)
            row["z_real_enclosed"] = zr
            row["z_control_all"] = float(trS.zscore(L[r2, c2], seq[r2, c2]).mean()) if r2.numel() else None
            if self.acc is not None:
                d = {k: cum[k] - self.acc[k] for k in cum}
                lt = d["loss_hyd"] + d["loss_wash"] + d["loss_ves"]
                row["year_losses_lineage_chains"] = {"bond_breaks": d["loss_hyd"], "free_washout": d["loss_wash"], "vesicle_washout": d["loss_ves"],
                                                     "share_vesicle_washout": d["loss_ves"] / lt if lt else None,
                                                     "fragments_lost_no_slot_all_chains": d["overflow"]}
                row["year_copies"] = d["cp_n"]; row["copies_per_chain_per_year"] = d["cp_n"] / max(1, row["chains"])
                row["copy_drop_monomer_share"] = d["lim_mono"] / d["lim_try"] if d["lim_try"] else None
                row["copy_drop_slot_share"] = d["lim_slot"] / d["lim_try"] if d["lim_try"] else None
                row["lambda_dt_gt3_share"] = d["sat3"] / d["sat_tot"] if d["sat_tot"] else None
                row["divisions_per_vesicle_per_year"] = d["n_div"] / max(1.0, float(VS.alive.sum()))
            self.acc = dict(cum)
            if VS.wash_mult is not None:
                row["ves"] = ves_stats(L, seq, chem, lab, VS.alive, VS.wash_mult, VS.w, torch.full((R,), dt), NC, trS, NV)
            self.rows.append(row)
            print("C", round(yr, 1), "step", self.k, "z", None if zr is None else round(zr, 3), round(time.time() - self.t0), "s", flush=True)
        if time.time() - self.t0 > WALL:
            return "stop"


lg = Log()
gC = torch.Generator(); gC.manual_seed(1)
t0 = time.time()
gk_chem.simulate(u, gC, (NC, SC, 300), override=ov, trait=trC, hook=lg)
md = np.array(lg.md); sis = np.array(lg.sis)
C = {"size": f"{NC} sites x {SC} slots, 64 vesicle slots, L_MAX 1024, Darwin demo, seed 1, no migration",
     "wall_s": round(time.time() - t0, 1), "steps": lg.k, "per_year": lg.rows}
if len(md) > 10:
    C["mother_daughter"] = {"n": len(md), "corr": float(np.corrcoef(md[:, 0], md[:, 1])[0, 1]),
                            "slope_daughter_on_mother": float(np.polyfit(md[:, 0], md[:, 1], 1)[0])}
if len(sis) > 10:
    C["sisters"] = {"n": len(sis), "corr": float(np.corrcoef(sis[:, 0], sis[:, 1])[0, 1])}
inf, enc = np.array(lg.inflow, float), np.array(lg.encl, float)
late = slice(len(inf) // 2, None)
gens = lg.rows[-1]["vesicle_gen_mean"] if lg.rows else None
spg = lg.k / gens if gens else None
C["new_ancestry_inflow"] = {"per_step_share_of_enclosed_scored_chains_second_half_of_run": float(inf[late].sum() / max(1, enc[late].sum())),
                            "steps_per_vesicle_generation": spg,
                            "per_vesicle_generation_share": float(inf[late].sum() / max(1, enc[late].sum())) * spg if spg else None,
                            "note": "enclosed scored chains (chemistry 0, >= 12 letters) that reached 12 letters by de-novo growth (initiation, no scored ancestor), as a share of the enclosed scored chains; all chains are enclosed in this setting"}
res["C_cpu_run"] = C
json.dump(res, open(OUT, "w"), indent=1)
print("C done", {k: v for k, v in C.items() if k != "per_year"}, flush=True)
