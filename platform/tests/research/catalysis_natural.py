"""Natural selection on catalytic activity? (coordinator request, 29 Oct 22:05). CPU only, no trait (trait off).
Short runs at 8 sites x 256 slots, 64 vesicle slots, L_MAX 1024, Darwin demo chemistry, random rationing:
  base  = fair rationing, p0 as in the preset (4.3e-6)
  p0_3  = fair rationing, p0 = 1e-3 (inside layer C: 1e-6 ... 1)
  Dp_3  = D' vesicle rules (K_max 2, no re-encapsulation, content division at 2 K_max, segregation), p0 = 1e-3
gk_cat3.apply_catalysis is wrapped (read-only) to record, per step, the expected catalysed events of every catalysed
triple (rate x nC x nA x nB x dt, before Poisson and thinning) and its compartment (X["lab"] of species A). Per
vesicle: catalytic activity = expected catalysed events per step inside it. Recorded: catalysed vs uncatalysed
condensation counts (engine counters), vesicles with activity > 0, corr(activity, divisions in the next interval) and
the mean activity of vesicles washed out vs surviving in the next step.
    python tests/catalysis_natural.py out.json [seconds per config] [config set: v1 | motif]"""
import json, os, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gokyona_env", "engine"))
import gk_chem, gk_cat3  # noqa: E402

torch.set_num_threads(4)
OUT = sys.argv[1]
SEC = float(sys.argv[2]) if len(sys.argv) > 2 else 120
PRE = json.load(open(os.path.join(HERE, "..", "..", "gokyona_env", "data", "presets.json")))["darwin_demo"]
N, S = 8, 256
gk_chem.configure(16, 64)
NV = gk_chem.NV
REC = {}
_orig = gk_cat3.apply_catalysis


def _wrapped(L, seq, chem, Ms, rows, X, trip, r_cm, r_cc, dt, gen, aux, poisson, alpha, bits, L_max, tmult=None):
    if trip.numel() and "lab" in X:
        ri, c, a, b = trip.unbind(1)
        nx = X["L"].shape[1]
        mono = b >= nx
        bs = b.clamp(max=nx - 1)
        nC = X["count"][ri, c].float(); nA = X["count"][ri, a].float()
        chA = X["chem"][ri, a]
        nB = torch.where(mono, Ms[rows[ri], chA] / alpha, X["count"][ri, bs].float())
        lam = torch.where(mono, r_cm[ri], r_cc[ri]) * nC * nA * nB * dt[ri]
        if tmult is not None:
            lam = lam * tmult
        lab_a = X["lab"][ri, a].long()
        act = torch.zeros(L.shape[0] * (NV + 1))
        act.index_add_(0, (rows[ri] * (NV + 1) + lab_a).cpu(), lam.float().cpu())
        REC["act"] = act.view(L.shape[0], NV + 1)
        REC["n_trip"] = int(trip.shape[0])
    else:
        REC["act"] = None
    ev, sat = _orig(L, seq, chem, Ms, rows, X, trip, r_cm, r_cc, dt, gen, aux, poisson, alpha, bits, L_max, tmult=tmult)
    if REC.get("act") is not None:
        # realised events per vesicle: each site's realised events split in proportion to the expected rates
        act = REC["act"]
        er = torch.zeros(L.shape[0]); er.index_add_(0, rows.cpu(), ev.float().cpu())
        tot = act.sum(1, keepdim=True).clamp(min=1e-300)
        REC["real"] = act / tot * er[:, None]
    return ev, sat


gk_cat3.apply_catalysis = _wrapped


class Hook:
    def __init__(self):
        self.t0 = time.time(); self.k = 0; self.prev = None; self.rows = []
        self.wash, self.surv = [], []
        self.div_pairs = []; self.c0 = None
        self.surv_pairs = []; self.share_active = []; self.phis = []; self.gprev = None

    def due(self):
        return True

    def __call__(self, loc):
        self.k += 1
        VS = loc["VS"]
        act = REC.get("act")
        a_v = act[:, 1:] if act is not None else torch.zeros(VS.alive.shape)
        alive = VS.alive.cpu()
        if self.prev is not None:
            pal, pid, pact, pdiv = self.prev
            same = alive & (VS.vid.cpu() == pid) | alive & (VS.vpar.cpu() == pid)
            gone = pal & ~alive                                # washed out / broken down (division keeps the slot)
            if bool(pal.any()):
                self.wash += pact[gone].tolist(); self.surv += pact[pal & same].tolist()
                self.surv_pairs += list(zip(pact[pal].tolist(), (~gone[pal]).float().tolist()))
            if getattr(VS, "vdiv", None) is not None:
                dv = (VS.vdiv.cpu() - pdiv)[pal & same]
                self.div_pairs += list(zip(pact[pal & same].tolist(), dv.tolist()))
        self.prev = (alive.clone(), VS.vid.cpu().clone(), a_v.clone(), VS.vdiv.cpu().clone() if getattr(VS, "vdiv", None) is not None else None)
        # catalysed share inside active vesicles: expected catalysed events vs uncatalysed growth of the vesicle's
        # chains (site growth this step split over the site's chains in proportion to chain number)
        g = loc["grow_row"].cpu().double()
        if self.gprev is not None and act is not None and float(loc["t"][0]) / 8766 > 10:
            dg = g - self.gprev
            L_, lab_ = loc["L"].cpu(), loc["lab"].cpu()
            nrow = (L_ > 0).sum(1).double().clamp(min=1)
            occ = L_ > 0
            rr_, cc_ = occ.nonzero(as_tuple=True)
            cnt = torch.zeros(L_.shape[0] * (NV + 1)).index_add_(0, rr_ * (NV + 1) + lab_[rr_, cc_].long(), torch.ones(rr_.numel())).view(L_.shape[0], NV + 1)
            cond_v = (dg / nrow)[:, None] * cnt.double()
            a1 = REC["real"][:, 1:].double() if REC.get("real") is not None else act[:, 1:].double(); c1 = cond_v[:, 1:]
            m = (a1 > 0) & alive
            if bool(m.any()):
                self.share_active.append(float((a1[m] / (a1[m] + c1[m]).clamp(min=1e-300)).mean()))
        self.gprev = g.clone()
        if self.k % 25 == 0:
            ph = gk_chem._ves_phi(VS).cpu()[alive]
            if ph.numel():
                self.phis.append([float(ph.min()), float(ph.quantile(0.1)), float(ph.median()), float(ph.quantile(0.9)), float(ph.max())])
            ce, co = float(loc["cat_events"].sum()), float(loc["cond_events"].sum())
            self.rows.append({"years": round(float(loc["t"][0]) / 8766, 2), "cat_events_cum": ce, "cond_events_cum": co,
                              "vesicles_with_activity_share": float(((a_v > 0) & alive).sum()) / max(1, int(alive.sum())),
                              "n_triples": REC.get("n_trip", 0), "vesicles_per_site": float(alive.sum()) / alive.shape[0],
                              "chains": int((loc["L"] > 0).sum()), "mean_len": float(loc["L"][loc["L"] > 0].float().mean())})
        if time.time() - self.t0 > SEC:
            return "stop"


res = {"size": f"{N} sites x {S} slots, 64 vesicle slots, L_MAX 1024, Darwin demo, trait off, seed 1", "configs": {}}
CONFIGS = {
    "v1": (("base", {}), ("p0_3", {"p0": 1e-3}),
           ("Dp_3", {"p0": 1e-3, "k_max": 2, "no_reencapsulation": True, "division_mode": "content", "segregate": True, "division_n0": "kmax"})),
    "motif": (("motif_p3", {"p0": 1e-3, "cat_motif": {"k": 2, "nx": 64}}),
              ("motif_p3_small", {"p0": 1e-3, "cat_motif": {"k": 2, "nx": 64}, "phi0_mult": 0.25, "period": 250.0}),
              ("motif_p3_Dp", {"p0": 1e-3, "cat_motif": {"k": 2, "nx": 64}, "k_max": 2, "no_reencapsulation": True, "division_mode": "content",
                               "segregate": True, "division_n0": "kmax"}),
              ("motif_p3_Dp_rehash", {"p0": 1e-3, "cat_motif": {"k": 2, "nx": 64}, "k_max": 2, "no_reencapsulation": True, "division_mode": "content",
                                      "segregate": True, "division_n0": "kmax", "salt_per_step": True})),
    "local": (("L_p2", {"p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                         "division_mode": "content", "local_conc": True}),
              ("L_p3", {"p0": 1e-3, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                         "division_mode": "content", "local_conc": True}),
              ("L_p2_rehash", {"p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                                "division_mode": "content", "local_conc": True, "salt_per_step": True})),
    "localmin": (("LM_p2", {"p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                            "division_mode": "content", "local_conc": True, "vesicle_min_size": 0.5}),
                 ("LM_p3", {"p0": 1e-3, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                            "division_mode": "content", "local_conc": True, "vesicle_min_size": 0.5}),
                 ("LM_p2_rehash", {"p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                                   "division_mode": "content", "local_conc": True, "vesicle_min_size": 0.5, "salt_per_step": True})),
    "motif2": (("m7_p2_small_content", {"p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                                         "division_mode": "content"}),
               ("m7_p3_small_content", {"p0": 1e-3, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                                         "division_mode": "content"}),
               ("m7_p2_small_content_rehash", {"p0": 1e-2, "L_cat": 7.0, "kappa": 1e5, "cat_motif": {"k": 3, "nx": 64}, "phi0_mult": 0.25, "period": 250.0,
                                                "division_mode": "content", "salt_per_step": True})),
}
SET = sys.argv[3] if len(sys.argv) > 3 else "v1"
for name, extra in CONFIGS[SET]:
    extra = dict(extra)
    if "phi0_mult" in extra:
        import gk_chem as _g
        _p = _g.decode(torch.tensor([PRE["u"]]))
        extra["phi0"] = float(_p["phi0"][0]) * extra.pop("phi0_mult")
    u = torch.tensor([PRE["u"]], dtype=torch.float32)
    ov = {k: float(v) for k, v in PRE["env"].items() if k != "avail_time"}
    ov.update({"mineral_type": float(round(ov.get("mineral_type", 0))), "endless": True, "ration": "random", "vesicle_stats": True})
    ov.update(extra)
    hk = Hook(); REC.clear()
    g = torch.Generator(); g.manual_seed(1)
    gk_chem.simulate(u, g, (N, S, 300), override=ov, trait=None, hook=hk)
    w, sv = np.array(hk.wash), np.array(hk.surv)
    dp = np.array(hk.div_pairs) if hk.div_pairs else np.zeros((0, 2))
    last = hk.rows[-1] if hk.rows else {}
    r = {"years": last.get("years"), "steps": hk.k,
         "catalysed_share_of_condensations": (last["cat_events_cum"] / max(1.0, last["cat_events_cum"] + last["cond_events_cum"])) if last else None,
         "cat_events": last.get("cat_events_cum"), "cond_events": last.get("cond_events_cum"),
         "vesicles_with_activity_share_late": float(np.mean([x["vesicles_with_activity_share"] for x in hk.rows[len(hk.rows) // 2:]])) if hk.rows else None,
         "triples_last": last.get("n_triples"),
         "activity_washed_mean": float(w.mean()) if w.size else None, "activity_surviving_mean": float(sv.mean()) if sv.size else None,
         "n_washed": int(w.size), "n_surviving": int(sv.size),
         "corr_activity_survival_next_step": (lambda q: float(np.corrcoef(q[:, 0], q[:, 1])[0, 1]) if len(q) > 10 and q[:, 0].std() > 0 and q[:, 1].std() > 0 else None)(np.array(hk.surv_pairs) if hk.surv_pairs else np.zeros((0, 2))),
         "catalysed_share_in_active_vesicles_realised": float(np.mean(hk.share_active)) if hk.share_active else None,
         "phi_alive_min_p10_med_p90_max_last": hk.phis[-1] if hk.phis else None,
         "max_factor_1_over_phi2": (1.0 / min(x[0] for x in hk.phis) ** 2) if hk.phis else None,
         "vesicles_per_site_last": last.get("vesicles_per_site"), "chains_last": last.get("chains"),
         "Ms_per_site_last": None,
         "corr_activity_divisions_next_step": float(np.corrcoef(dp[:, 0], dp[:, 1])[0, 1]) if len(dp) > 10 and dp[:, 0].std() > 0 and dp[:, 1].std() > 0 else None,
         "series": hk.rows}
    res["configs"][name] = r
    print(name, {k: v for k, v in r.items() if k != "series"}, flush=True)
json.dump(res, open(OUT, "w"), indent=1)
