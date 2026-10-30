"""Live statistics of Gokyona Environment: chains, lineages, vesicles, ladder-level proxies, and the Darwin panel
(real half vs control half, selective sweeps, lineage tree, rho(w, z) against a permutation null).
numpy + torch only. The live ladder lamps are proxies of the paper's batch detectors (see README)."""
import math
import numpy as np
import torch
import words as wd
import gk_cat3 as cat
import gk_chem

L_H = 6
NPERM = 199
PREFIX = 24          # lineage identity for dominance / sweeps: the first 24 letters (copies start at letter 0)
MAXID = 1 << 15      # vesicle ancestry table per site (id modulo MAXID; ample between two reports)


def install_ancestry(gk_chem):
    """Record parent ids at every division (read-only hook; draws no random numbers)."""
    orig = gk_chem.Vesicles._divide
    if getattr(orig, "_gk", False):
        return

    def _divide(self, t, div, L, lab, pair):
        old_alive = self.alive.clone()
        orig(self, t, div, L, lab, pair)
        new = div | (self.alive & ~old_alive)
        if not hasattr(self, "par_tab"):
            self.par_tab = torch.full((self.R, MAXID), -1, dtype=torch.int64, device=self.dev)
        rr, cc = new.nonzero(as_tuple=True)
        if rr.numel():
            self.par_tab[rr, self.vid[rr, cc] % MAXID] = self.vpar[rr, cc]
    _divide._gk = True
    gk_chem.Vesicles._divide = _divide


def _rank(x):
    o = np.argsort(x, kind="mergesort")
    r = np.empty(len(x)); r[o] = np.arange(len(x))
    xs = np.asarray(x)[o]
    i = 0
    while i < len(xs):                                  # average ranks for ties
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[o[i:j + 1]] = (i + j) / 2
        i = j + 1
    return r


def spearman(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    a, b = _rank(x), _rank(y)
    return float(np.corrcoef(a, b)[0, 1])


def null99(x, y, rng):
    if len(x) < 10 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    a, b = _rank(x), _rank(y)
    vals = [np.corrcoef(a, rng.permutation(b))[0, 1] for _ in range(NPERM)]
    return float(np.percentile(vals, 99))


class Tracker:
    """Host-side history across reports (picklable for checkpoints)."""

    def __init__(self, N):
        self.N = N
        self.half = {"real": {"dom": None, "dom_score": None, "dom_share": 0.0, "dom_series": []},
                     "control": {"dom": None, "dom_score": None, "dom_share": 0.0, "dom_series": []}}
        self.sweeps = []
        self.nodes = {}
        self.parent = {}            # vesicle id -> parent id (accumulated)
        self.prev_ves = None        # {(row, vid): z} at the previous report
        self.prev_dom_keys = None
        self.dom_persist = 0
        self.ctrl_jumps = []        # |delta score| of control-half dominant changes (null for sweeps)
        self.rng = np.random.default_rng(7)
        self.vstat = {"real": {}, "control": {}}
        self.zr, self.zc = [], []          # mean trait z per report, real / control half
        self.resp_run = 0                  # consecutive reports passing the current trait-response rule
        self.resp_vals = []
        self.plateaus = []                 # levels (real - control) of the accepted shifts
        self.l4hist = {"real": [], "control": []}   # rho(w, z) > null99 per report
        self.first_seen = {}               # lineage key -> first sim_h it was observed (any half)
        self.dom_changes = []              # changes of the dominant lineage identity per half


def compute(loc, trait, tr, sim_h):
    """Statistics of the current state. loc: locals() of gk_chem._simulate at a hook call."""
    L, seq, chem = loc["L"], loc["seq"], loc["chem"]
    N, S = loc["N"], loc["S"]
    R = L.shape[0]
    dev = L.device
    ves_on = loc.get("ves_on", False)
    lab = loc["lab"] if ves_on else None
    occ = L > 0
    rows = torch.arange(R, device=dev)[:, None].expand(R, S)
    ctrl = (rows % N) >= N // 2
    Lf = L.float()
    tot = float(Lf[occ].sum()) or 1.0
    lmax_now = int(L.max()) if bool(occ.any()) else 0
    edges = [1]
    while edges[-1] < max(lmax_now, 2):
        edges.append(edges[-1] * 2)
    Lo = L[occ].cpu().numpy()
    counts = np.histogram(Lo, bins=edges + [edges[-1] * 2])[0].tolist() if len(Lo) else []
    st = {"chains": {"n": int(occ.sum()), "mean_len": float(Lf[occ].mean()) if bool(occ.any()) else 0.0, "max_len": lmax_now,
                     "len_hist": {"edges": edges + [edges[-1] * 2], "counts": counts},
                     "nuc_share": float(Lf[occ & (chem == 0)].sum()) / tot,
                     "enclosed_share": (float(Lf[occ & (lab > 0)].sum()) / tot) if ves_on else 0.0}}
    # ---- lineages (chemistry 0, L >= L_H) ----
    long_ = occ & (chem == 0) & (L >= L_H)
    if loc.get("prom") is not None:                   # protected domain: chains without a start mark get no lineage
        long_ = long_ & loc["prom"]
    rr, cc = long_.nonzero(as_tuple=True)
    if rr.numel():
        Lk = L[rr, cc].clamp(max=PREFIX)
        keys = cat.species_key(chem[rr, cc], Lk, seq[rr, cc] & wd.mask(Lk))
    else:
        keys = torch.zeros(0, dtype=torch.long, device=dev)
    if rr.numel() and getattr(trait, "env", None) and trait.env.get("gradient"):
        scores = trait.zscore(L[rr, cc], seq[rr, cc], rows=rr)    # local environment of each chain's site
    else:
        scores = trait.zscore(L[rr, cc], seq[rr, cc]) if rr.numel() else torch.zeros(0, device=dev)   # real landscape (neutral marker when off)
    ishalf_c = ctrl[rr, cc]
    top = []
    dist = 0
    halves = {}
    if rr.numel():
        uk, inv, cnt = torch.unique(keys, return_inverse=True, return_counts=True)
        dist = int(uk.numel())
        order = torch.argsort(cnt, descending=True)[:8]
        first = torch.full((uk.numel(),), keys.numel(), dtype=torch.long, device=dev).scatter_reduce(
            0, inv, torch.arange(keys.numel(), device=dev), reduce="amin")
        for oi in order.tolist():
            i = int(first[oi])
            r_, c_ = int(rr[i]), int(cc[i])
            top.append({"len": int(L[r_, c_]), "count": int(cnt[oi]), "score": round(float(scores[i]), 3),
                        "seq": wd.to_str(seq[r_, c_].cpu().tolist(), int(L[r_, c_]), 128),
                        "key": "%06x" % (int(keys[i]) & 0xFFFFFF), "half": "control" if bool(ishalf_c[i]) else "real"})
        for hname, hm in (("real", ~ishalf_c), ("control", ishalf_c)):
            if not bool(hm.any()):
                halves[hname] = None
                continue
            kh, sh, lh = keys[hm], scores[hm], L[rr, cc][hm].float()
            ukh, invh, cnth = torch.unique(kh, return_inverse=True, return_counts=True)
            j = int(torch.argmax(cnth))
            members = invh == j
            sites_h = (rr[hm] % N)[members]
            halves[hname] = {"score_mean": float((sh * lh).sum() / lh.sum()), "score_max": float(sh.max()),
                             "dom_key": int(ukh[j]), "dom_share": float(cnth[j]) / float(kh.numel()),
                             "dom_score": float(sh[members][0]), "dom_len": int(lh[members][0]),
                             "dom_sites": int(torch.unique(sites_h).numel())}
            # first-seen times of all lineages; ages of the current top 5 of this half
            ukn = ukh.cpu().numpy().tolist()
            for kk in ukn:
                if kk not in tr.first_seen:
                    tr.first_seen[kk] = sim_h
            o5 = torch.argsort(cnth, descending=True)[:5].tolist()
            halves[hname]["top5"] = [{"key": "%06x" % (ukn[i] & 0xFFFFFF), "share": round(float(cnth[i]) / float(kh.numel()), 4),
                                      "first_seen_h": tr.first_seen[ukn[i]], "age_years": round((sim_h - tr.first_seen[ukn[i]]) / 8766, 2)}
                                     for i in o5]
    st["lineages"] = {"distinct_ge_LH": dist, "top": top}
    # ---- sweeps and tree ----
    for hname in ("real", "control"):
        h, H = halves.get(hname), tr.half[hname]
        if h is None:
            continue
        if H["dom"] is not None and h["dom_key"] != H["dom"] and h["dom_share"] >= 0.3:
            ds = h["dom_score"] - (H["dom_score"] or 0.0)
            if hname == "control":
                tr.ctrl_jumps.append(abs(ds))
            thr = float(np.percentile(tr.ctrl_jumps, 99)) if len(tr.ctrl_jumps) >= 5 else 0.2
            beyond = ds > thr
            tr.sweeps.append({"sim_h": round(sim_h, 2), "half": hname, "key": "%06x" % (h["dom_key"] & 0xFFFFFF),
                              "score": round(h["dom_score"], 3), "delta": round(ds, 3), "share": round(h["dom_share"], 3),
                              "beyond_null": bool(beyond), "null_thr": round(thr, 3)})
            tr.nodes["%06x" % (h["dom_key"] & 0xFFFFFF)] = {"id": "%06x" % (h["dom_key"] & 0xFFFFFF),
                                                             "parent": "%06x" % (H["dom"] & 0xFFFFFF), "born_h": round(sim_h, 2),
                                                             "score": round(h["dom_score"], 3), "share": round(h["dom_share"], 3),
                                                             "len": h["dom_len"], "half": hname}
        elif H["dom"] is None:
            tr.nodes["%06x" % (h["dom_key"] & 0xFFFFFF)] = {"id": "%06x" % (h["dom_key"] & 0xFFFFFF), "parent": None,
                                                             "born_h": round(sim_h, 2), "score": round(h["dom_score"], 3),
                                                             "share": round(h["dom_share"], 3), "len": h["dom_len"], "half": hname}
        if H["dom"] is not None and h["dom_key"] != H["dom"]:
            tr.dom_changes.append({"sim_h": round(sim_h, 1), "half": hname, "old": "%06x" % (H["dom"] & 0xFFFFFF),
                                   "new": "%06x" % (h["dom_key"] & 0xFFFFFF), "share": round(h["dom_share"], 4),
                                   "sites": h["dom_sites"], "new_first_seen_h": tr.first_seen.get(h["dom_key"])})
        H["dom"], H["dom_score"], H["dom_share"] = h["dom_key"], h["dom_score"], h["dom_share"]
    # ---- vesicles: rho(w, z) per half (z = mean score of enclosed long chains, w = descendants next report) ----
    darwin_h = {"real": {}, "control": {}}
    if ves_on:
        VS = loc["VS"]
        al = VS.alive
        st["vesicles"] = {"per_site": float(al.sum()) / R,
                          "gen_max": int(VS.vgen[al].max()) if bool(al.any()) else 0,
                          "gen_mean": float(VS.vgen[al].float().mean()) if bool(al.any()) else 0.0}
        vr, vc = al.nonzero(as_tuple=True)
        vid = VS.vid[vr, vc].cpu().numpy(); vrow = vr.cpu().numpy()
        par = VS.par_tab.cpu().numpy() if hasattr(VS, "par_tab") else None
        # z per vesicle
        if rr.numel() and lab is not None:
            lb = lab[rr, cc].long()
            inside = lb > 0
            NV1 = gk_chem.NV + 1
            comp = rr * NV1 + lb
            zsum = torch.zeros(R * NV1, device=dev).index_add_(0, comp[inside], (scores * L[rr, cc].float())[inside])
            zl = torch.zeros(R * NV1, device=dev).index_add_(0, comp[inside], L[rr, cc].float()[inside])
            z = (zsum / zl.clamp(min=1)).view(R, NV1)[:, 1:]
            zv = z[vr, vc].cpu().numpy()
        else:
            zv = np.zeros(len(vid))
        cur = {(int(rw), int(a)): float(zz) for rw, a, zz in zip(vrow, vid, zv)}
        if tr.prev_ves and par is not None:
            pk = np.array([k[0] * (1 << 40) + k[1] for k in tr.prev_ves], dtype=np.int64)
            pidx = {int(v): i for i, v in enumerate(pk)}
            wcount = np.zeros(len(pk))
            rows_a, x = vrow.astype(np.int64), vid.astype(np.int64)
            live = np.ones(len(x), bool)
            for _ in range(4000):
                if not live.any():
                    break
                code = rows_a * (1 << 40) + x
                hit = live & np.isin(code, pk)
                for c in code[hit]:
                    wcount[pidx[int(c)]] += 1
                live &= ~hit
                nxt = par[rows_a, x % MAXID]
                live &= (nxt >= 0) & (nxt != x)
                x = np.where(live, nxt, x)
            ks = list(tr.prev_ves)
            for hname, hc in (("real", False), ("control", True)):
                sel = [i for i, k in enumerate(ks) if ((k[0] % N) >= N // 2) == hc]
                zz = np.array([tr.prev_ves[ks[i]] for i in sel]); ww = wcount[sel]
                rho = spearman(ww, zz)
                thr = null99(ww, zz, tr.rng)
                darwin_h[hname] = {"rho_wz": None if math.isnan(rho) else round(rho, 3),
                                   "null99": None if math.isnan(thr) else round(thr, 3), "n_ancestors": len(sel),
                                   "L4": bool(len(sel) >= 10 and not math.isnan(rho) and not math.isnan(thr) and rho > thr)}
        tr.prev_ves = cur
    else:
        st["vesicles"] = {"per_site": 0.0, "gen_max": 0, "gen_mean": 0.0}
    for hname in ("real", "control"):
        if darwin_h[hname]:
            tr.l4hist[hname].append(bool(darwin_h[hname].get("L4")))
            last5 = tr.l4hist[hname][-5:]
            darwin_h[hname]["L4_last_report"] = darwin_h[hname]["L4"]
            darwin_h[hname]["L4"] = len(last5) == 5 and sum(last5) >= 3       # persistent: >= 3 of the last 5 reports
            darwin_h[hname]["L4_last5"] = int(sum(last5))
    for hname in ("real", "control"):
        h = halves.get(hname) or {}
        tr.vstat[hname].update({k: round(v, 4) if isinstance(v, float) else v for k, v in h.items() if k in ("score_mean", "score_max", "dom_score")})
        tr.vstat[hname].update(darwin_h[hname])
    # clades by descent (founder tags; census every census_h): top tag share among the half's lineage chains.
    # Rule (coordinator 29 Oct 11:29; fixed before any clade results were used): every tag that first reaches >= 30% of
    # its half's tagged chains is logged (both halves). A real-half event counts towards "sustained" only if
    #  (1) share >= 30%;
    #  (2) share > the 95th percentile of the control half's top-clade share over all reports so far (>= 20 control
    #      reports needed, otherwise "pending" until there are 20; then judged with the share it had at the event);
    #  (3) the clade's mean z (length-weighted, at the event) > the real half's mean z at its census (the first report
    #      after that census), i.e. a fitter clade, not drift.
    if not hasattr(tr, "clade_sweeps"):
        tr.clade_sweeps, tr.clade_seen = [], set()
    if not hasattr(tr, "clade_ctl_hist"):
        tr.clade_ctl_hist, tr.census_z = [], {}
    clades = {"real": None, "control": None}
    if "ftag" in loc and rr.numel():
        ftg = loc["ftag"][rr, cc]
        Lf = L[rr, cc].float()
        kmax = int(ftg.max()) >> 32 if bool((ftg > 0).any()) else -1
        if kmax >= 0 and kmax not in tr.census_z and halves.get("real"):
            tr.census_z[kmax] = halves["real"]["score_mean"]
        for hname, hm in (("control", ishalf_c), ("real", ~ishalf_c)):
            tg = ftg[hm & (ftg > 0)]
            nh = int(hm.sum())
            if tg.numel() == 0 or nh == 0:
                continue
            ut, ct = torch.unique(tg, return_counts=True)
            j = int(torch.argmax(ct))
            mem = hm & (ftg == ut[j])
            sites_c = int(torch.unique((rr[mem] % N)).numel())
            zc = float((scores[mem] * Lf[mem]).sum() / Lf[mem].sum())
            share = float(ct[j]) / nh
            clades[hname] = {"tag": "%x" % int(ut[j]), "census": int(ut[j]) >> 32, "top1_share": round(share, 4),
                             "top1_sites": sites_c, "n_clades": int(ut.numel()), "z_mean": round(zc, 3)}
            if hname == "control":
                tr.clade_ctl_hist.append(share)
            key = (hname, int(ut[j]))
            if share >= 0.3 and key not in tr.clade_seen:
                tr.clade_seen.add(key)
                ev = {"sim_h": round(sim_h, 1), "half": hname, "tag": "%x" % int(ut[j]), "share": round(share, 4),
                      "sites": sites_c, "z_mean": round(zc, 3)}
                if hname == "real":
                    cz = tr.census_z.get(int(ut[j]) >> 32)
                    ev.update({"census_z_real": None if cz is None else round(cz, 3), "verdict": "pending"})
                tr.clade_sweeps.append(ev)
    # judge pending real-half events once >= 20 control reports exist
    if len(tr.clade_ctl_hist) >= 20:
        p95 = float(np.percentile(tr.clade_ctl_hist, 95))
        for ev in tr.clade_sweeps:
            if ev["half"] == "real" and ev.get("verdict") == "pending":
                fitter = ev.get("census_z_real") is not None and ev["z_mean"] > ev["census_z_real"]
                ev["ctl_p95"] = round(p95, 4)
                ev["verdict"] = "counts" if (ev["share"] > p95 and fitter) else (
                    "below control p95" if ev["share"] <= p95 else "not fitter")
    # regional clade takeovers (spatial gradient; coordinator 29 Oct 20:04): in every window of 4 adjacent sites of a
    # half, the top founder tag's share of the window's tagged lineage chains; an event when a tag first reaches >= 30%
    # in some window (counted once per tag). A real-half event counts if its share exceeds the 95th percentile of the
    # control windows' top shares over all reports so far (>= 20 reports) and the clade's mean local z exceeds the
    # window's mean local z at the tag's census (first report after it); else "pending" / "below control p95" / "not fitter"
    regional = None
    if "ftag" in loc and rr.numel() and getattr(trait, "env", None) and trait.env.get("gradient"):
        if not hasattr(tr, "reg_events"):
            tr.reg_events, tr.reg_seen, tr.reg_ctl_hist, tr.reg_census_z = [], set(), [], {}
        ftg = loc["ftag"][rr, cc]
        site = rr % N
        half_n = N // 2
        kcen = int(ftg.max()) >> 32 if bool((ftg > 0).any()) else -1
        new_census = kcen >= 0 and kcen not in tr.reg_census_z
        if new_census:
            tr.reg_census_z[kcen] = {}
        ctl_tops, real_tops = [], []
        for h0, hname in ((half_n, "control"), (0, "real")):
            for w0 in range(h0, h0 + half_n - 3):
                m = (site >= w0) & (site < w0 + 4)
                if new_census and hname == "real" and bool(m.any()):
                    tr.reg_census_z[kcen][w0] = float(scores[m].mean())
                tg = ftg[m & (ftg > 0)]
                if tg.numel() < 20:
                    continue
                ut, ct = torch.unique(tg, return_counts=True)
                j = int(torch.argmax(ct)); share = float(ct[j]) / int(m.sum())
                (ctl_tops if hname == "control" else real_tops).append((share, w0, int(ut[j])))
                if share >= 0.3 and (hname, int(ut[j])) not in tr.reg_seen:
                    tr.reg_seen.add((hname, int(ut[j])))
                    mem = m & (ftg == ut[j])
                    ev = {"sim_h": round(sim_h, 1), "half": hname, "window": [w0 - h0, w0 - h0 + 3], "tag": "%x" % int(ut[j]),
                          "share": round(share, 4), "z_mean": round(float(scores[mem].mean()), 3)}
                    if hname == "real":
                        cz = tr.reg_census_z.get(int(ut[j]) >> 32, {}).get(w0)
                        ev.update({"census_z_window": None if cz is None else round(cz, 3), "verdict": "pending"})
                    tr.reg_events.append(ev)
        if ctl_tops:
            tr.reg_ctl_hist.append(max(c[0] for c in ctl_tops))
        if len(tr.reg_ctl_hist) >= 20:
            p95r = float(np.percentile(tr.reg_ctl_hist, 95))
            for ev in tr.reg_events:
                if ev["half"] == "real" and ev.get("verdict") == "pending":
                    fit = ev.get("census_z_window") is not None and ev["z_mean"] > ev["census_z_window"]
                    ev["ctl_p95"] = round(p95r, 4)
                    ev["verdict"] = "counts" if (ev["share"] > p95r and fit) else ("below control p95" if ev["share"] <= p95r else "not fitter")
        # top real clade: site range and centroid (site index along the line, 0 = w 0 pole)
        rm = (~ishalf_c) & (ftg > 0)
        top = None
        if bool(rm.any()):
            ut, ct = torch.unique(ftg[rm], return_counts=True)
            j = int(torch.argmax(ct))
            ms = site[rm & (ftg == ut[j])].float()
            top = {"tag": "%x" % int(ut[j]), "share": round(float(ct[j]) / int((~ishalf_c).sum()), 4),
                   "sites": int(torch.unique(ms).numel()), "centroid": round(float(ms.mean()), 3),
                   "range": [int(ms.min()), int(ms.max())]}
        regional = {"events": tr.reg_events[-50:], "counted": sum(1 for e in tr.reg_events if e.get("verdict") == "counts"),
                    "real_events": sum(1 for e in tr.reg_events if e["half"] == "real"),
                    "control_events": sum(1 for e in tr.reg_events if e["half"] == "control"),
                    "max_window_share_real": round(max([c[0] for c in real_tops] or [0]), 4),
                    "max_window_share_control": round(max([c[0] for c in ctl_tops] or [0]), 4),
                    "ctl_p95": round(float(np.percentile(tr.reg_ctl_hist, 95)), 4) if len(tr.reg_ctl_hist) >= 20 else None,
                    "top_real_clade": top}
    clade_sweeps = tr.clade_sweeps
    real_sweeps = [s for s in tr.sweeps if s["half"] == "real" and s["beyond_null"] and s["delta"] > 0]
    real_clade = [c for c in clade_sweeps if c["half"] == "real" and c.get("verdict") == "counts"]
    # trait response (fixed 29 Oct 08:45 before further runs): real - control mean z above 3 x the SD of the control's
    # report-to-report values (last 20 reports, >= 5) for >= 5 consecutive reports -> one shift (plateau = mean of
    # those reports); a further shift needs >= 5 consecutive reports above the previous plateau + 3 SD.
    hr, hc = halves.get("real"), halves.get("control")
    resp = {"lamp": False, "shifts": len(tr.plateaus), "plateaus": [round(x, 3) for x in tr.plateaus], "diff": None, "sd_ctl": None, "run": tr.resp_run}
    if hr and hc:
        tr.zr.append(hr["score_mean"]); tr.zc.append(hc["score_mean"])
        diff = hr["score_mean"] - hc["score_mean"]
        resp["diff"] = round(diff, 3)
        if len(tr.zc) >= 5:
            sd = float(np.std(tr.zc[-20:], ddof=1))
            resp["sd_ctl"] = round(sd, 3)
            base = tr.plateaus[-1] if tr.plateaus else 0.0
            if diff > base + 3 * sd:
                tr.resp_run += 1
                tr.resp_vals.append(diff)
            else:
                tr.resp_run, tr.resp_vals = 0, []
            if tr.resp_run >= 5:
                tr.plateaus.append(float(np.mean(tr.resp_vals)))
                tr.resp_run, tr.resp_vals = 0, []
        resp.update({"shifts": len(tr.plateaus), "plateaus": [round(x, 3) for x in tr.plateaus], "run": tr.resp_run,
                     "lamp": len(tr.plateaus) >= 1})
    if len(tr.plateaus) >= 2 or len(real_sweeps) >= 2 or len(real_clade) >= 2:
        status = "sustained evolution"
    elif tr.plateaus:
        status = "selection on the trait"
    elif real_sweeps or real_clade:
        status = "one sweep"
    else:
        status = "no selection"
    st["darwin"] = {"status": status,
                    "status_note": ("selection on the trait = real - control mean z > 3 x control SD for >= 5 consecutive reports "
                                    "(one shift); sustained evolution = >= 2 successive shifts (each a new plateau passing the same rule "
                                    "above the previous one), >= 2 successive lineage sweeps beyond the control-half null, or >= 2 counted "
                                    "real-half clade sweeps (>= 30% share, above the control's 95th-percentile top-clade share, "
                                    "fitter than the half at its census)"),
                    "trait_response": resp,
                    "dominance": {hn: ({"key": "%06x" % (halves[hn]["dom_key"] & 0xFFFFFF), "share": round(halves[hn]["dom_share"], 4),
                                        "sites": halves[hn]["dom_sites"], "top5": halves[hn]["top5"]} if halves.get(hn) else None)
                                  for hn in ("real", "control")},
                    "dominant_changes": len(tr.dom_changes), "dominant_changes_last": tr.dom_changes[-20:],
                    "clades": clades, "clade_sweeps": clade_sweeps[-50:],
                    "clade_counted": len(real_clade),
                    "regional": regional,
                    "clade_ctl_p95": round(float(np.percentile(tr.clade_ctl_hist, 95)), 4) if len(tr.clade_ctl_hist) >= 20 else None,
                    "clade_ctl_reports": len(tr.clade_ctl_hist),
                    "real": {"score_mean": tr.vstat["real"].get("score_mean"), "score_max": tr.vstat["real"].get("score_max"),
                             "score_dominant": tr.vstat["real"].get("dom_score"), "rho_wz": tr.vstat["real"].get("rho_wz"),
                             "null99": tr.vstat["real"].get("null99"), "L4": tr.vstat["real"].get("L4", False),
                             "L4_last_report": tr.vstat["real"].get("L4_last_report"), "L4_last5": tr.vstat["real"].get("L4_last5")},
                    "control": {"score_mean": tr.vstat["control"].get("score_mean"), "score_max": tr.vstat["control"].get("score_max"),
                                "score_dominant": tr.vstat["control"].get("dom_score"), "rho_wz": tr.vstat["control"].get("rho_wz"),
                                "null99": tr.vstat["control"].get("null99"), "L4": tr.vstat["control"].get("L4", False),
                                "L4_last_report": tr.vstat["control"].get("L4_last_report"), "L4_last5": tr.vstat["control"].get("L4_last5")},
                    "sweeps": tr.sweeps[-50:], "tree": {"nodes": list(tr.nodes.values())[-60:]}}
    rd = st["darwin"]["real"]
    st["trait_stats"] = {"score_mean": rd["score_mean"], "score_max": rd["score_max"], "score_dominant": rd["score_dominant"],
                         "rho_wz": rd["rho_wz"], "null99": rd["null99"]}
    # ---- ladder lamps (live proxies) ----
    e1 = not bool(torch.isnan(loc["t_first_e1"]).all())
    site_l8 = (L >= 8).any(1).float().mean().item() >= 0.5
    cat_ev = float(loc["cat_events"].sum()) if "cat_events" in loc else 0.0
    dom_keys = {h["dom_key"] for h in halves.values() if h}
    tr.dom_persist = tr.dom_persist + 1 if (tr.prev_dom_keys and dom_keys & tr.prev_dom_keys) else 0
    tr.prev_dom_keys = dom_keys
    e3 = False
    if ves_on and rr.numel():
        enc_long = torch.zeros(R, dtype=torch.bool, device=dev).index_fill_(0, rr[(lab[rr, cc] > 0)], True)
        e3 = float(enc_long.float().mean()) >= 0.5
    st["levels"] = {"E1": e1, "E2": e1 and site_l8, "L1": cat_ev > 0, "L2": tr.dom_persist >= 5, "E3": e3,
                    "L3_loose": None, "L3_strict": None,
                    "L4": bool(st["darwin"]["real"]["L4"]) and status != "no selection"}
    return st
