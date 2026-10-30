"""Gokyona Environment engine: one open-ended run of the bottom-up chemistry (gk_chem = chem4c4/4c5 physics with
multi-word sequences), with live status, history, checkpoints and pause/resume/stop via control.json.
    python -m gokyona_env.engine.run --run-dir <data-dir>/runs/<name> [--device auto|xpu|cuda|mps|cpu]
"""
import argparse, datetime, json, math, os, pickle, sys, time
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import gk_chem  # noqa: E402
import devices  # noqa: E402
import darwin  # noqa: E402
from trait import Trait  # noqa: E402

PKG = os.path.dirname(HERE)
BANNER = ("Model, not life. Length beyond 64 is unvalidated. Paper result: without a sequence trait, "
          "selection stops after one sweep.")
TIME_NOTE = "Time scale follows the model's rate constants from the literature; it is not a calibrated Earth clock."
N_STEPS = 300                     # as the paper's screening runs: dt_target = t_end / 300


def atomic_json(obj, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: None if isinstance(o, float) and math.isnan(o) else str(o))
    os.replace(tmp, path)


class Hook:
    def __init__(self, run_dir, cfg, trait_score, dev, N, wall0=0.0, tracker=None, preset_env=None):
        self.dir, self.cfg, self.trait, self.dev, self.N = run_dir, cfg, trait_score, dev, N
        self.report_every = float(cfg.get("report_every_s", 5))
        self.ckpt_every = float(cfg.get("checkpoint_every_s", 240))
        self.t_start = time.time() - wall0
        self.last_report = 0.0
        self.last_ckpt = time.time()
        self.ctrl_path = os.path.join(run_dir, "control.json")
        self.ctrl_mtime = os.path.getmtime(self.ctrl_path) if os.path.exists(self.ctrl_path) else 0.0
        self.tr = tracker or darwin.Tracker(N)
        self.prev = None
        self.state = "running"
        self.last_steps, self.last_wall, self.last_sim = 0, time.time(), 0.0
        self.step_ms = None
        self.speed = None
        self.resume_info = None

    def due(self):
        return time.time() - self.last_report >= self.report_every

    def _control(self):
        if not os.path.exists(self.ctrl_path):
            return None
        m = os.path.getmtime(self.ctrl_path)
        if m == self.ctrl_mtime:
            return None
        self.ctrl_mtime = m
        try:
            return json.load(open(self.ctrl_path)).get("cmd")
        except (OSError, ValueError):
            return None

    def status(self, loc, st=None):
        t = float(loc["t"][0])
        p = loc["p"]
        period = float(p["period"][0])
        now = time.time()
        steps = int(loc["steps"])
        if steps > self.last_steps:
            self.step_ms = (now - self.last_wall) / (steps - self.last_steps) * 1000
            self.speed = (t - self.last_sim) / max(now - self.last_wall, 1e-9)
        self.last_steps, self.last_wall, self.last_sim = steps, now, t
        wall = now - self.t_start
        alloc, peak = devices.allocated(self.dev)
        cpn = float(loc["cp_n"][0]) if "cp_n" in loc else 0.0
        st = st or {}
        nlong = max(1, st.get("lineages", {}).get("distinct_ge_LH", 0))
        eta = None
        if self.speed and self.speed > 0:
            target = 8766.0 * (math.floor(t / 8766.0) + 1)
            eta = (target - t) / self.speed / 3600.0
        out = {"name": self.cfg["name"], "state": self.state, "error": None, "preset": self.cfg.get("preset"),
               "env": self.cfg.get("env"), "trait": self.cfg.get("trait"), "lmax": gk_chem.L_MAX, "W": gk_chem.wd.W,
               "lmax_auto_max": self.cfg.get("lmax_auto_max"), "sites": self.N, "slots": loc["S"],
               "time": {"sim_h": t, "sim_days": t / 24.0, "sim_years": t / 8766.0, "cycles": t / period if period > 0 else None,
                        "period_h": period, "wall_s": wall, "speed_sim_h_per_wall_s": self.speed, "eta_wall_h_to_1y": eta,
                        "eta_target_years": math.floor(t / 8766.0) + 1},
               "perf": {"steps": steps, "step_ms": self.step_ms, "gpu_temp_c": devices.gpu_temp(),
                        "vram_gib": None if alloc is None else alloc / 2 ** 30, "vram_peak_gib": None if peak is None else peak / 2 ** 30,
                        "vram_budget_gib": self.cfg.get("vram_budget_gb"), "device": devices.describe(self.dev),
                        "hw_note": self.cfg.get("hw_note")},
               "banner": BANNER, "time_note": TIME_NOTE}
        out.update(st)
        if "vesicles" in out:
            nch = max(1, st.get("chains", {}).get("n", 0))
            out["vesicles"]["molecular_generations"] = cpn / nch          # copies made per chain present (approximate)
            vs = loc.get("VS")
            if vs is not None:
                out["vesicles"]["divisions"] = float(getattr(vs, "n_div", torch.zeros(1)).sum())
        out.setdefault("chains", {})
        out["copies_total"] = cpn
        sat = loc.get("cp_sat")
        out["perf"]["copy_start_saturation"] = None if sat is None else float(sat[0])
        if "mig" in loc:
            out["perf"]["migration"] = {k_: v_ for k_, v_ in loc["mig"].items() if k_ != "log"}
        for k_ in ("lim_try", "lim_mono", "lim_slot"):
            if k_ in loc:
                out["perf"][k_] = float(loc[k_][0])
        out["perf"]["dt_h"] = float(loc["dt"][0]) if "dt" in loc and torch.is_tensor(loc["dt"]) else None
        return out

    def __call__(self, loc):
        self.last_report = time.time()
        sim_h = float(loc["t"][0])
        st = darwin.compute(loc, self.trait, self.tr, sim_h)
        stat = self.status(loc, st)
        atomic_json(stat, os.path.join(self.dir, "status.json"))
        d = st["darwin"]
        line = {"wall_s": stat["time"]["wall_s"], "sim_h": sim_h, "mean_len": st["chains"]["mean_len"],
                "max_len": st["chains"]["max_len"], "distinct": st["lineages"]["distinct_ge_LH"],
                "vesicles": st["vesicles"]["per_site"], "gen_max": st["vesicles"]["gen_max"], "n_chains": st["chains"]["n"],
                "score_mean": d["real"]["score_mean"], "score_max": d["real"]["score_max"], "score_dominant": d["real"]["score_dominant"],
                "rho_wz": d["real"]["rho_wz"], "score_mean_real": d["real"]["score_mean"], "score_mean_ctl": d["control"]["score_mean"],
                "score_dom_real": d["real"]["score_dominant"], "score_dom_ctl": d["control"]["score_dominant"],
                "rho_wz_real": d["real"]["rho_wz"], "rho_wz_ctl": d["control"]["rho_wz"],
                "null99_real": d["real"]["null99"], "null99_ctl": d["control"]["null99"],
                "step_ms": stat["perf"]["step_ms"], "gpu_temp_c": stat["perf"]["gpu_temp_c"]}
        for hn in ("real", "control"):
            dm = (d.get("dominance") or {}).get(hn)
            line[f"top1_share_{hn}"] = dm["share"] if dm else None
            line[f"top1_sites_{hn}"] = dm["sites"] if dm else None
            line[f"top1_key_{hn}"] = dm["key"] if dm else None
            line[f"top1_age_y_{hn}"] = dm["top5"][0]["age_years"] if dm and dm["top5"] else None
        line["dominant_changes"] = d.get("dominant_changes")
        for hn in ("real", "control"):
            cl = (d.get("clades") or {}).get(hn)
            line[f"clade_top1_share_{hn}"] = cl["top1_share"] if cl else None
            line[f"clade_top1_sites_{hn}"] = cl["top1_sites"] if cl else None
        line["clade_sweeps"] = len(d.get("clade_sweeps") or [])
        line["clade_counted"] = d.get("clade_counted")
        if "diag" in loc:                                  # cumulative diagnostics (gk_chem), real half
            line["diag"] = [round(float(x), 4) for x in loc["diag"]]
            line["lim"] = [float(loc[k].sum()) for k in ("lim_try", "lim_mono", "lim_slot")]
            line["breaks"] = [float(loc["loss_hyd"].sum()), float(loc["guard_n"].sum())]
            line["copyc"] = [float(loc[k].sum()) for k in ("sat_tot", "sat3", "cp_sat", "cp_n", "cp_len", "cp_ext")]
            line["Ms_mean"] = round(float(loc["Ms"][:, 0].float().mean()), 1)
            if loc.get("moli"):
                line["mo_n"] = float(loc["mo_n"].sum())
            if loc.get("lprim"):
                line["lp_starts"] = float(loc["lp_starts"].sum()); line["lp_expected"] = float(loc["lp_expected"].sum())
            if loc.get("r_phys") is not None:
                line["r_phys_mean"] = float(loc["r_phys"].mean())
            if loc.get("ckin"):
                line["ck_done"] = float(loc["ck_done"].sum()); line["ck_growing"] = int((loc["ctarget"] > 0).sum())
            line["grow_cum"] = float(loc["grow_row"].sum())
            sh = loc["sat_hist"].sum(0).double()           # lambda dt histogram, 0.1-decade bins from 1e-6
            if float(sh.sum()) > 0:
                cs = torch.cumsum(sh, 0); i = int(torch.searchsorted(cs, cs[-1] / 2))
                line["lamdt_median_cum"] = round(10 ** (-6 + 0.1 * (i + 0.5)), 6)
        if getattr(self.trait, "env", None):              # environment-dependent landscape: z_A, z_B per half, w(t)
            Le, ce = loc["L"], loc["chem"]
            be = (ce == 0) & (Le >= 12)
            re_, cce = be.nonzero(as_tuple=True)
            line["env_w"] = round(float(self.trait.w_of(float(loc["t"][0]))), 4)
            if re_.numel():
                za, zb = self.trait.zscore_ab(Le[re_, cce], loc["seq"][re_, cce])
                ctl_e = (re_ % self.N) >= self.N // 2
                for hn, hm in (("real", ~ctl_e), ("control", ctl_e)):
                    if bool(hm.any()):
                        line[f"zA_{hn}"] = round(float(za[hm].mean()), 4); line[f"zB_{hn}"] = round(float(zb[hm].mean()), 4)
                        line[f"pref_{hn}"] = round(float((zb[hm] - za[hm]).mean()), 4)
                        line[f"generalist_{hn}"] = round(float(torch.minimum(za[hm], zb[hm]).mean()), 4)
            if self.trait.env.get("gradient") and re_.numel():
                # per-site cline (B = 1: row = site): pref (z_B - z_A), local z_env, activated monomers, w_s
                sw = self.trait.site_w(float(loc["t"][0])).cpu()
                site_e = (re_ % self.N).cpu()
                zl = self.trait.zscore(Le[re_, cce], loc["seq"][re_, cce], rows=re_).cpu()
                prf = (zb - za).cpu()
                ps, zs = [], []
                for s_ in range(self.N):
                    mm = site_e == s_
                    ps.append(round(float(prf[mm].mean()), 4) if bool(mm.any()) else None)
                    zs.append(round(float(zl[mm].mean()), 4) if bool(mm.any()) else None)
                line["site_w"] = [round(float(x), 3) for x in sw]
                line["site_pref"] = ps
                line["site_zlocal"] = zs
                line["site_Ms"] = [round(float(x), 1) for x in loc["Ms"][:, 0].cpu()]
                hN = self.N // 2
                for hn, sl in (("real", slice(0, hN)), ("control", slice(hN, self.N))):
                    xs = [(float(sw[i]), ps[i]) for i in range(self.N)[sl] if ps[i] is not None]
                    if len(xs) >= 3:
                        xa = torch.tensor([a for a, _ in xs]); ya = torch.tensor([b for _, b in xs])
                        if float(xa.std()) > 0 and float(ya.std()) > 0:
                            line[f"cline_corr_{hn}"] = round(float(torch.corrcoef(torch.stack([xa, ya]))[0, 1]), 4)
                            line[f"cline_slope_{hn}"] = round(float(((xa - xa.mean()) * (ya - ya.mean())).sum() / ((xa - xa.mean()) ** 2).sum()), 4)
                rg = d.get("regional") or {}
                line["regional_counted"] = rg.get("counted"); line["regional_real_events"] = rg.get("real_events")
                line["regional_max_share_real"] = rg.get("max_window_share_real"); line["regional_max_share_control"] = rg.get("max_window_share_control")
                tp = rg.get("top_real_clade") or {}
                line["top_clade_centroid"] = tp.get("centroid"); line["top_clade_sites"] = tp.get("sites"); line["top_clade_share"] = tp.get("share")
            if loc.get("env_val") is not None:
                line["env_value"] = round(float(loc["env_val"]), 4)
            if loc.get("imp_rate", 0) > 0:
                line["impacts"] = int(loc["imp_n"]); line["impact_loss"] = float(loc["imp_loss"].sum())
            if "cat_n" in loc:
                line["cat_events"] = float(loc["cat_n"].sum()); line["cat_loss"] = float(loc["cat_loss"].sum())
        VSx = loc.get("VS")
        if VSx is not None and getattr(VSx, "n0", None) is not None:   # single-genotype vesicle statistics (real half)
            Lv, cv, lv = loc["L"], loc["chem"], loc["lab"]
            Rv = Lv.shape[0]; NVv = VSx.alive.shape[1]
            realr = ((torch.arange(Rv, device=Lv.device) % self.N) < self.N // 2)
            scm = (cv == 0) & (Lv >= 12) & (lv > 0) & realr[:, None]
            rv, cvv = scm.nonzero(as_tuple=True)
            alr = VSx.alive & realr[:, None]
            if rv.numel():
                zv = self.trait.zscore(Lv[rv, cvv], loc["seq"][rv, cvv])
                comp = rv * NVv + lv[rv, cvv].long() - 1
                nv = torch.zeros(Rv * NVv, device=Lv.device).index_add_(0, comp, torch.ones_like(zv))
                s1 = torch.zeros(Rv * NVv, device=Lv.device).index_add_(0, comp, zv)
                zb = s1 / nv.clamp(min=1)
                within = float(((zv - zb[comp]) ** 2).sum()); tot = float(((zv - zv.mean()) ** 2).sum())
                line["ves_within_share"] = round(within / tot, 4) if tot > 0 else None
                has = (nv > 0).view(Rv, NVv) & alr
                line["ves_zbar_sd"] = round(float(zb.view(Rv, NVv)[has].std()), 4) if int(has.sum()) > 1 else None
                na = nv.view(Rv, NVv)[alr]
                line["ves_occupied_share"] = round(float((na > 0).float().mean()), 4) if na.numel() else None
                line["ves_scored_mean"] = round(float(na.mean()), 3) if na.numel() else None
                line["ves_chain_share_1_2"] = round(float(na[(na >= 1) & (na <= 2)].sum()) / max(1.0, float(na.sum())), 4) if na.numel() else None
                line["ves_dist"] = [round(float(x), 4) for x in ((na == 0).float().mean(), (na == 1).float().mean(), (na == 2).float().mean(),
                                                                   ((na >= 3) & (na <= 4)).float().mean(), (na >= 5).float().mean())] if na.numel() else None
                # divisions since the last report vs z-bar at the last report; mother-daughter z-bar
                zb2 = zb.view(Rv, NVv)
                cur = {}
                ar, ac = has.nonzero(as_tuple=True)
                for a1, a2, zz, vi, vp in zip(ar.tolist(), ac.tolist(), zb2[ar, ac].tolist(), VSx.vid[ar, ac].tolist(), VSx.vpar[ar, ac].tolist()):
                    cur[(a1, vi)] = (zz, vp)
                prev = getattr(self, "_vprev", None)
                if prev is not None:
                    md = [(prev[(k[0], v[1])][0], v[0]) for k, v in cur.items() if k not in prev and v[1] >= 0 and (k[0], v[1]) in prev]
                    if len(md) > 10:
                        line["mother_daughter_corr"] = round(float(np.corrcoef(np.array(md).T)[0, 1]), 4)
                        line["mother_daughter_n"] = len(md)
                    dv = (VSx.vdiv - self._vdiv_prev)
                    # same vesicle line at a slot: the slot's current id equals, or descends from, its id at the last
                    # report (the first daughter keeps the mother's slot with a new id; vdiv counts at that slot)
                    ps_ = getattr(self, "_vprev_slot", {})
                    pr_ = []
                    for a1, a2, vi, vp in zip(ar.tolist(), ac.tolist(), VSx.vid[ar, ac].tolist(), VSx.vpar[ar, ac].tolist()):
                        q = ps_.get((a1, a2))
                        if q is not None and (vi == q[1] or vp == q[1] or float(dv[a1, a2]) > 0):
                            pr_.append((q[0], float(dv[a1, a2])))
                    if len(pr_) > 10 and np.std([p_[1] for p_ in pr_]) > 0:
                        line["div_vs_zbar_corr"] = round(float(np.corrcoef(np.array(pr_).T)[0, 1]), 4)
                    yrs_ = (float(loc["t"][0]) - self._t_prev) / 8766.0
                    if yrs_ > 0 and int(alr.sum()):
                        line["div_per_vesicle_year"] = round(float(dv[alr].sum()) / int(alr.sum()) / yrs_, 3)
                self._vprev, self._vdiv_prev, self._t_prev = cur, VSx.vdiv.clone(), float(loc["t"][0])
                self._vprev_slot = {(a1, a2): (zz, vi) for a1, a2, zz, vi in zip(ar.tolist(), ac.tolist(), zb2[ar, ac].tolist(), VSx.vid[ar, ac].tolist())}
            line["seg_counts"] = [float(x) for x in VSx.seg_counts]
        if "diag" in loc:                                  # scored share (protected domain) and mean chain z, all vs scored
            Lp, cp = loc["L"], loc["chem"]
            pp = loc["prom"] if loc.get("prom") is not None else torch.ones_like(Lp, dtype=torch.bool)
            ctlr = ((torch.arange(Lp.shape[0], device=Lp.device) % self.N) >= self.N // 2)[:, None]
            base = (cp == 0) & (Lp >= 12)
            rr_, cc_ = base.nonzero(as_tuple=True)
            zz = self.trait.zscore(Lp[rr_, cc_], loc["seq"][rr_, cc_]) if rr_.numel() else None
            for hn, hm in (("real", ~ctlr), ("control", ctlr)):
                b = base & hm
                line[f"scored_share_{hn}"] = round(float((b & pp).sum()) / max(1, int(b.sum())), 4)
                if zz is not None:                         # mean chain z (real landscape, chains >= 12), all vs scored
                    m_all = hm.expand_as(base)[rr_, cc_]
                    m_sc = m_all & pp[rr_, cc_]
                    line[f"z_all_{hn}"] = round(float(zz[m_all].mean()), 4) if bool(m_all.any()) else None
                    line[f"z_scored_{hn}"] = round(float(zz[m_sc].mean()), 4) if bool(m_sc.any()) else None
        mg = stat["perf"].get("migration")
        if mg:
            line["mig_moves"] = mg.get("moves")
        with open(os.path.join(self.dir, "history.jsonl"), "a") as fh:
            fh.write(json.dumps(line, default=lambda o: None) + "\n")
        cmd = self._control()
        if time.time() - self.last_ckpt >= self.ckpt_every or cmd in ("stop", "pause"):
            self.checkpoint(loc)
        if cmd == "pause":
            self.state = "paused"
            atomic_json(dict(stat, state="paused"), os.path.join(self.dir, "status.json"))
            while True:
                time.sleep(1.0)
                c = self._control()
                if c == "resume":
                    self.state = "running"
                    self.last_wall = time.time()
                    break
                if c == "stop":
                    cmd = "stop"
                    break
        if cmd == "stop":
            self.state = "stopped"
            atomic_json(dict(stat, state="stopped"), os.path.join(self.dir, "status.json"))
            return "stop"
        return None

    def checkpoint(self, loc):
        cpu = lambda x: x.detach().cpu() if torch.is_tensor(x) else x
        state = {k: cpu(loc[k]) for k in ("L", "seq", "chem", "pair", "busy", "dlen", "F", "M", "Ms", "t", "t_first_e1", "cp_n")}
        state["ftag"], state["census_k"] = cpu(loc["ftag"]), int(loc.get("census_k", -1))
        if loc.get("prom") is not None:
            state["prom"] = cpu(loc["prom"])
        state["steps"] = int(loc["steps"])
        if loc.get("ves_on"):
            state["lab"] = cpu(loc["lab"])
            VS = loc["VS"]
            state["VS"] = {k: cpu(v) for k, v in VS.__dict__.items() if torch.is_tensor(v) and not k.startswith("_")}
        state["gen"] = loc["gen"].get_state()
        state["aux"] = loc["aux"].get_state()
        state["wall_s"] = time.time() - self.t_start
        state["tracker"] = pickle.dumps(self.tr)
        state["W"] = gk_chem.wd.W
        if "mig" in loc:
            atomic_json({k: v for k, v in loc["mig"].items()}, os.path.join(self.dir, "migration_log.json"))
        atomic_json({"dom_changes": self.tr.dom_changes, "first_seen_n": len(self.tr.first_seen)},
                    os.path.join(self.dir, "lineage_turnover.json"))
        tmp = os.path.join(self.dir, "checkpoint.pt.tmp")
        torch.save(state, tmp)
        os.replace(tmp, os.path.join(self.dir, "checkpoint.pt"))
        self.last_ckpt = time.time()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--device", default="auto")
    ap.add_argument("--max-wall-s", type=float, default=0.0, help="stop after this wall time (0 = until stopped)")
    a = ap.parse_args()
    run_dir = os.path.abspath(a.run_dir)
    cfg = json.load(open(os.path.join(run_dir, "config.json")))
    status_path = os.path.join(run_dir, "status.json")
    atomic_json({"name": cfg["name"], "state": "starting", "banner": BANNER, "time_note": TIME_NOTE}, status_path)
    try:
        dev = devices.pick(a.device)
        torch.set_num_threads(max(1, min(4, os.cpu_count() or 1)))
        lmax = int(cfg.get("lmax", 64))
        W = max(1, math.ceil(lmax / 64))
        gk_chem.configure(W, int(cfg.get("vesicle_slots", 8)))
        presets = json.load(open(os.path.join(PKG, "data", "presets.json")))
        pre = presets[cfg.get("preset", "engine12538")]
        u = torch.tensor([pre["u"]], dtype=torch.float32, device=dev)
        env = dict(pre["env"]); env.update(cfg.get("env") or {})
        ov = {k: float(v) for k, v in env.items() if k != "avail_time"}
        ov["mineral_type"] = float(round(ov.get("mineral_type", 0)))
        ov["endless"] = True
        if cfg.get("multi_start"):
            ov["multi_start"] = True
        if float(cfg.get("migration_m", 0) or 0) > 0:
            ov["migration"] = float(cfg["migration_m"])
            ov["migration_mode"] = cfg.get("migration_mode", "random")
        if cfg.get("half_pool_mix"):
            ov["half_pool_mix"] = True
        if cfg.get("local_pool_mix"):
            ov["local_pool_mix"] = True
        ov["census_h"] = float(cfg.get("census_years", 100)) * 8766.0
        # monomer rationing: new runs get "random" (fair; the server and page write it explicitly); a config without the
        # key (runs created before 29 Oct 2026 evening) keeps the legacy chem4c4 slot order so that it resumes unchanged
        ov["ration"] = cfg.get("ration", "slot")
        if cfg.get("protected_domain"):
            ov["protected_domain"] = True
        _tc = cfg.get("trait") or {}
        envc = (_tc.get("env_landscape") or {}) if _tc.get("enabled") else (cfg.get("env_chem") or {})
        if envc.get("couple_chemistry") or cfg.get("env_chem"):   # the environment variable acts on the chemistry
            src = dict(envc if envc.get("couple_chemistry") else cfg["env_chem"])
            var = src.get("variable", "T")
            v0 = float(src.get("v0", ov.get(var)))
            ov["env_chem"] = {k: src[k] for k in ("schedule", "period_years") if k in src}
            ov["env_chem"].update({"variable": var, "v0": v0, "v1": float(src["v1"])})
        icf = cfg.get("impact") or {}
        if float(icf.get("rate_per_year", 0) or 0) > 0:      # meteor / comet impacts (global)
            ov["impact"] = dict(icf, rate_per_hour=float(icf["rate_per_year"]) / 8766.0)
            if (cfg.get("trait") or {}).get("env_landscape") is not None:
                cfg["trait"]["env_landscape"] = dict(cfg["trait"]["env_landscape"], impact=icf)
        if cfg.get("uv_gradient"):                         # per-site UV on the latitude axis
            ov["uv_gradient"] = {"uv_pole": float(cfg["uv_gradient"].get("uv_pole", ov.get("UV"))),
                                 "uv_eq": float(cfg["uv_gradient"]["uv_eq"])}
        for k_ in ("k_max", "no_reencapsulation", "division_mode", "division_membrane_req", "segregate", "vesicle_stats", "division_n0", "division_trait_s", "access_trait_s", "cat_motif", "salt_per_step", "local_conc", "vesicle_min_size", "growth_uncapped", "copy_kinetic", "monomer_competition", "copy_physical", "long_primers", "growth_physical", "mineral_oligo"):
            if cfg.get(k_) not in (None, False, 0, "external"):
                ov[k_] = cfg[k_]
        cat = cfg.get("catastrophe") or {}
        if float(cat.get("rate_per_year", 0) or 0) > 0:      # per-site catastrophes (events per simulated year)
            ov["catastrophe_rate"] = float(cat["rate_per_year"]) / 8766.0
            ov["catastrophe_kill"] = float(cat.get("kill", 0.95))
        if cfg.get("promoter_required"):
            ov["promoter_required"] = True
        if cfg.get("domain_break_guard"):
            ov["domain_break_guard"] = int(cfg["domain_break_guard"]) if not isinstance(cfg["domain_break_guard"], bool) else 24
        if cfg.get("dt_max_h"):
            ov["dt_max"] = float(cfg["dt_max_h"])
        N, S = int(cfg.get("sites", 64)), int(cfg.get("slots", 512))
        tcfg = cfg.get("trait") or {}
        trait_dyn = Trait(tcfg, N, dev) if tcfg.get("enabled") else None
        trait_score = trait_dyn or Trait(dict(tcfg, enabled=False), N, dev)
        guard = None
        try:                                           # optional GPU thermal guard (Linux xe hwmon sysfs); skipped elsewhere
            import thermal
            if thermal.find_hwmon() is not None and dev.type != "cpu":
                thermal.gpu_lock(path=os.path.join(os.path.dirname(os.path.dirname(run_dir)), ".gpu.lock"))
                guard = thermal.Guard(log_path=os.path.join(run_dir, "gpu_temp.log"), job="gokyona " + cfg["name"], lock=False)
        except Exception as e:  # noqa: BLE001
            print("thermal guard not used:", repr(e)[:200], flush=True)
        resume, tracker, wall0 = None, None, 0.0
        cp = os.path.join(run_dir, "checkpoint.pt")
        if os.path.exists(cp):
            resume = torch.load(cp, weights_only=False)
            if resume.get("W") != W:
                raise RuntimeError(f"checkpoint has W={resume.get('W')}, config asks W={W}")
            tracker = pickle.loads(resume["tracker"])
            wall0 = resume.get("wall_s", 0.0)
            print("resuming at sim_h", float(resume["t"][0]), flush=True)
        darwin.install_ancestry(gk_chem)
        hwp = os.path.join(os.path.dirname(os.path.dirname(run_dir)), "hardware.json")
        if os.path.exists(hwp):
            try:
                hw = json.load(open(hwp))
                cfg.setdefault("hw_note", hw.get("hw_note")); cfg.setdefault("lmax_auto_max", hw.get("lmax_max"))
            except (OSError, ValueError):
                pass
        hook = Hook(run_dir, cfg, trait_score, dev, N, wall0, tracker)
        if guard is not None or a.max_wall_s:
            base_call = hook.__call__
            t0 = time.time()

            class _H:
                def due(self_):
                    return hook.due()

                def __call__(self_, loc):
                    if guard is not None:
                        guard.check()
                    r = base_call(loc)
                    if a.max_wall_s and time.time() - t0 >= a.max_wall_s:
                        hook.checkpoint(loc)
                        atomic_json(dict(json.load(open(status_path)), state="stopped"), status_path)
                        return "stop"
                    return r
            hk = _H()
        else:
            hk = hook
        g = torch.Generator(device=dev)
        g.manual_seed(int(cfg.get("seed", 1)))
        gk_chem.simulate(u, g, (N, S, N_STEPS), override=ov, trait=trait_dyn, hook=hk, resume=resume)
        st = json.load(open(status_path))
        if st.get("state") == "running":
            st["state"] = "finished"
            atomic_json(st, status_path)
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        try:
            st = json.load(open(status_path))
        except (OSError, ValueError):
            st = {"name": cfg.get("name")}
        st.update({"state": "error", "error": repr(e)[:500], "banner": BANNER, "time_note": TIME_NOTE})
        atomic_json(st, status_path)
        sys.exit(1)


if __name__ == "__main__":
    main()
