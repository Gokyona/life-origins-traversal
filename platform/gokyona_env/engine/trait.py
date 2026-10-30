"""Sequence -> trait rule of Gokyona Environment (OFF by default; with it off the engine is the 4c-4/4c-5 physics).

Landscapes (score(seq), with a fixed seed):
  "domain" (default, coordinator 29 Oct 07:46): a functional domain = the first D = 24 letters (copies start at
     letter 0, so copy errors reach it at rate eps per letter). score = mean of the 6-mer table over the windows
     inside the domain (19 for a full domain; chains shorter than 24 use what they have; chains < 12 score 0),
     standardised: z = score / SD_theory(n windows) for random sequences (computed once), so s is in SD units.
     "epistatic" (option): the window value depends on the 6-mer AND its position (a 19 x 64 table; NK-style,
     rugged landscape).
  "additive-mean" (reference): mean of the 64-value table over all 6-mer windows of the whole chain (unstandardised).
Modes (owner's decision, 29 Oct 2026):
  "self"   - a template's copy rate is multiplied by exp(s * score(template));
  "shared" - chains of length >= L_rep (chemistry 0) act as replicases for all templates in the same
             compartment (vesicle, or the site's outside): copy rate x (1 + s * max(0, mean score of the
             replicases present)); the benefit stays local only if enclosed.
  "membrane" - the mean z of the scored chains enclosed in a vesicle raises its membrane uptake/synthesis by
             exp(s * z_mean) (compartment-level selection; coordinator 29 Oct 08:03);
  "stability" - a chain's bond-breaking rate is multiplied by exp(-s * z) (survival).
  "division"  - a vesicle's division rate is multiplied by exp(s * z_mean) (reproduction; explicit assumption);
  "retention" - a vesicle's washout rate 1/tau is multiplied by exp(-s * z_mean) (survival; explicit assumption).
Pair score (default ON): z of a chain = mean of its own domain score and its complement's (copies are the
complementary strand, so the trait is then inherited every generation, not every second one).
Optional fidelity link: eps_eff = eps * exp(-s_f * score(template)).
Built-in control: sites [0, N/2) of each run use the real landscape; sites [N/2, N) a time-shuffled landscape
(the table is permuted afresh at every step: same score distribution, no heritable link to sequence).
"""
import torch
import words as wd


def _nv():
    import gk_chem
    return gk_chem.NV


class Trait:
    def __init__(self, cfg, N, dev):
        self.enabled = bool(cfg.get("enabled", False))
        self.mode = cfg.get("mode", "self")
        self.s = float(cfg.get("s", 0.5))
        self.k = int(cfg.get("k", 6))
        self.seed = int(cfg.get("seed", 2026))
        self.L_rep = int(cfg.get("L_rep", 24))
        self.fid = bool(cfg.get("fidelity_link", False))
        self.pair = bool(cfg.get("pair_score", True))
        self._all = None
        self.s_f = float(cfg.get("s_f", 0.0))
        self.N, self.dev = N, dev
        self.landscape = cfg.get("landscape", "domain")
        self.epistatic = bool(cfg.get("epistatic", False))
        self.D = int(cfg.get("domain_len", 24))
        g = torch.Generator().manual_seed(self.seed)
        nw = self.D - self.k + 1
        if self.landscape == "domain" and self.epistatic:
            self.table = torch.randn(nw, 1 << self.k, generator=g).to(dev)      # position-specific
        else:
            self.table = torch.randn(1 << self.k, generator=g).to(dev)
        self.sd = None
        if self.landscape == "domain":
            self.sd = self._sd_for(self.table)
        self._shuffle_gen = torch.Generator().manual_seed(self.seed + 1)
        self._step_tab = self.table
        self._cache = None
        # Gokyona environment-dependent landscape (option, coordinator 29 Oct 18:50): two independent landscapes A
        # (the table above, "cold") and B ("hot"); the selected score is z_env = ((1-w) z_A + w z_B) / norm(w), with
        # norm(w) = sqrt((1-w)^2 + w^2 + 2 w (1-w) rho) so that its SD over random sequences stays 1 (rho = correlation
        # of z_A and z_B over random domains, Monte Carlo). w(t) follows a schedule: constant (0), season (sinusoid,
        # period 1 simulated year), trend (triangle 0 -> 1 -> 0, period 200 years), catastrophe (0, or 1 for
        # shift_years after a catastrophe at the site). The control half uses per-step shuffles of both tables.
        self.env = dict(cfg.get("env_landscape") or {}) if cfg.get("env_landscape") else None
        self._t_h, self._cat_last, self.w_now = 0.0, None, 0.0
        if self.env:
            if self.landscape != "domain":
                raise ValueError("env_landscape needs the domain landscape")
            import envsched
            if self.env.get("variable", "T") not in envsched.VARIABLES:
                raise ValueError("env_landscape variable must be one of " + ", ".join(envsched.VARIABLES))
            gB = torch.Generator().manual_seed(self.seed + 7777)
            self.tableB = (torch.randn(nw, 1 << self.k, generator=gB) if self.table.dim() == 2 else torch.randn(1 << self.k, generator=gB)).to(dev)
            self.sdB = self._sd_for(self.tableB)
            self._step_tabB = self.tableB
            self._shuffle_genB = torch.Generator().manual_seed(self.seed + 7778)
            self.rho = self._rho_ab()
        self._prom = None                  # Gokyona protected domain: chains without a start mark are not scored

    def _sd_for(self, table):
        """SD of the domain score of random sequences, per number of windows n (Monte Carlo, once)."""
        nw = self.D - self.k + 1
        gm = torch.Generator().manual_seed(self.seed + 99)
        bits = torch.randint(0, 2, (20000, self.D), generator=gm)
        idx = torch.zeros(20000, nw, dtype=torch.long)
        for i in range(self.k):
            idx |= bits[:, i:i + nw] << i
        tab = table.cpu()
        v = tab[torch.arange(nw)[None, :], idx] if tab.dim() == 2 else tab[idx]
        sd = torch.ones(nw + 1)
        for n in range(1, nw + 1):
            sd[n] = v[:, :n].mean(1).std().clamp(min=1e-6)
        return sd.to(self.dev)

    def _rho_ab(self):
        """Correlation of the (pair) scores z_A and z_B over random 64-letter chains."""
        g = torch.Generator().manual_seed(self.seed + 4242)
        n = 20000
        import words as _wd
        L = torch.full((n,), 64, dtype=torch.long)
        seq = torch.zeros(n, _wd.W, dtype=torch.long)
        seq[:, 0] = torch.randint(-(1 << 62), 1 << 62, (n,), generator=g) * 2 + torch.randint(0, 2, (n,), generator=g)
        L, seq = L.to(self.dev), seq.to(self.dev) & _wd.mask(L.to(self.dev))
        za, zb = self.zscore_ab(L, seq)
        return float(torch.corrcoef(torch.stack([za, zb]))[0, 1])

    def w_of(self, t_h, rows=None):
        """Schedule value w(t) (scalar), or per row [n] when rows are given (catastrophe shift variant)."""
        import envsched
        e = self.env
        sch = e.get("schedule", "constant")
        w = envsched.w_total(e, t_h, e.get("impact"), getattr(self, "_impact_last", None))
        if e.get("gradient") and rows is not None:
            # spatial gradient (coordinator 29 Oct 20:04): w_s = s / (n_half - 1) along a line of sites in each half
            # (mirrored in the control half); moving bands: w_s(t) = clip(w_s + band_amp * tri(t) + a_s sin(2 pi t / 1 y),
            # 0, 1), tri = the trend triangle 0 -> 1 -> 0 with period band_period_years (band_amp 0 = static); a_s =
            # season_amp (1 - w_s) = latitude-dependent season (season_amp 0 = none; coordinator 29 Oct 20:13)
            half = max(1, self.N // 2)
            s_loc = (rows % self.N) % half
            base = s_loc.float() / max(1, half - 1)
            w_lat = base
            amp = float(e.get("band_amp", 0) or 0)
            if amp:
                import envsched
                base = base + amp * envsched.w_at({"schedule": "trend", "period_years": float(e.get("band_period_years", 100.0))}, t_h)
            sa = float(e.get("season_amp", 0) or 0)
            if sa:                                 # latitude-dependent season: a_s = season_amp (1 - w_s), largest at the
                import math                        # cold pole, ~0 at the hot end; a_s sin(2 pi t / season_period)
                base = base + sa * (1 - w_lat) * math.sin(2 * math.pi * t_h / (float(e.get("season_period_years", 1.0)) * 8766.0))
            return base.clamp(0.0, 1.0)
        if rows is None or sch != "catastrophe" or not float(e.get("shift_years", 0) or 0) or self._cat_last is None:
            return w
        recent = (t_h - self._cat_last.to(rows.device)[rows]) < float(e["shift_years"]) * 8766.0
        return torch.where(recent, torch.ones(rows.numel(), device=rows.device), torch.full((rows.numel(),), w, device=rows.device))

    def site_w(self, t_h):
        """Per-site w over the N sites of a run (gradient), or the scalar schedule value repeated."""
        rows = torch.arange(self.N, device=self.dev)
        w = self.w_of(t_h, rows)
        return w if torch.is_tensor(w) else torch.full((self.N,), float(w), device=self.dev)

    def zscore_ab(self, L, seq, ctrl_rows=None):
        """Pair scores under landscape A and B (real tables; control rows use this step's shuffles)."""
        za = self._pair(L, seq, ctrl_rows, self.table, self._step_tab, self.sd)
        zb = self._pair(L, seq, ctrl_rows, self.tableB, self._step_tabB, self.sdB)
        return za, zb

    def _pair(self, L, seq, ctrl_rows, tab, step_tab, sd):
        z = self.score(L, seq, ctrl_rows, table=tab, step_table=step_tab, sd=sd)
        if self.pair:
            comp = torch.bitwise_not(seq) & wd.mask(L)
            z = 0.5 * (z + self.score(L, comp, ctrl_rows, table=tab, step_table=step_tab, sd=sd))
        return z

    def fidelity_on(self):
        return self.enabled and self.fid and self.s_f != 0

    # ---- scores ----
    def score(self, L, seq, ctrl_rows=None, table=None, step_table=None, sd=None):
        """Scores of chains (L [n], seq [n, W]); rows in ctrl use the step's shuffled table if given."""
        n = L.numel()
        out = torch.zeros(n, device=L.device)
        if n == 0:
            return out
        k = self.k
        dom = self.landscape == "domain"
        Le = L.clamp(max=self.D) if dom else L
        lmax = int(Le.max())
        if lmax < k:
            return out
        bits = wd.unpack(seq, lmax)                                   # [n, lmax]
        nw = lmax - k + 1
        idx = torch.zeros(n, nw, dtype=torch.long, device=L.device)
        for i in range(k):
            idx |= bits[:, i:i + nw] << i
        nvalid = (Le - k + 1).clamp(min=0)
        valid = torch.arange(nw, device=L.device)[None, :] < nvalid[:, None]
        tab = self.table if table is None else table
        pos = torch.arange(nw, device=L.device)[None, :].expand(n, nw)
        look = (lambda T: T[pos, idx]) if tab.dim() == 2 else (lambda T: T[idx])
        v = look(tab)
        if ctrl_rows is not None and bool(ctrl_rows.any()):
            v = torch.where(ctrl_rows[:, None], look(self._step_tab if step_table is None else step_table), v)
        cnt = valid.sum(1).clamp(min=1)
        sc = torch.where(valid, v, torch.zeros_like(v)).sum(1) / cnt
        if dom:
            sdt = self.sd if sd is None else sd
            sc = sc / sdt[nvalid.clamp(max=sdt.numel() - 1)]
            return sc * (L >= 12)
        return sc * (L >= k)

    def zscore(self, L, seq, ctrl_rows=None, rows=None):
        """Trait value of chains: own score, or the mean with the complement's (pair score); with the environment-
        dependent landscape the mix z_env(t) (rows: per-row w for the catastrophe shift variant)."""
        if self.env:
            za, zb = self.zscore_ab(L, seq, ctrl_rows)
            w = self.w_of(self._t_h, rows)
            norm = ((1 - w) ** 2 + w ** 2 + 2 * w * (1 - w) * self.rho) ** 0.5
            return ((1 - w) * za + w * zb) / norm
        z = self.score(L, seq, ctrl_rows)
        if self.pair:
            comp = torch.bitwise_not(seq) & wd.mask(L)
            z = 0.5 * (z + self.score(L, comp, ctrl_rows))
        return z

    def all_scores(self, L, seq, chem):
        """[R, S] trait values of all nucleotide-like chains (>= 12 letters for the domain landscape), once per step."""
        R, S = L.shape
        z = torch.zeros(R, S, device=L.device)
        m = (chem == 0) & (L >= (12 if self.landscape == "domain" else self.k))
        if self._prom is not None:
            m = m & self._prom
        rr, cc = m.nonzero(as_tuple=True)
        if rr.numel():
            z[rr, cc] = self.zscore(L[rr, cc], seq[rr, cc], self.is_ctrl(rr), rows=rr if self.env else None)
        return z

    def stability(self, L, seq, chem):
        """First trait call of a step (bond breaking): new control table, all scores; multiplier of k_break."""
        self.new_step()
        if not self.enabled:
            self._all = None
            return 1.0
        self._all = self.all_scores(L, seq, chem)
        if self.mode == "stability" and self.s != 0:
            return torch.exp(-self.s * self._all)
        return 1.0

    def vesicle_mults(self, L, seq, chem, lab, R):
        """(membrane uptake, division probability, washout rate) multipliers per vesicle [R, NV] or None."""
        if not self.enabled or self._all is None or self.mode not in ("membrane", "division", "retention"):
            return None, None, None
        zm, has = self.vesicle_z(L, chem, lab, R)
        if self.mode == "membrane":
            return self._fill(torch.exp(self.s * zm), has), None, None
        if self.mode == "division":
            return None, self._fill(torch.exp(self.s * zm), has), None
        return None, None, self._fill(torch.exp(-self.s * zm), has)

    def _fill(self, mult, has):
        """Multiplier for vesicles without scored chains: 1 (the study's rule), or with the protected domain on (F1,
        coordinator 29 Oct 13:07) the mean multiplier of the scored vesicles of the same half in this step, so that
        having or lacking scored chains carries no expected advantage (exp is convex: E[exp(-s z)] > 1 for random z)."""
        if self._prom is None:
            return mult * has + (~has)
        out = torch.where(has, mult, torch.zeros_like(mult))
        ctl = self.is_ctrl(torch.arange(mult.shape[0], device=mult.device))[:, None].expand_as(has)
        for hm in (ctl, ~ctl):
            sel = has & hm
            fill = mult[sel].mean() if bool(sel.any()) else torch.ones((), device=mult.device)
            out = torch.where(~has & hm, fill, out)
        return out

    def vesicle_z(self, L, chem, lab, R):
        """Mean trait value of the scored chains enclosed in each vesicle [R, NV], and whether it has any."""
        NV = _nv()
        inside = (lab > 0) & (chem == 0) & (L >= 12)
        if self._prom is not None:
            inside = inside & self._prom
        rr, cc = inside.nonzero(as_tuple=True)
        tot = torch.zeros(R * NV, device=L.device)
        num = torch.zeros(R * NV, device=L.device)
        if rr.numel():
            comp = rr * NV + lab[rr, cc].long() - 1
            tot.index_add_(0, comp, self._all[rr, cc])
            num.index_add_(0, comp, torch.ones(rr.numel(), device=L.device))
        return (tot / num.clamp(min=1)).view(R, NV), (num.view(R, NV) > 0)

    def membrane(self, L, seq, chem, lab, R):
        if not self.enabled or self.mode != "membrane" or self._all is None:
            return None
        NV = _nv()
        inside = (lab > 0) & (chem == 0) & (L >= 12)
        rr, cc = inside.nonzero(as_tuple=True)
        tot = torch.zeros(R * NV, device=L.device)
        num = torch.zeros(R * NV, device=L.device)
        if rr.numel():
            comp = rr * NV + lab[rr, cc].long() - 1
            tot.index_add_(0, comp, self._all[rr, cc])
            num.index_add_(0, comp, torch.ones(rr.numel(), device=L.device))
        zm = (tot / num.clamp(min=1)).view(R, NV)
        return torch.exp(self.s * zm) * (num.view(R, NV) > 0) + (num.view(R, NV) == 0)

    def is_ctrl(self, rows):
        return (rows % self.N) >= self.N // 2

    def new_step(self):
        nk = self.table.shape[-1]
        if self.table.dim() == 2:
            perms = torch.stack([torch.randperm(nk, generator=self._shuffle_gen) for _ in range(self.table.shape[0])]).to(self.dev)
            self._step_tab = torch.gather(self.table, 1, perms)
        else:
            self._step_tab = self.table[torch.randperm(nk, generator=self._shuffle_gen).to(self.dev)]
        if self.env:
            if self.tableB.dim() == 2:
                pB = torch.stack([torch.randperm(nk, generator=self._shuffle_genB) for _ in range(self.tableB.shape[0])]).to(self.dev)
                self._step_tabB = torch.gather(self.tableB, 1, pB)
            else:
                self._step_tabB = self.tableB[torch.randperm(nk, generator=self._shuffle_genB).to(self.dev)]
            self.w_now = self.w_of(self._t_h)
        self._cache = None

    # ---- hooks called by gk_chem ----
    def rate(self, r_m, trr, tcc, L, seq, chem, lab):
        if self._all is None or trr.numel() == 0:
            self._cache = (trr, tcc, torch.zeros(trr.numel(), device=r_m.device))
            return r_m
        sc_t = self._all[trr, tcc]
        self._cache = (trr, tcc, sc_t)
        if self.mode not in ("self", "shared") or self.s == 0:
            return r_m
        if self.mode == "self":
            return r_m * torch.exp(self.s * sc_t)
        # shared replicase: mean score of chains >= L_rep (chemistry 0) per compartment (row, label)
        R, S = L.shape
        rep = (L >= self.L_rep) & (chem == 0)
        rr, cc = rep.nonzero(as_tuple=True)
        if rr.numel() == 0:
            return r_m
        sc_r = self._all[rr, cc]
        nc = _nv() + 1 if lab is not None else 1
        comp = rr * nc + (lab[rr, cc].long() if lab is not None else 0)
        tot = torch.zeros(R * nc, device=r_m.device).index_add_(0, comp, sc_r)
        num = torch.zeros(R * nc, device=r_m.device).index_add_(0, comp, torch.ones_like(sc_r))
        mean = tot / num.clamp(min=1)
        tcomp = trr * nc + (lab[trr, tcc].long() if lab is not None else 0)
        boost = 1 + self.s * mean[tcomp].clamp(min=0) * (num[tcomp] > 0)
        return r_m * boost

    def eps(self, eps_t, trr, tcc):
        if self._cache is None:
            return eps_t
        return eps_t * torch.exp(-self.s_f * self._cache[2])
