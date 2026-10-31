"""P3-1b (Review 4, R4.3): copy of lot-s4/chem4c5.py (md5 c747dc68, not modified) on chem4c4b (the P3-1b switches
MULTI_START, SLOT_ROUND_FIX, POISSON_INV, default off). Set the switches on the chem4c4b module (chem4c5b._c), not here:
`from chem4c4b import *` copies their values at import.
--- chem4c5 description follows ---
Stage 4c-5: vesicle lineage record and the L3 / L4 detectors (spec docs/stage4c5_spec.md, commit
05017cd, md5 f11bba4b; L3 strict revised 28 Oct 23:00 (paired random-copy condition), 23:30 (only
chains made after the ancestor's checkpoint count) and 29 Oct 00:10 (local ancestor 3 generations up,
>= 90% sequence identity)). No new chemistry or physics: chem4c4 v2 (md5 b79360f6) runs unchanged; this
module wraps its Vesicles class with recording hooks that only read state and draw no random
numbers (validation test 1: bit-identical to chem4c4).

Record (spec §1), kept on the GPU per site for the ancestry and on the host per checkpoint:
  * id table per site (MAXID ids): parent id, generation, birth checkpoint, end checkpoint and end
    cause (1 washout, 2 breakdown outside the window, 3 dissolved below the CVC); formation births
    have parent -1; a second daughter without a free slot is never given an id (division event with
    d2 = -1, the "missing daughter slot" loss);
  * per vesicle slot per checkpoint (host): id, parent, generation, membrane a, volume share phi,
    membrane gained since the previous checkpoint (cumulative uptake: growth from non-vesicle
    membrane, exchange and dissolution losses; formation and the division halving excluded),
    divisions in the slot since the previous checkpoint, enclosed letters, enclosed letters of
    chains of length >= L_H, osmotic load Pi (as the dynamics compute it), the two dominant enclosed
    lineage keys (whole-sequence species keys of chemistry-0 chains >= L_H; a sequence and its
    complement count as one lineage, keyed by the smaller key) with their enclosed letters; a
    host table key -> (length, sequence) for the stored keys;
  * division events per site at checkpoint resolution: parent id, daughter ids, interval;
  * per chain slot the step at which the chain in it was made (int32, -1 = present at the start): set
    when the slot's (length, sequence, chemistry) changes (a copy or initiation into a free slot, a
    primed copy growing in its primer's slot, a fragment, elongation); per run the step of each
    checkpoint; per checkpoint, vesicle and dominant key the letters of that key in chains made
    after checkpoint m, for every m (top_ln [..., 2, N_CHECK]); at each run's last checkpoint a host
    snapshot of its enclosed chemistry-0 chains >= L_H (site, vesicle slot, length, sequence, birth
    step) and the per-site checkpoint steps (final_chains; host-side analysis only).
Detectors (spec §2), computed from the record at the end (Lineage.detect):
  * joint window: longest run of checkpoints ending at the last one with E3 strict at the
    checkpoint (>= half the sites hold a vesicle enclosing a chain >= L_H) and the L2 lineage
    present (largest site share of the L2 detector >= L2_SHARE);
  * evaluable: L2_long and E3_strict, and the window holds >= 30 (vesicle, interval) pairs and a
    generation span >= 3;
  * L3 loose: Spearman rho(g, Pi) >= 0.1 with permutation p < 0.01 (999 permutations of Pi among the
    vesicles of a site and checkpoint); g = gain / mean membrane over the interval, Pi at its start;
  * L3 strict (local, spec 05017cd): L3 loose and >= half the sites hold a vesicle, alive at the last
    checkpoint, that encloses >= 0.5 x its ancestor's dominant-key letters in chains made after the
    ancestor's last recorded checkpoint whose sequence matches the key with >= 90% identity (either
    strand, compared over the shorter length from position 0, the copy start; chemistry 0, >= L_H).
    The ancestor is 3 generations up (great-grandparent); its key and letters are those of the last
    checkpoint at which it was alive. If the great-grandparent was never alive at a checkpoint (born
    and divided within one interval), the nearest ancestor further up that was is used (depth
    reported). Chains that already existed do not count, however often re-enclosed. And (paired_inherit) the share of sites passing exceeds that of
    the paired random-copy run (same environment, constants, seed, window) by >= L3_RC_MARGIN; detect()
    gives the unpaired verdict, the caller applies paired_inherit (L3_strict and L4 can only lose);
  * L4: L3 strict, heritability (parent -> daughter regression of z, slope > 0, p < 0.01) and
    selection: Spearman rho(w, z) over >= 10 ancestors above the 99th percentile of 999 permutations
    of z among the ancestors; z = enclosed lineage letters per volume share; w = number of
    descendants alive at the last checkpoint; L4 sequence version reported (top rising keys).
    (z uses the letters of the two dominant enclosed lineages, chemistry 0, >= L_H.)
Records are written to disk only on request (lineage_dir), for evaluable runs.
Speed change (16), 29 Oct: the per-run detectors run in DET_WORKERS forked processes; the L3-loose and
rho(w, z) nulls use the same permutation draws as the loop version and compute the null correlations from
permuted ranks; paired reruns build the record only. Checked bit-identical (all outputs and details) against
the loop version on 22 runs (4 joint points x 2, test-4 and test-3 fixtures); 2.00 -> 0.32 s per run.
"""
import math, os, pickle
import numpy as np
import torch
from scipy.stats import spearmanr, linregress, rankdata
import chem4c4b as _c
from chem4c4b import *  # noqa: F401,F403  (the 4c-4 model, unchanged)
from chem4c4b import _e3_empty  # noqa: F401

MAXID = 2048              # ids per site in the ancestry table (overflow flagged)
L3_RHO, L3_P, L3_PAIRS, L3_GENS, NPERM = 0.1, 0.01, 30, 3, 999
L3_RC_MARGIN = 0.25       # inheritance share above the paired random-copy run (spec §2, revised bd0f418)
PERTURB_SETS = None       # Review 9: [(name, {constant: value})]: paired_inherit re-evaluates local_inherit on the records in memory
PERTURB_OUT = None        # set by paired_inherit: {name: (share main, share main', share c)} per run (nan outside the mask)
L3_UP, L3_IDENT, L3_FRAC = 3, 0.9, 0.5   # local inheritance: generations up, sequence identity, letter fraction (05017cd)
_Base = _c.Vesicles


class Vesicles(_Base):
    """chem4c4.Vesicles plus the lineage record; reads state only, draws no random numbers."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        R, dev = self.R, self.dev
        self.t_par = torch.full((R, MAXID), -2, dtype=torch.int32, device=dev)
        self.t_gen = torch.zeros(R, MAXID, dtype=torch.int16, device=dev)
        self.t_birth = torch.full((R, MAXID), -1, dtype=torch.int16, device=dev)
        self.t_end = torch.full((R, MAXID), -1, dtype=torch.int16, device=dev)
        self.t_cause = torch.zeros(R, MAXID, dtype=torch.int8, device=dev)
        self.overflow = torch.zeros(R, dtype=torch.bool, device=dev)
        self.cg = torch.zeros(R, NV, device=dev)                 # membrane gained since last checkpoint
        self.ndiv_slot = torch.zeros(R, NV, dtype=torch.int16, device=dev)
        self.k_now = torch.zeros(R, dtype=torch.int16, device=dev)   # current checkpoint interval
        self.cont = torch.ones(R, NV, dtype=torch.bool, device=dev)  # slot lineage unbroken since last checkpoint
        self.chk = []                                             # host records per checkpoint
        self.div_events = []                                      # host: (rows, parent, d1, d2, k)
        self.key_table = {}
        self.nstep = 0                                            # hook step counter (chain birth record)
        self.ck_step = torch.zeros(R, N_CHECK, dtype=torch.int32, device=dev)
        self.final_chains = {}                                    # run -> host snapshot at its last checkpoint

    def _births(self, L, seq, chem):
        """Chain birth record: slots whose (L, seq, chem) changed since the last call get the current step."""
        ch = (L != self.pL) | (seq != self.pseq) | (chem != self.pchem)
        self.cb = torch.where(ch & (L > 0), torch.full_like(self.cb, self.nstep), self.cb)
        self.cb = torch.where(L == 0, torch.full_like(self.cb, -1), self.cb)
        self.pL, self.pseq, self.pchem = L.clone(), seq.clone(), chem.clone()

    # ---- id table helpers ----
    def _register(self, mask, parent, gen):
        rr, cc = mask.nonzero(as_tuple=True)
        if rr.numel() == 0:
            return
        ids = self.vid[rr, cc]
        ok = ids < MAXID
        self.overflow.index_fill_(0, rr[~ok], True)
        rr, ids = rr[ok], ids[ok].long()
        self.t_par[rr, ids] = parent[mask][ok].int()
        self.t_gen[rr, ids] = gen[mask][ok].to(torch.int16)
        self.t_birth[rr, ids] = self.k_now[rr]

    def _end(self, mask, cause):
        rr, cc = mask.nonzero(as_tuple=True)
        if rr.numel() == 0:
            return
        ids = self.vid[rr, cc].long()
        ok = (ids >= 0) & (ids < MAXID)
        self.t_end[rr[ok], ids[ok]] = self.k_now[rr[ok]]
        self.t_cause[rr[ok], ids[ok]] = cause

    # ---- hooks ----
    def init_chains(self, L, seq, chem, pair, dlen, lab, special):
        super().init_chains(L, seq, chem, pair, dlen, lab, special)
        ex = (special.get("init_vesicles") or {}).get("extra")
        if ex:                               # extra chains placed in given vesicle slots (validation fixtures)
            j0 = int((L[0] > 0).sum())
            for slot in ex["slots"]:
                for c in range(int(ex["n"])):
                    L[:, j0] = int(ex["L"]); seq[:, j0] = int(ex["seq"]); chem[:, j0] = int(ex.get("chem", 0)); lab[:, j0] = int(slot) + 1
                    j0 += 1
        if self.alive.any():
            self._register(self.alive, torch.full_like(self.vid, -1), torch.zeros_like(self.vgen))
        self.cb = torch.full(L.shape, -1, dtype=torch.int32, device=self.dev)   # chains present at the start
        self.pL, self.pseq, self.pchem = L.clone(), seq.clone(), chem.clone()

    def _form(self, t, nnew, L, lab, pair):
        before = self.alive.clone()
        super()._form(t, nnew, L, lab, pair)
        new = self.alive & ~before
        self._register(new, torch.full_like(self.vid, -1), torch.zeros_like(self.vgen))
        self._pre = torch.where(new, torch.full_like(self._pre, A0_VES), self._pre)
        self.cont &= ~new

    def _divide(self, t, div, L, lab, pair):
        old_id, old_alive, a_m = self.vid.clone(), self.alive.clone(), self.va.clone()
        super()._divide(t, div, L, lab, pair)
        born = self.alive & ~old_alive
        # parents end by division (cause 4), daughters registered with parent = old id
        rr, cc = div.nonzero(as_tuple=True)
        if rr.numel():
            ids = old_id[rr, cc].long()
            ok = (ids >= 0) & (ids < MAXID)
            self.t_end[rr[ok], ids[ok]] = self.k_now[rr[ok]]
            self.t_cause[rr[ok], ids[ok]] = 4
        self._register(div, self.vpar, self.vgen)
        self._register(born, self.vpar, self.vgen)
        self.cont &= ~born                    # a second daughter starts a new slot lineage
        self.ndiv_slot = self.ndiv_slot + div.to(torch.int16)
        # gain bookkeeping: the halving is not uptake
        self._pre = torch.where(div, a_m / 2, self._pre)
        self._pre = torch.where(born, self.va, self._pre)
        # division events (host, aggregated to checkpoint interval)
        if rr.numel():
            d1 = self.vid[rr, cc]
            par = old_id[rr, cc]
            d2 = torch.full_like(d1, -1)
            br, bc = born.nonzero(as_tuple=True)
            if br.numel():
                key_b = br.long() * (1 << 32) + self.vpar[br, bc].long()
                key_d = rr.long() * (1 << 32) + par.long()
                srt, idx = torch.sort(key_b)
                pos = torch.searchsorted(srt, key_d).clamp(max=srt.numel() - 1)
                hit = srt[pos] == key_d
                d2 = torch.where(hit, self.vid[br[idx[pos]], bc[idx[pos]]], d2)
            self.div_events.append(torch.stack([rr.long(), par.long(), d1.long(), d2.long(), self.k_now[rr].long()], 1).cpu())

    def step(self, t, dt, ph, L, seq, chem, pair, dlen, lab, unpair, Ms):
        self.nstep += 1
        self._births(L, seq, chem)
        alive0 = self.alive.clone()
        self._pre = self.va.clone()
        vid0 = self.vid.clone()
        win = (self.win_run & (t < self.close_t)).to(self.dev)[self.rr]
        super().step(t, dt, ph, L, seq, chem, pair, dlen, lab, unpair, Ms)
        gone = alive0 & ~self.alive & (self.vid == vid0)          # same id, no longer alive (not divided)
        if bool(gone.any()):
            small = gone & (self._pre < 0.5 * A0_VES)
            out = gone & ~win[:, None] & ~small
            self._end_ids(vid0, small, 3)
            self._end_ids(vid0, out, 2)
            self._end_ids(vid0, gone & ~small & ~out, 1)
        self.cg = torch.where(self.alive, self.cg + (self.va - self._pre), torch.zeros_like(self.cg))
        self.cont &= self.alive | ~alive0          # emptied slot -> lineage broken

    def _end_ids(self, vid0, mask, cause):
        rr, cc = mask.nonzero(as_tuple=True)
        if rr.numel() == 0:
            return
        ids = vid0[rr, cc].long()
        ok = (ids >= 0) & (ids < MAXID)
        self.t_end[rr[ok], ids[ok]] = self.k_now[rr[ok]]
        self.t_cause[rr[ok], ids[ok]] = cause

    def checkpoint(self, runs2, kidx, rows2, L, lab, chem):
        super().checkpoint(runs2, kidx, rows2, L, lab, chem)
        self._births(L, self.seq, chem)
        kr = self.k_now[rows2].long().clamp(max=N_CHECK - 1)
        self.ck_step[rows2, kr] = self.nstep
        Lr, lr, cr = L[rows2], lab[rows2], chem[rows2]
        n = rows2.numel()
        inside = lr > 0
        j = (lr.long() - 1).clamp(min=0)
        Lf = Lr.float()
        rid = torch.arange(n, device=self.dev)[:, None].expand_as(Lr)
        enc = torch.zeros(n * NV, device=self.dev).index_add_(0, (rid * NV + j)[inside], Lf[inside]).view(n, NV)
        lg = inside & (Lr >= L_H)
        encl = torch.zeros(n * NV, device=self.dev).index_add_(0, (rid * NV + j)[lg], Lf[lg]).view(n, NV)
        outl = torch.where((lr == 0) & (Lr > 0), Lf, torch.zeros_like(Lf)).sum(1)
        va, al = self.va[rows2], self.alive[rows2]
        phi = (self.phi0[rows2][:, None] * (va / A0_VES) ** 1.5).clamp(max=1.0) * al
        outvol = (1 - phi.sum(1)).clamp(min=0.05)
        Pi = enc / phi.clamp(min=1e-30) / (outl / outvol + 1.0)[:, None]
        # dominant lineage keys (chemistry 0, L >= L_H, a sequence and its complement as one lineage)
        seqr = self.seq[rows2]
        m0 = lg & (cr == 0)
        k1 = _c.cat4c2.species_key(cr, Lr, seqr)
        k2 = _c._comp_key(Lr, seqr)
        key = torch.where(m0, torch.minimum(k1, k2), torch.zeros_like(k1))
        top_k = torch.zeros(n, NV, 2, dtype=torch.int64, device=self.dev)
        top_l = torch.zeros(n, NV, 2, device=self.dev)
        top_ln = torch.zeros(n, NV, 2, N_CHECK, device=self.dev)
        rr, cc = m0.nonzero(as_tuple=True)
        if rr.numel():
            g = torch.stack([rr, j[rr, cc], key[rr, cc]], 1)
            ug, inv = torch.unique(g, dim=0, return_inverse=True)
            let = torch.zeros(ug.shape[0], device=self.dev).index_add_(0, inv, Lf[rr, cc])
            # letters in chains made after checkpoint m (m < current checkpoint; later m: none)
            cbr = self.cb[rows2][rr, cc].long()
            cks = self.ck_step[rows2][rr].long()                                  # [nc, N_CHECK]
            past = torch.arange(N_CHECK, device=self.dev)[None, :] < kr[rr][:, None]
            newm = ((cbr[:, None] > cks) & past).float() * Lf[rr, cc][:, None]
            letn = torch.zeros(ug.shape[0], N_CHECK, device=self.dev).index_add_(0, inv, newm)
            o = torch.argsort(-let, stable=True)
            ug, let, letn = ug[o], let[o], letn[o]
            o2 = torch.argsort(ug[:, 0] * NV + ug[:, 1], stable=True)
            ug, let, letn = ug[o2], let[o2], letn[o2]
            grp = ug[:, 0] * NV + ug[:, 1]
            first = torch.searchsorted(grp, grp)
            rank = torch.arange(grp.numel(), device=self.dev) - first
            sel = rank < 2
            top_k[ug[sel, 0], ug[sel, 1], rank[sel]] = ug[sel, 2]
            top_l[ug[sel, 0], ug[sel, 1], rank[sel]] = let[sel]
            top_ln[ug[sel, 0], ug[sel, 1], rank[sel]] = letn[sel]
            # host table key -> (L, seq) for the stored keys
            ks = set(ug[sel, 2].tolist())
            new = [x for x in ks if x not in self.key_table]
            if new:
                kk = torch.tensor(new, device=self.dev)
                hitm = m0 & torch.isin(key, kk)
                hr, hc = hitm.nonzero(as_tuple=True)
                for x, l_, s_ in zip(key[hr, hc].tolist(), Lr[hr, hc].tolist(), seqr[hr, hc].tolist()):
                    self.key_table.setdefault(x, (l_, s_))
        rec = {"k": kidx.cpu(), "runs": runs2.cpu(), "vid": self.vid[rows2].cpu(), "par": self.vpar[rows2].cpu(),
               "gen": self.vgen[rows2].cpu(), "alive": al.cpu(), "a": va.cpu(), "phi": phi.cpu(),
               "gain": self.cg[rows2].cpu(), "ndiv": self.ndiv_slot[rows2].cpu(), "enc": enc.cpu(), "encl": encl.cpu(),
               "Pi": Pi.cpu(), "top_k": top_k.cpu(), "top_l": top_l.cpu(), "top_ln": top_ln.cpu(), "site_letters": Lf.sum(1).cpu(), "out_letters": outl.cpu(),
               "cont": (self.cont[rows2] & al).cpu()}
        self.chk.append(rec)
        lastm = kr == N_CHECK - 1
        if bool(lastm.any()):
            # host snapshot of the enclosed chemistry-0 chains >= L_H of runs at their last checkpoint
            m = m0 & lastm[:, None]
            rr_, cc_ = m.nonzero(as_tuple=True)
            run_r = self.rr[rows2][rr_].cpu().numpy()
            site = (rows2[rr_] % self.N).cpu().numpy().astype(np.int16)
            vs = j[rr_, cc_].cpu().numpy().astype(np.int8)
            Lc = Lr[rr_, cc_].cpu().numpy().astype(np.int16)
            sq = seqr[rr_, cc_].cpu().numpy()
            bs = self.cb[rows2][rr_, cc_].cpu().numpy()
            cks = self.ck_step[rows2].cpu().numpy()
            runs_l = self.rr[rows2][lastm].cpu().numpy()
            for b in np.unique(runs_l):
                sel = run_r == b
                rows_b = (self.rr[rows2] == int(b)).cpu().numpy()
                self.final_chains[int(b)] = {"site": site[sel], "slot": vs[sel], "L": Lc[sel], "seq": sq[sel], "born": bs[sel],
                                             "ck_step": cks[rows_b].reshape(-1, N_CHECK)}
        self.cont[rows2] = True
        self.cg[rows2] = 0
        self.ndiv_slot[rows2] = 0
        self.k_now[rows2] = self.k_now[rows2] + 1

    def outputs(self, t, t_end, l2_t, next_chk2, L, chem, lab, out):
        res = super().outputs(t, t_end, l2_t, next_chk2, L, chem, lab, out)
        if Lineage._record_only:              # paired reruns: only the record is used (local_inherit)
            res["lineage"] = Lineage(self, out, res)
        else:
            res.update(Lineage(self, out, res).detect(Lineage._override))
        return res


_POP8 = np.array([bin(i).count("1") for i in range(256)], np.int64)    # popcount per byte


def _spearman(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(spearmanr(x, y).correlation)


def _spearman_null(x, y, rng):
    """[_spearman(x, rng.permutation(y)) for NPERM draws], computed as scipy does (ranks of x and of the
    permuted y, np.corrcoef of the n x 2 rank array) without re-ranking: the ranks of a permutation are the
    permuted ranks. Same draws (rng.permutation(len(y)) consumes the draws of rng.permutation(y))."""
    n = len(x)
    if n < 3 or np.std(x) == 0 or np.std(y) == 0:
        for _ in range(NPERM):
            rng.permutation(n)
        return [float("nan")] * NPERM
    rx, ry = rankdata(x), rankdata(y)
    return [float(np.corrcoef(np.column_stack((rx, ry[rng.permutation(n)])), rowvar=False)[1, 0]) for _ in range(NPERM)]


DET_WORKERS = int(os.environ.get("LOT_DET_WORKERS", "4"))   # processes for the per-run detectors (speed change (16))
DET_TIMEOUT = 1800.0                                         # s; on timeout the runs are evaluated serially
_DET = None


def _det_worker(a):
    return _DET._detect_one(*a)


class Lineage:
    """Assembles the host record per run and evaluates the L3 / L4 detectors (spec §2)."""
    _override = None      # fixture tests: evaluate from this checkpoint on, without the L2/E3 gate
    _record_only = False  # paired reruns: build the record, skip the detectors

    def __init__(self, V, out, res):
        self.V, self.out, self.res = V, out, res
        B, N = V.B, V.N
        self.B, self.N = B, N
        K = N_CHECK
        shp = (B, N, K, NV)
        f = lambda dt: np.zeros(shp, dt)
        self.r = {k: f(np.float32) for k in ("a", "phi", "gain", "enc", "encl", "Pi")}
        for k in ("vid", "par", "gen", "ndiv"):
            self.r[k] = np.full(shp, -1 if k in ("vid", "par") else 0, np.int64)
        self.r["alive"] = np.zeros(shp, bool)
        self.r["cont"] = np.zeros(shp, bool)
        self.r["top_k"] = np.zeros(shp + (2,), np.int64)
        self.r["top_l"] = np.zeros(shp + (2,), np.float32)
        self.r["top_ln"] = np.zeros(shp + (2, K), np.float32)    # letters of the key in chains made after checkpoint m
        self.r["site_letters"] = np.zeros((B, N, K), np.float32)
        self.r["out_letters"] = np.zeros((B, N, K), np.float32)
        self.rec = np.zeros((B, K), bool)
        for c in V.chk:
            runs, k = c["runs"].numpy(), c["k"].numpy()
            n = len(runs)
            for name in self.r:
                v = c[name].numpy()
                v = v.reshape((n, N) + v.shape[1:])
                self.r[name][runs, :, k] = v
            self.rec[runs, k] = True
        self.par_t = V.t_par.cpu().numpy().reshape(B, N, MAXID)
        self.gen_t = V.t_gen.cpu().numpy().reshape(B, N, MAXID)
        self.end_t = V.t_end.cpu().numpy().reshape(B, N, MAXID)
        self.cause_t = V.t_cause.cpu().numpy().reshape(B, N, MAXID)
        self.ovf = V.overflow.cpu().numpy().reshape(B, N).any(1)
        ev = torch.cat(V.div_events).numpy() if V.div_events else np.zeros((0, 5), np.int64)
        self.div = ev   # rows (global row = run * N + site), parent, d1, d2, k

    # ---- pieces ----
    def joint_window(self, b):
        sh = self.out["L2_share_matrix"][b].numpy() if torch.is_tensor(self.out["L2_share_matrix"]) else self.out["L2_share_matrix"][b]
        l2k = sh.max(0) >= L2_SHARE                                  # [K]
        e3k = self.V.e3_enc[b].float().mean(0).numpy() >= 0.5         # [K]
        ok = l2k & e3k & self.rec[b]
        K = len(ok)
        k0 = K
        while k0 > 0 and ok[k0 - 1]:
            k0 -= 1
        return k0 if k0 < K else None                                   # first checkpoint of the window

    def pairs(self, b, k0, K=N_CHECK):
        """(slot lineage, interval) pairs inside [k0, K-1]: g = gain / mean a, Pi at the interval start."""
        r = self.r
        g, P, site, chk = [], [], [], []
        for k in range(max(k0, 0), K - 1):
            # slot lineage over the interval: alive at both ends and never emptied or refilled by a newborn
            # (the first daughter keeps its slot; the division halving is excluded from the gain)
            same = r["alive"][b, :, k] & r["alive"][b, :, k + 1] & r["cont"][b, :, k + 1]
            s_, v_ = np.nonzero(same)
            if len(s_) == 0:
                continue
            a0, a1 = r["a"][b, s_, k, v_], r["a"][b, s_, k + 1, v_]
            g.append(r["gain"][b, s_, k + 1, v_] / np.maximum((a0 + a1) / 2, 1e-9))
            P.append(r["Pi"][b, s_, k, v_]); site.append(s_); chk.append(np.full(len(s_), k))
        if not g:
            return None
        return np.concatenate(g), np.concatenate(P), np.concatenate(site), np.concatenate(chk)

    def l3_loose(self, pr, rng):
        if pr is None:
            return float("nan"), 1.0, 0
        g, P, s, k = pr
        rho = _spearman(g, P)
        if not np.isfinite(rho):
            return rho, 1.0, len(g)
        grp = s * 100 + k
        order = np.argsort(grp, kind="stable")
        g_, P_, grp_ = g[order], P[order], grp[order]
        n = len(P_)
        bounds = np.flatnonzero(np.diff(grp_)) + 1
        multi = [(a, b) for a, b in zip(np.r_[0, bounds], np.r_[bounds, n]) if b - a > 1]
        # permutation draws as the loop version (rng.permutation of each segment of > 1, in order): an
        # in-place shuffle of the segment's slice of arange(n) consumes the same draws (speed change (16))
        idx = np.empty((NPERM + 1, n), np.int64)
        base = np.arange(n)
        idx[0] = base
        for i in range(1, NPERM + 1):
            row = base.copy()
            for a, b in multi:
                rng.shuffle(row[a:b])
            idx[i] = row
        # Spearman of a permutation of P = Pearson of rank(g) and the permuted rank(P); observed (row 0) and
        # null rows in one product, so equal vectors give equal values
        rg, rP = rankdata(g_), rankdata(P_)
        rg_c, rP_c = rg - rg.mean(), rP - rP.mean()
        vals = rP_c[idx] @ rg_c
        ge = int((vals[1:] >= vals[0]).sum())
        return rho, (1 + ge) / (NPERM + 1), len(g)

    def ancestor_of(self, b, site, vid, anc_ids):
        """Walk the parent table from vid up to an id in anc_ids (or -1); returns (ancestor, depth)."""
        par = self.par_t[b, site]
        d, x = 0, int(vid)
        while x >= 0 and x < MAXID and d < 4096:
            if x in anc_ids:
                return x, d
            x = int(par[x]); d += 1
        return -1, d

    def lineage_letters(self, b, s, k, v, key, after=None):
        """Enclosed letters of lineage `key` in vesicle v at checkpoint k (only in chains made after
        checkpoint `after` if given)."""
        tk = self.r["top_k"][b, s, k, v]
        tl = self.r["top_l"][b, s, k, v] if after is None else self.r["top_ln"][b, s, k, v, :, after]
        return float(tl[tk == key].sum()) if key != 0 else 0.0

    def _vid_last(self, b, s):
        """vid -> (last checkpoint alive, slot) for site s of run b."""
        r = self.r
        out = {}
        al = r["alive"][b, s]
        for k in range(N_CHECK):
            for v in np.nonzero(al[k])[0]:
                out[int(r["vid"][b, s, k, v])] = (k, v)
        return out

    def local_inherit(self, b):
        """Local inheritance test (spec 05017cd): number of sites with >= 1 vesicle alive at the last
        checkpoint that encloses >= L3_FRAC x its ancestor's dominant-key letters in chains made after
        the ancestor's last recorded checkpoint with >= L3_IDENT identity to the key (either strand)."""
        r, N, last = self.r, self.N, N_CHECK - 1
        fc = self.V.final_chains.get(b)
        if fc is None:
            return 0, {"vesicles_tested": 0}
        kt = self.V.key_table
        ok_sites, tested, depths, passed_v = 0, 0, [], 0
        for s in range(N):
            seen = self._vid_last(b, s)
            par = self.par_t[b, s]
            cs = fc["site"] == s
            c_slot, c_L, c_seq, c_born = fc["slot"][cs], fc["L"][cs].astype(np.int64), fc["seq"][cs], fc["born"][cs]
            ck = fc["ck_step"][s]
            site_ok = False
            for v in np.nonzero(r["alive"][b, s, last])[0]:
                x, d = int(r["vid"][b, s, last, v]), 0
                anc = None
                while 0 <= x < MAXID:
                    x = int(par[x]); d += 1
                    if x < 0:
                        break
                    if d >= L3_UP and x in seen and seen[x][0] < last:
                        anc = x
                        break
                if anc is None:
                    continue
                ka, va = seen[anc]
                key = int(r["top_k"][b, s, ka, va, 0])
                la = float(r["top_l"][b, s, ka, va, 0])
                if key == 0 or la <= 0 or key not in kt:
                    continue
                tested += 1; depths.append(d)
                kl, ks = kt[key]
                m = (c_slot == v) & (c_born > ck[ka])
                if not m.any():
                    continue
                n = np.minimum(c_L[m], int(kl))
                mask = np.where(n >= 64, -1, (np.int64(1) << n) - 1)
                xo = (c_seq[m] ^ np.int64(ks)) & mask
                mis = _POP8[xo.astype(np.int64).view(np.uint8).reshape(-1, 8)].sum(1)
                ident = np.maximum(n - mis, mis) / np.maximum(n, 1)            # either strand (BITS = 1: complement = NOT)
                good = ident >= L3_IDENT
                if c_L[m][good].sum() >= L3_FRAC * la:
                    site_ok = True; passed_v += 1
            ok_sites += site_ok
        return ok_sites, {"vesicles_tested": tested, "vesicles_passed": passed_v,
                          "ancestor_depth_mean": float(np.mean(depths)) if depths else None,
                          "ancestor_depth_max": int(max(depths)) if depths else None}

    def inherit(self, b, k0):
        """Persistence test from window start k0 to the last checkpoint: number of sites with a
        descendant >= 3 generations below a k0 ancestor that encloses the ancestor's dominant lineage
        in chains made after checkpoint k0, with >= 0.5 x the ancestor's letters of it at k0;
        ancestors' z and w (descendants alive at the last checkpoint)."""
        r, N, K = self.r, self.N, N_CHECK
        last = K - 1
        z = r["top_l"][b].sum(-1) / np.maximum(r["phi"][b], 1e-9)       # lineage letters (2 dominant keys) per volume share
        ok_sites, anc_z, anc_w, anc_list = 0, [], [], []
        for s in range(N):
            av = np.nonzero(r["alive"][b, s, k0])[0]
            anc = {int(r["vid"][b, s, k0, v]): v for v in av}
            if not anc:
                continue
            lv = np.nonzero(r["alive"][b, s, last])[0]
            w = {a: 0 for a in anc}
            passed = False
            for v in lv:
                a, d = self.ancestor_of(b, s, int(r["vid"][b, s, last, v]), set(anc))
                if a < 0:
                    continue
                w[a] += 1
                va = anc[a]
                key = int(r["top_k"][b, s, k0, va, 0])
                la = float(r["top_l"][b, s, k0, va, 0])
                gd = int(r["gen"][b, s, last, v]) - int(r["gen"][b, s, k0, va])
                if key != 0 and gd >= 3 and self.lineage_letters(b, s, last, v, key, after=k0) >= 0.5 * la and la > 0:
                    passed = True
            ok_sites += passed
            for a, va in anc.items():
                anc_z.append(z[s, k0, va]); anc_w.append(w[a]); anc_list.append((s, a))
        return ok_sites, anc_z, anc_w, z

    def l3_strict_and_l4(self, b, k0, rng):
        r, N, K = self.r, self.N, N_CHECK
        last = K - 1
        _, anc_z, anc_w, z = self.inherit(b, k0)          # window-start ancestors: L4 selection (w, z)
        ok_sites, loc = self.local_inherit(b)             # L3 strict inheritance (local, spec 05017cd)
        strict_inh = ok_sites >= 0.5 * N
        # heritability: parent z at k, daughter z at k+1 (daughters born in interval k -> k+1)
        ev = self.div
        sel = (ev[:, 0] // N == b) & (ev[:, 4] > k0) & (ev[:, 4] <= last)
        zp, zd = [], []
        for row, par, d1, d2, k in ev[sel]:
            # a division logged in interval k happened after checkpoint k-1 and before checkpoint k
            s = int(row % N)
            if k < 1:
                continue
            pv = np.nonzero(r["alive"][b, s, k - 1] & (r["vid"][b, s, k - 1] == par))[0]
            if len(pv) == 0:
                continue
            for dd in (d1, d2):
                if dd < 0:
                    continue
                dv = np.nonzero(r["alive"][b, s, k] & (r["vid"][b, s, k] == dd))[0]
                if len(dv):
                    zp.append(z[s, k - 1, pv[0]]); zd.append(z[s, k, dv[0]])
        herit = None
        if len(zp) >= 3 and np.std(zp) > 0:
            lr = linregress(zp, zd)
            herit = (float(lr.slope), float(lr.pvalue))
        anc_z, anc_w = np.array(anc_z), np.array(anc_w, float)
        rho_w = _spearman(anc_w, anc_z) if len(anc_z) >= 3 else float("nan")
        thr = float("nan")
        if np.isfinite(rho_w) and len(anc_z) >= 10:
            null = _spearman_null(anc_w, anc_z, rng)
            null = np.array([x for x in null if np.isfinite(x)])
            thr = float(np.percentile(null, 99)) if len(null) else float("nan")
        sel_ok = np.isfinite(rho_w) and np.isfinite(thr) and len(anc_z) >= 10 and rho_w > thr
        her_ok = herit is not None and herit[0] > 0 and herit[1] < 0.01
        return {"inherit_sites": ok_sites, "inherit_share": ok_sites / N, "inherit_ok": bool(strict_inh), "inherit_local": loc, "heritability": herit, "n_division_pairs": len(zp),
                "n_ancestors": int(len(anc_z)), "rho_wz": rho_w, "rho_wz_null99": thr, "selection_ok": bool(sel_ok),
                "heritability_ok": bool(her_ok)}

    def key_risers(self, b, k0):
        r = self.r
        last = N_CHECK - 1
        def freq(k):
            al = r["alive"][b, :, k]
            ks = r["top_k"][b, :, k, :, 0][al]
            ks = ks[ks != 0]
            u, c = np.unique(ks, return_counts=True)
            return dict(zip(u.tolist(), (c / max(al.sum(), 1)).tolist()))
        f0, f1 = freq(k0), freq(last)
        d = {k: f1.get(k, 0) - f0.get(k, 0) for k in set(f0) | set(f1)}
        top = sorted(d.items(), key=lambda x: -x[1])[:3]
        return [(int(k), round(v, 4), self.V.key_table.get(k)) for k, v in top]

    # ---- per-run detection ----
    def _detect_one(self, b, k0):
        """Detectors of run b with window start k0 (numpy / scipy only: safe in a forked worker)."""
        rng = np.random.default_rng(1000 + b)
        pr = self.pairs(b, k0)
        npairs = 0 if pr is None else len(pr[0])
        gens = self.r["gen"][b, :, k0:][self.r["alive"][b, :, k0:]]
        span = int(gens.max() - gens.min()) if len(gens) else 0
        ev_ok = npairs >= L3_PAIRS and span >= L3_GENS
        rho, p, _ = self.l3_loose(pr, rng)
        loose = ev_ok and np.isfinite(rho) and rho >= L3_RHO and p < L3_P
        d = self.l3_strict_and_l4(b, k0, rng)
        strict = loose and d["inherit_ok"]
        l4 = strict and d["selection_ok"] and d["heritability_ok"]
        d["key_risers"] = self.key_risers(b, k0)
        return npairs, span, ev_ok, rho, p, loose, strict, l4, d

    def _run_tasks(self, tasks):
        """Evaluate (b, k0) tasks, in DET_WORKERS forked processes when there are several (each run
        has its own RNG, so the results do not depend on the split); serial fallback on any error."""
        global _DET
        if DET_WORKERS > 1 and len(tasks) > 1:
            import multiprocessing as mp
            _DET = self
            try:
                with mp.get_context("fork").Pool(min(DET_WORKERS, len(tasks))) as pool:
                    res = pool.map_async(_det_worker, tasks, chunksize=1).get(DET_TIMEOUT)
                return [(b, k0, v) for (b, k0), v in zip(tasks, res)]
            except Exception as e:                                   # noqa: BLE001
                print("detector pool failed, serial:", repr(e)[:200], flush=True)
            finally:
                _DET = None
        return [(b, k0, self._detect_one(b, k0)) for b, k0 in tasks]

    def detect(self, window_override=None):
        B = self.B
        # candidates: E1, L2_loose and copies >= L_H (L2_strict needs the paired rerun: the runner gates on it)
        l2long = (self.out["E1"] & self.out["L2_candidate"] & (self.out["copy_mean_len"] >= L_H)).numpy() \
            if window_override is None else np.ones(B, bool)
        e3s = self.res["E3_strict"].numpy() if window_override is None else np.ones(B, bool)
        tasks = []
        for b in range(B):
            if window_override is not None:
                k0 = window_override
            else:
                if not (l2long[b] and e3s[b]):
                    continue
                k0 = self.joint_window(b)
                if k0 is None:
                    continue
            tasks.append((b, k0))
        return self._finish(tasks)

    def detect_at(self, k0s):
        """Detectors without the L2/E3 gate (paired controls, spec a9227f6 choice 19): run b is scored
        over the window starting at k0s[b] (None: not scored); evaluability (pairs, generations) as usual."""
        return self._finish([(b, int(k)) for b, k in enumerate(k0s) if k is not None])

    def _finish(self, tasks):
        B = self.B
        o = {k: np.zeros(B, bool) for k in ("L3_loose", "L3_strict", "L4", "L3_evaluable")}
        o.update({k: np.full(B, np.nan) for k in ("rho_g_Pi", "p_g_Pi", "rho_w_z", "n_pairs", "gen_span", "joint_window_k0")})
        details = [None] * B
        for b, k0, v in self._run_tasks(tasks):
            npairs, span, ev_ok, rho, p, loose, strict, l4, d = v
            o["n_pairs"][b], o["gen_span"][b], o["joint_window_k0"][b] = npairs, span, k0
            o["L3_evaluable"][b] = ev_ok
            o["rho_g_Pi"][b], o["p_g_Pi"][b] = rho, p
            o["L3_loose"][b], o["L3_strict"][b], o["L4"][b] = loose, strict, l4
            o["rho_w_z"][b] = d["rho_wz"]
            details[b] = d
        out = {k: torch.as_tensor(v) for k, v in o.items()}
        out["L3_details"] = details
        out["lineage_overflow"] = torch.as_tensor(self.ovf)
        out["lineage"] = self                                          # host record (validation / point study)
        return out


# the 4c-5 simulate: chem4c4's, with the recording Vesicles class switched in for this call only
# (chem4c4.simulate itself stays unpatched; paired_control / paired_control_l2 rerun without record)
def simulate(*a, window_override=None, record_only=False, **k):
    old = _c.Vesicles
    _c.Vesicles = Vesicles
    Lineage._override = window_override
    Lineage._record_only = record_only
    try:
        return _c.simulate(*a, **k)
    finally:
        _c.Vesicles = old
        Lineage._override = None
        Lineage._record_only = False


def paired_inherit(u, seed, fidelity, out, override=None, window_override=None, mask=None, reuse_main=False):
    """L3 strict, copied not re-enclosed (spec §2, revised bd0f418). For the runs in `mask` (default:
    L3 loose with the unpaired inheritance test passed) the conditions are rerun as a batch of their
    own, once as they are (main') and with copy_random (c), both with seed `seed`, so they differ in
    the copied letters only (as paired_control_l2). Both are scored with the persistence test over
    the main run's window (k0 to the last checkpoint). Paired condition: share(main') - share(c) >=
    L3_RC_MARGIN. reuse_main: `out` itself is main' (all B runs in the mask, `out` from this seed,
    override and window). Returns L3_strict and L4 of `out` with the condition applied, and the
    shares of main' and c (nan outside the mask)."""
    B = u.shape[0]
    det = out["L3_details"]
    if mask is None:
        mask = [bool(out["L3_loose"][b]) and det[b] is not None and det[b]["inherit_ok"] for b in range(B)]
    idx = [b for b in range(B) if mask[b] and det[b] is not None]
    sm, sc = np.full(B, np.nan), np.full(B, np.nan)
    if idx:
        def rerun(rnd):
            g = torch.Generator(device=u.device)
            g.manual_seed(seed)
            ov = dict(override or {})
            if rnd:
                ov["copy_random"] = True
            return simulate(u[idx], g, fidelity, override=ov, window_override=window_override, record_only=True)
        full = reuse_main and len(idx) == B
        m = out if full else rerun(False)
        c = rerun(True)
        for j, b in enumerate(idx):
            k0 = int(out["joint_window_k0"][b])
            lm, lc = m["lineage"], c["lineage"]
            sm[b] = lm.local_inherit(b if full else j)[0] / lm.N
            sc[b] = lc.local_inherit(j)[0] / lc.N
        if PERTURB_SETS:                     # Review 9: detector perturbation on the same records (no simulation)
            global PERTURB_OUT
            PERTURB_OUT = {}
            lo = out["lineage"]
            for name, cs in PERTURB_SETS:
                saved = {k_: globals()[k_] for k_ in cs}
                globals().update(cs)
                try:
                    a0, a1, a2 = np.full(B, np.nan), np.full(B, np.nan), np.full(B, np.nan)
                    for j, b in enumerate(idx):
                        a0[b] = lo.local_inherit(b)[0] / lo.N
                        a1[b] = lm.local_inherit(b if full else j)[0] / lm.N
                        a2[b] = lc.local_inherit(j)[0] / lc.N
                finally:
                    globals().update(saved)
                PERTURB_OUT[name] = (a0, a1, a2)
        del c
    if PERTURB_SETS and not idx:
        globals()["PERTURB_OUT"] = {}
    ok = torch.as_tensor(np.nan_to_num(sm - sc, nan=-1.0) >= L3_RC_MARGIN - 1e-9)
    return out["L3_strict"] & ok, out["L4"] & ok, sm, sc


def control_k0(out):
    """Window start for scoring a paired control: the main run's joint window, else the whole record."""
    k = out["joint_window_k0"].numpy()
    return [int(x) if np.isfinite(x) else 0 for x in k]


def score_control(u, seed, fidelity, override, k0s):
    """Paired control (beta = 0 or random copy; same constants and seed) scored WITHOUT the L2/E3 gate
    wherever evaluable (spec a9227f6, choice 19), over the window starting at k0s[b]; L3 strict / L4
    with the paired random-copy condition. Returns numpy arrays, with the control's own gate status."""
    g = torch.Generator(device=u.device)
    g.manual_seed(seed)
    oc = simulate(u, g, fidelity, override=dict(override), record_only=True)
    oc.update(oc["lineage"].detect_at(k0s))
    st, l4, sm, sc = paired_inherit(u, seed, fidelity, oc, override=dict(override))
    n = lambda k: oc[k].numpy() if torch.is_tensor(oc[k]) else np.asarray(oc[k])
    res = {"L3_evaluable": n("L3_evaluable"), "L3_loose": n("L3_loose"), "L3_strict": st.numpy(), "L4": l4.numpy(),
           "rho_g_Pi": n("rho_g_Pi"), "p_g_Pi": n("p_g_Pi"), "inherit_share": sm, "inherit_share_rc": sc,
           "E1": n("E1"), "L2_candidate": n("L2_candidate"), "E3_strict": n("E3_strict"),
           "copy_len_ok": (oc["copy_mean_len"] >= L_H).numpy()}
    del oc
    return res


paired_control = _c.paired_control
paired_control_l2 = _c.paired_control_l2
