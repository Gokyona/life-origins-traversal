"""Stage 4c-3 copy of cat4c2.py (md5 13a08af5) with deterministic tie-breaking in abundant():
the species sort is stable and the top-nx score has the species key as a final tie-breaker, so
the species chosen when counts and lengths tie no longer depend on the GPU's tie order (the
kill-and-resume test of 4c-3 showed different results across processes in copying runs).
--- 4c-2 description follows ---
Stage 4c-2 catalysis machinery used by chem4c2.simulate (spec docs/stage4c2_spec.md).

Catalysis acts among the abundant species of a site: every refresh, each site's NX most
abundant distinct sequences with at least MIN_COUNT copies form its species set X. Catalysed
ligations are all triples (C, A, B) with C, A in X and B in X or an activated monomer letter
of A's chemistry, same chemistry for A and B, L_A + L_B <= L_MAX, and hash(C, A, B) < p_cat(L_C)
(raf4c2.catalyses). Sequences present as a single copy neither catalyse nor get catalytically
ligated (a stated approximation; they can still grow, break and be produced as products).

Rates (mass action, per hour): chain + chain  r_cc * n_C * n_A * n_B,
                               chain + monomer r_cm * n_C * n_A * M*_letter,
with r_cm = k_c * f_min * kappa * Arr_cat(T) * g_dry * c_eff^2 / V^2 and r_cc = r_cm * phi,
phi = M* / (M + M*) (share of activated ends). Events per step are Poisson, thinned to the
instances available, and assigned to distinct slots; A's slot receives AB (sequence of B
appended after A), B's slot is freed (or one activated monomer is used).
"""
import torch
import gk_raf as raf
import words as wd

NX = 16
MIN_COUNT = 2
I64MAX = (1 << 63) - 1


def _words(seq, L):
    """Accept a plain low-word value (same shape as L) or a word tensor [..., W]."""
    return wd.from_low(seq) if seq.dim() == L.dim() else seq


def species_key(chem, L, seq):
    seq = _words(seq, L)
    z = raf._mix(chem + raf._s64(0x5BD1E995))
    z = raf._mix(z ^ L)
    z = raf._mix(z ^ seq[..., 0])
    for w in range(1, seq.shape[-1]):
        z = raf._mix(z ^ seq[..., w])
    return z


def fold(seq):
    """64-bit value of a word sequence for the catalysis hash: the sequence itself when W = 1."""
    z = seq[..., 0]
    for w in range(1, seq.shape[-1]):
        z = raf._mix(z ^ seq[..., w])
    return z


def _mask(nbits):
    return wd.mask(nbits)


def abundant(L, seq, chem, nx=NX, min_count=MIN_COUNT):
    """For rows [Ra,S]: the nx most abundant species with >= min_count copies.
    Returns dict of [Ra,nx] tensors (valid, key, chem, L, seq, count) and sidx [Ra,S]
    (species index of each slot in X, -1 if none)."""
    Ra, S = L.shape
    dev = L.device
    occ = L > 0
    k = torch.where(occ, species_key(chem, L, seq), torch.full_like(L, I64MAX))
    ks, order = torch.sort(k, dim=1, stable=True)
    pos = torch.arange(S, device=dev).expand(Ra, S)
    boundary = torch.ones_like(ks, dtype=torch.bool)
    boundary[:, 1:] = ks[:, 1:] != ks[:, :-1]
    start = torch.cummax(torch.where(boundary, pos, torch.zeros_like(pos)), dim=1).values
    cnt = torch.zeros(Ra, S, dtype=torch.long, device=dev).scatter_add_(1, start, torch.ones_like(start))
    Ls = torch.gather(L, 1, order)
    ok = boundary & (ks != I64MAX) & (cnt >= min_count)
    # (count, length) first, then the low 40 bits of the species key: practically no ties
    score = torch.where(ok, (cnt * 128 + Ls) * (1 << 40) + (ks & ((1 << 40) - 1)), torch.full_like(cnt, -1))
    nx = min(nx, S)
    val, idx = torch.topk(score, nx, dim=1)
    valid = val >= 0
    rep = torch.gather(order, 1, idx)
    X = {"valid": valid, "key": torch.gather(ks, 1, idx), "count": torch.where(valid, torch.gather(cnt, 1, idx), 0),
         "L": torch.where(valid, torch.gather(L, 1, rep), 0), "seq": torch.gather(seq, 1, rep[..., None].expand(-1, -1, seq.shape[-1])),
         "chem": torch.gather(chem, 1, rep)}
    rank_at = torch.full((Ra, S), -1, dtype=torch.long, device=dev)
    ranks = torch.arange(nx, device=dev).expand(Ra, nx)
    rank_at.scatter_(1, idx, torch.where(valid, ranks, torch.full_like(ranks, -1)))
    sidx_sorted = torch.gather(rank_at, 1, start)
    sidx_sorted = torch.where(ks != I64MAX, sidx_sorted, torch.full_like(sidx_sorted, -1))
    sidx = torch.full((Ra, S), -1, dtype=torch.long, device=dev).scatter_(1, order, sidx_sorted)
    X["sidx"] = sidx
    return X


def rematch(X, L, seq, chem):
    """Between refreshes: keep the species list of X but recount copies and slot positions on
    the current state (species that vanished get count 0 and drop out of the rates)."""
    occ = L > 0
    k = torch.where(occ, species_key(chem, L, seq), torch.full_like(L, I64MAX))
    eq = (k.unsqueeze(2) == X["key"].unsqueeze(1)) & X["valid"].unsqueeze(1)     # [Ra,S,nx]
    has = eq.any(2)
    sidx = torch.where(has, eq.float().argmax(2), torch.full_like(L, -1))
    count = eq.sum(1)
    Y = dict(X)
    Y["count"] = count
    Y["sidx"] = sidx
    return Y


def catalytic_triples(X, p0, L_cat, lam, salt, alpha, bits, L_max, chunk=4096):
    """All catalysed triples among X. p0, lam float64 [Ra]; L_cat long [Ra].
    Returns long tensors (row, c, a, b) with b >= nx meaning monomer letter b - nx."""
    Ra, nx = X["L"].shape
    dev = X["L"].device
    nb = nx + alpha
    out = []
    ci = torch.arange(nx, device=dev).view(nx, 1, 1)
    ai = torch.arange(nx, device=dev).view(1, nx, 1)
    bi = torch.arange(nb, device=dev).view(1, 1, nb)
    live = ((X["valid"] & (X["L"] >= L_cat.view(-1, 1))).any(1)).nonzero().squeeze(1)   # rows with a catalyst
    for r0 in range(0, live.numel(), chunk):
        sl = live[r0:r0 + chunk]
        v, L, sq, ch = X["valid"][sl], X["L"][sl], fold(X["seq"][sl]), X["chem"][sl]
        n = v.shape[0]
        Cv = (v & (L >= L_cat[sl].view(n, 1))).view(n, nx, 1, 1)
        Av = v.view(n, 1, nx, 1)
        is_mono = (bi >= nx).view(1, 1, 1, nb)
        bspec = bi.clamp(max=nx - 1).view(nb)
        chA = ch.view(n, 1, nx, 1)
        LA = L.view(n, 1, nx, 1)
        chB = torch.where(is_mono, chA, ch[:, bspec].view(n, 1, 1, nb))
        LB = torch.where(is_mono, torch.ones_like(LA), L[:, bspec].view(n, 1, 1, nb))
        Bv = torch.where(is_mono, torch.ones_like(Av), v[:, bspec].view(n, 1, 1, nb) & (chB == chA))
        ok = Cv & Av & Bv & (LA + LB <= L_max)
        if not bool(ok.any()):
            continue
        r_, c_, a_, b_ = ok.nonzero(as_tuple=True)
        C = (ch[r_, c_], L[r_, c_], sq[r_, c_])
        A = (ch[r_, a_], L[r_, a_], sq[r_, a_])
        mono = b_ >= nx
        bs = b_.clamp(max=nx - 1)
        B = (ch[r_, a_], torch.where(mono, torch.ones_like(b_), L[r_, bs]),
             torch.where(mono, fold(wd.from_low(b_ - nx)), sq[r_, bs]))
        rr = sl[r_]
        hit = raf.catalyses(C, A, B, p0[rr], L_cat[rr], lam[rr], salt)
        out.append(torch.stack([rr[hit], c_[hit], a_[hit], b_[hit]], 1))
    if not out:
        return torch.zeros(0, 4, dtype=torch.long, device=dev)
    return torch.cat(out)


def apply_catalysis(L, seq, chem, Ms, rows, X, trip, r_cm, r_cc, dt, gen, aux, poisson, alpha, bits, L_max, tmult=None):
    """Fire catalysed ligations for the active rows (indices into the full state).
    r_cm, r_cc, dt: float32 [Ra]. Modifies L, seq, chem, Ms in place. Returns events per
    active row and the number of rows whose demand was thinned (saturation)."""
    Ra, nx = X["L"].shape
    dev = L.device
    ev_row = torch.zeros(Ra, device=dev)
    if trip.numel() == 0:
        return ev_row, 0
    ri, c, a, b = trip.unbind(1)
    mono = b >= nx
    bs = b.clamp(max=nx - 1)
    nC = X["count"][ri, c].float()
    nA = X["count"][ri, a].float()
    chA = X["chem"][ri, a]
    ms_letter = Ms[rows[ri], chA] / alpha
    nB = torch.where(mono, ms_letter, X["count"][ri, bs].float())
    lam = torch.where(mono, r_cm[ri], r_cc[ri]) * nC * nA * nB * dt[ri]
    if tmult is not None:                                   # Gokyona local_conc: per-triple compartment factor
        lam = lam * tmult
    k = poisson(lam, gen, aux)
    # thin to available instances: species demand (A and species-B roles) and monomers
    gA = ri * nx + a
    gB = ri * nx + bs
    dem = torch.zeros(Ra * nx, device=dev).index_add_(0, gA, k)
    dem = dem.index_add_(0, gB[~mono], k[~mono])
    avail = X["count"].reshape(-1).float()
    fac = torch.clamp(avail / torch.clamp(dem, min=1e-9), max=1.0)
    gm = rows[ri] * 2 + chA
    dem_m = torch.zeros(Ms.numel(), device=dev).index_add_(0, gm[mono], k[mono])
    fac_m = torch.clamp(Ms.reshape(-1) / torch.clamp(dem_m, min=1e-9), max=1.0)
    f = torch.minimum(fac[gA], torch.where(mono, fac_m[gm], fac[gB]))
    saturated = int(((f < 1) & (k > 0)).sum())
    k = torch.floor(k * f).long()
    if int(k.sum()) == 0:
        return ev_row, saturated
    t_ev = torch.repeat_interleave(torch.arange(k.numel(), device=dev), k)
    E = t_ev.numel()
    # instance ranks: uses (A role for every event, B role for species-B events) grouped by species
    use_g = torch.cat([gA[t_ev], gB[t_ev][~mono[t_ev]]])
    use_e = torch.cat([torch.arange(E, device=dev), torch.arange(E, device=dev)[~mono[t_ev]]])
    use_role = torch.cat([torch.zeros(E, dtype=torch.long, device=dev),
                          torch.ones(int((~mono[t_ev]).sum()), dtype=torch.long, device=dev)])
    o = torch.argsort(use_g, stable=True)
    gsorted = use_g[o]
    first = torch.searchsorted(gsorted, gsorted, right=False)
    rank = torch.empty_like(o)
    rank[o] = torch.arange(o.numel(), device=dev) - first
    # slots of each species, grouped by g = row * nx + sidx
    S = L.shape[1]
    sidx = X["sidx"]
    rr, ss = (sidx >= 0).nonzero(as_tuple=True)
    gslot = rr * nx + sidx[rr, ss]
    os_ = torch.argsort(gslot, stable=True)
    gslot_s = gslot[os_]
    slot_flat = rows[rr[os_]] * S + ss[os_]
    start = torch.searchsorted(gslot_s, use_g)
    slot_of_use = slot_flat[start + rank]
    sA = torch.full((E,), -1, dtype=torch.long, device=dev)
    sB = torch.full((E,), -1, dtype=torch.long, device=dev)
    roleA = use_role == 0
    sA[use_e[roleA]] = slot_of_use[roleA]
    sB[use_e[~roleA]] = slot_of_use[~roleA]
    Lf, sf, cf = L.view(-1), seq.view(-1, seq.shape[-1]), chem.view(-1)
    LA = Lf[sA]
    is_m = mono[t_ev]
    letter = (b[t_ev] - nx).clamp(min=0)
    LB = torch.where(is_m, torch.ones_like(LA), Lf[sB.clamp(min=0)])
    seqB = torch.where(is_m[:, None], wd.from_low(letter), sf[sB.clamp(min=0)])
    newseq = sf[sA] | wd.shl(seqB & _mask(LB * bits), LA * bits)
    keepB = ~is_m
    sBk = sB[keepB]
    Lf[sBk] = 0
    sf[sBk] = 0
    sf[sA] = newseq
    Lf[sA] = LA + LB
    gmE = gm[t_ev][is_m]
    Ms.view(-1).index_add_(0, gmE, -torch.ones(gmE.numel(), device=dev))
    ev_row.index_add_(0, ri[t_ev], torch.ones(E, device=dev))
    return ev_row, saturated


HASH_DEVICE = None        # set to the GPU by chem4c2 so the RAF hashing runs there


def site_raf(Xr, p0, L_cat, lam, salt, alpha, bits):
    """maxRAF of one site's species set (CPU tensors of one row of X). Ligations only.
    Returns dict: size, components, member keys (non-food species in the maxRAF), catalyst
    lengths, nucleotide share of member mass."""
    v = Xr["valid"]
    present = set(int(x) for x in Xr["key"][v])
    if int(v.sum()) == 0:
        return {"size": 0, "present_keys": present}
    L, sq, ch, key = Xr["L"][v], Xr["seq"][v], Xr["chem"][v], Xr["key"][v]
    # long members (> 64 letters, multi-word engine only) are not split in the RAF analysis
    nxs = L.numel()
    # food: monomer letters of both chemistries
    fch = torch.tensor([c for c in (0, 1) for _ in range(alpha)], dtype=torch.long)
    fL = torch.ones(2 * alpha, dtype=torch.long)
    fsq = wd.from_low(torch.tensor([x for _ in (0, 1) for x in range(alpha)], dtype=torch.long))
    allL = torch.cat([L, fL]); allsq = torch.cat([sq, fsq]); allch = torch.cat([ch, fch])
    allkey = torch.cat([key, species_key(fch, fL, fsq)])
    food = torch.zeros(nxs + 2 * alpha, dtype=torch.bool)
    food[nxs:] = True
    # every split of every member with length >= 2
    P = ((L >= 2) & (L <= 64)).nonzero().squeeze(1)
    if P.numel() == 0:
        return {"size": 0, "present_keys": present}
    Lp = L[P]
    jmax = int(Lp.max())
    j = torch.arange(1, jmax).view(1, -1)
    okj = j < Lp.view(-1, 1)
    pr, jj = okj.nonzero(as_tuple=True)
    pidx = P[pr]
    jv = j[0, jj]
    pre_key = species_key(ch[pidx], jv, sq[pidx] & _mask(jv * bits))
    suf_key = species_key(ch[pidx], L[pidx] - jv, wd.shr(sq[pidx], jv * bits) & _mask((L[pidx] - jv) * bits))
    match = lambda kk: (kk.view(-1, 1) == allkey.view(1, -1))
    ma, mb = match(pre_key), match(suf_key)
    has = ma.any(1) & mb.any(1)
    if not bool(has.any()):
        return {"size": 0, "present_keys": present}
    ra = ma[has].float().argmax(1)
    rb = mb[has].float().argmax(1)
    rab = pidx[has]
    R = ra.numel()
    cands = ((L >= L_cat)).nonzero().squeeze(1)
    if cands.numel() == 0:
        return {"size": 0, "present_keys": present}
    cr = torch.arange(R).repeat_interleave(cands.numel())
    cc = cands.repeat(R)
    hd = HASH_DEVICE or "cpu"
    d = lambda t: t.to(hd)
    fq, fa = fold(sq), fold(allsq)
    hit = raf.catalyses((d(ch[cc]), d(L[cc]), d(fq[cc])), (d(allch[ra[cr]]), d(allL[ra[cr]]), d(fa[ra[cr]])),
                        (d(allch[rb[cr]]), d(allL[rb[cr]]), d(fa[rb[cr]])),
                        torch.tensor(float(p0), dtype=torch.float64, device=hd), int(L_cat),
                        torch.tensor(float(lam), dtype=torch.float64, device=hd), salt).cpu()
    in_raf = raf.max_raf(nxs + 2 * alpha, food, ra, rb, rab, cr[hit], cc[hit])
    # catalysed reactions of the maxRAF as key tuples (reactant, reactant, product, catalyst)
    hr, hc = cr[hit], cc[hit]
    keep = in_raf[hr]
    cat_tuples = set(zip(allkey[ra[hr[keep]]].tolist(), allkey[rb[hr[keep]]].tolist(),
                         allkey[rab[hr[keep]]].tolist(), key[hc[keep]].tolist()))
    size = int(in_raf.sum())
    if size == 0:
        return {"size": 0, "present_keys": present}
    members = torch.unique(torch.cat([ra[in_raf], rb[in_raf], rab[in_raf]]))
    members = members[members < nxs]
    used_cat = torch.unique(cc[hit][in_raf[cr[hit]]])
    mass = L[members].float()
    return {"size": size, "components": raf.raf_components(ra, rb, rab, in_raf),
            "member_keys": set(int(x) for x in key[members]),
            "member_counts": {int(key[m]): int(Xr["count"][v][m]) for m in members.tolist()},
            "product_counts": {int(key[m]): int(Xr["count"][v][m]) for m in torch.unique(rab[in_raf]).tolist()},
            "member_L": [int(x) for x in L[members]], "member_chem": [int(x) for x in ch[members]],
            "present_keys": present, "cat_tuples": cat_tuples,
            "catalyst_L": [int(x) for x in L[used_cat]],
            "nuc_share": float((mass * (ch[members] == 0).float()).sum() / mass.sum())}


def present_keys(L, seq, chem, min_count=MIN_COUNT):
    """Per row: set of species keys with >= min_count copies (no cap). L etc. [Ra,S] -> list of sets."""
    Ra, S = L.shape
    occ = L > 0
    k = torch.where(occ, species_key(chem, L, seq), torch.full_like(L, I64MAX))
    rows = torch.arange(Ra, device=L.device).view(Ra, 1).expand(Ra, S)
    pair = torch.stack([rows[occ], k[occ]], 1)
    if pair.numel() == 0:
        return [set() for _ in range(Ra)]
    uk, cnt = torch.unique(pair, dim=0, return_counts=True)
    uk = uk[cnt >= min_count].cpu()
    out = [set() for _ in range(Ra)]
    for r, kk in uk.tolist():
        out[r].add(kk)
    return out


def subset(Xr, keys):
    """Row X restricted to species whose key is in `keys`."""
    v = Xr["valid"] & torch.tensor([int(k) in keys for k in Xr["key"].tolist()], dtype=torch.bool)
    return dict(Xr, valid=v)


def tuples_raf(tuples, alpha):
    """maxRAF of a reaction set given as (a, b, ab, catalyst) key tuples; food = monomer letters
    of both chemistries. Returns (number of reactions in the maxRAF, the RAF's tuples)."""
    if not tuples:
        return 0, set()
    fch = torch.tensor([c for c in (0, 1) for _ in range(alpha)], dtype=torch.long)
    food_keys = set(species_key(fch, torch.ones_like(fch), torch.tensor([x for _ in (0, 1) for x in range(alpha)])).tolist())
    rx = sorted({(a, b, ab) for a, b, ab, _ in tuples})
    keys = sorted({k for t in tuples for k in t} | food_keys)
    ix = {k: i for i, k in enumerate(keys)}
    rid = {r: i for i, r in enumerate(rx)}
    t = lambda z: torch.tensor(z, dtype=torch.long)
    ra, rb, rab = t([ix[r[0]] for r in rx]), t([ix[r[1]] for r in rx]), t([ix[r[2]] for r in rx])
    cr = t([rid[(a, b, ab)] for a, b, ab, _ in tuples])
    cc = t([ix[c] for *_, c in tuples])
    food = torch.zeros(len(keys), dtype=torch.bool)
    for k in food_keys:
        food[ix[k]] = True
    in_raf = raf.max_raf(len(keys), food, ra, rb, rab, cr, cc)
    kept = {rx[i] for i in in_raf.nonzero().squeeze(1).tolist()}
    return len(kept), {tp for tp in tuples if tp[:3] in kept}
