"""Stage 4c-4: compartment-aware catalysis for chem4c4 (spec docs/stage4c4_spec.md §3: catalysed
ligations happen only between chains with the same compartment label). cat4c3.py (md5 8dba0b4f)
is unchanged and still used for everything else (detectors, RAF analysis). The three functions
below are copies of cat4c3.abundant, catalytic_triples and rematch with one change each: the
species identity includes the compartment label (so the same sequence inside vesicle j and
outside are different species), and a catalysed triple needs catalyst, A and a species-B in the
same compartment (an activated-monomer B can be used in any compartment). The catalysis hash
still uses the true chemistry, length and sequence, so which reactions are catalysed is the
same chemistry in every compartment.
"""
import torch
import gk_raf as raf
import words as wd
from gk_cat3 import NX, MIN_COUNT, I64MAX, species_key, _mask, fold  # noqa: F401

LAB_SALT = raf._s64(0x9E3779B97F4A7C15)


def lab_key(chem, L, seq, lab):
    return raf._mix(species_key(chem, L, seq) ^ (lab.long() * LAB_SALT))


def abundant_lab(L, seq, chem, lab, nx=NX, min_count=MIN_COUNT):
    """cat4c3.abundant with the compartment label in the species identity; X["lab"] added."""
    Ra, S = L.shape
    dev = L.device
    occ = L > 0
    k = torch.where(occ, lab_key(chem, L, seq, lab), torch.full_like(L, I64MAX))
    ks, order = torch.sort(k, dim=1, stable=True)
    pos = torch.arange(S, device=dev).expand(Ra, S)
    boundary = torch.ones_like(ks, dtype=torch.bool)
    boundary[:, 1:] = ks[:, 1:] != ks[:, :-1]
    start = torch.cummax(torch.where(boundary, pos, torch.zeros_like(pos)), dim=1).values
    cnt = torch.zeros(Ra, S, dtype=torch.long, device=dev).scatter_add_(1, start, torch.ones_like(start))
    Ls = torch.gather(L, 1, order)
    ok = boundary & (ks != I64MAX) & (cnt >= min_count)
    score = torch.where(ok, (cnt * 128 + Ls) * (1 << 40) + (ks & ((1 << 40) - 1)), torch.full_like(cnt, -1))
    nx = min(nx, S)
    val, idx = torch.topk(score, nx, dim=1)
    valid = val >= 0
    rep = torch.gather(order, 1, idx)
    X = {"valid": valid, "key": torch.gather(ks, 1, idx), "count": torch.where(valid, torch.gather(cnt, 1, idx), 0),
         "L": torch.where(valid, torch.gather(L, 1, rep), 0), "seq": torch.gather(seq, 1, rep[..., None].expand(-1, -1, seq.shape[-1])),
         "chem": torch.gather(chem, 1, rep), "lab": torch.gather(lab.long(), 1, rep)}
    rank_at = torch.full((Ra, S), -1, dtype=torch.long, device=dev)
    ranks = torch.arange(nx, device=dev).expand(Ra, nx)
    rank_at.scatter_(1, idx, torch.where(valid, ranks, torch.full_like(ranks, -1)))
    sidx_sorted = torch.gather(rank_at, 1, start)
    sidx_sorted = torch.where(ks != I64MAX, sidx_sorted, torch.full_like(sidx_sorted, -1))
    sidx = torch.full((Ra, S), -1, dtype=torch.long, device=dev).scatter_(1, order, sidx_sorted)
    X["sidx"] = sidx
    return X


def rematch_lab(X, L, seq, chem, lab):
    """cat4c3.rematch with the compartment label in the species identity."""
    occ = L > 0
    k = torch.where(occ, lab_key(chem, L, seq, lab), torch.full_like(L, I64MAX))
    eq = (k.unsqueeze(2) == X["key"].unsqueeze(1)) & X["valid"].unsqueeze(1)
    has = eq.any(2)
    sidx = torch.where(has, eq.float().argmax(2), torch.full_like(L, -1))
    Y = dict(X)
    Y["count"] = eq.sum(1)
    Y["sidx"] = sidx
    return Y


def catalytic_triples_lab(X, p0, L_cat, lam, salt, alpha, bits, L_max, chunk=4096):
    """cat4c3.catalytic_triples restricted to catalyst, A and species-B of one compartment."""
    Ra, nx = X["L"].shape
    dev = X["L"].device
    nb = nx + alpha
    out = []
    bi = torch.arange(nb, device=dev).view(1, 1, nb)
    live = ((X["valid"] & (X["L"] >= L_cat.view(-1, 1))).any(1)).nonzero().squeeze(1)
    for r0 in range(0, live.numel(), chunk):
        sl = live[r0:r0 + chunk]
        v, L, sq, ch, lb = X["valid"][sl], X["L"][sl], fold(X["seq"][sl]), X["chem"][sl], X["lab"][sl]
        n = v.shape[0]
        Cv = (v & (L >= L_cat[sl].view(n, 1))).view(n, nx, 1, 1)
        Av = v.view(n, 1, nx, 1)
        is_mono = (bi >= nx).view(1, 1, 1, nb)
        bspec = bi.clamp(max=nx - 1).view(nb)
        chA = ch.view(n, 1, nx, 1)
        LA = L.view(n, 1, nx, 1)
        lbC = lb.view(n, nx, 1, 1)
        lbA = lb.view(n, 1, nx, 1)
        lbB = lb[:, bspec].view(n, 1, 1, nb)
        chB = torch.where(is_mono, chA, ch[:, bspec].view(n, 1, 1, nb))
        LB = torch.where(is_mono, torch.ones_like(LA), L[:, bspec].view(n, 1, 1, nb))
        Bv = torch.where(is_mono, torch.ones_like(Av), v[:, bspec].view(n, 1, 1, nb) & (chB == chA) & (lbB == lbA))
        ok = Cv & Av & Bv & (LA + LB <= L_max) & (lbC == lbA)
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


def catalytic_triples_motif(X, p0, L_cat, lam, salt, alpha, bits, L_max, k_junc=2, chunk=4096):
    """Gokyona option cat_motif (coordinator 29 Oct 22:13; default off): motif-based catalysis within one compartment.
    Same reactions (ligation A + B -> AB, B a species or an activated monomer), same compartment rule and the same
    p_cat(L_C) = p0 (1 - exp(-(L_C - L_cat + 1)/lam)) as catalytic_triples_lab, but the hash that decides whether C
    catalyses A + B is taken over
      the catalyst's motif = its first L_cat letters, pair-symmetric (the motif or its complement, whichever is
        smaller), so a copy (the complementary strand) keeps the catalyst's identity, and
      the junction = the last k letters of A and the first k letters of B (fewer if B is shorter; a monomer = 1
        letter), so a catalyst acts on every pair of substrates with matching ends,
    with chemistry kept in each part. Inputs as catalytic_triples_lab; X["seq"] holds the full words."""
    Ra, nx = X["L"].shape
    dev = X["L"].device
    nb = nx + alpha
    out = []
    bi = torch.arange(nb, device=dev).view(1, 1, nb)
    live = ((X["valid"] & (X["L"] >= L_cat.view(-1, 1))).any(1)).nonzero().squeeze(1)
    for r0 in range(0, live.numel(), chunk):
        sl = live[r0:r0 + chunk]
        v, L, ch, lb = X["valid"][sl], X["L"][sl], X["chem"][sl], X["lab"][sl]
        sw = X["seq"][sl]                                   # [n, nx, W]
        n = v.shape[0]
        lc = L_cat[sl].view(n, 1)
        mot = sw[..., 0] & wd.mask64(lc * bits)
        motc = torch.bitwise_not(mot) & wd.mask64(lc * bits)
        motif = torch.minimum(mot, motc)                    # pair-symmetric catalyst motif [n, nx]
        kA = torch.clamp(L, max=k_junc)
        a_end = wd.shr(sw, (L - kA).clamp(min=0) * bits)[..., 0] & wd.mask64(kA * bits)
        b_start = sw[..., 0] & wd.mask64(kA * bits)         # first k letters (k_B = min(k, L_B))
        Cv = (v & (L >= lc)).view(n, nx, 1, 1)
        Av = v.view(n, 1, nx, 1)
        is_mono = (bi >= nx).view(1, 1, 1, nb)
        bspec = bi.clamp(max=nx - 1).view(nb)
        chA = ch.view(n, 1, nx, 1)
        LA = L.view(n, 1, nx, 1)
        lbC = lb.view(n, nx, 1, 1)
        lbA = lb.view(n, 1, nx, 1)
        lbB = lb[:, bspec].view(n, 1, 1, nb)
        chB = torch.where(is_mono, chA, ch[:, bspec].view(n, 1, 1, nb))
        LB = torch.where(is_mono, torch.ones_like(LA), L[:, bspec].view(n, 1, 1, nb))
        Bv = torch.where(is_mono, torch.ones_like(Av), v[:, bspec].view(n, 1, 1, nb) & (chB == chA) & (lbB == lbA))
        ok = Cv & Av & Bv & (LA + LB <= L_max) & (lbC == lbA)
        if not bool(ok.any()):
            continue
        r_, c_, a_, b_ = ok.nonzero(as_tuple=True)
        mono = b_ >= nx
        bs = b_.clamp(max=nx - 1)
        kB = torch.where(mono, torch.ones_like(b_), kA[r_, bs])
        bpart = torch.where(mono, b_ - nx, b_start[r_, bs])
        C = (ch[r_, c_], lc.view(-1)[r_], motif[r_, c_])
        A = (ch[r_, a_], kA[r_, a_], a_end[r_, a_])
        Bt = (ch[r_, a_], kB, bpart)
        rr = sl[r_]
        h = raf.hash_u01([C[0], C[1], C[2], A[0], A[1], A[2], Bt[0], Bt[1], Bt[2]], salt)
        hit = h < raf.p_cat(L[r_, c_], p0[rr], L_cat[rr], lam[rr])
        out.append(torch.stack([rr[hit], c_[hit], a_[hit], b_[hit]], 1))
    if not out:
        return torch.zeros(0, 4, dtype=torch.long, device=dev)
    return torch.cat(out)
