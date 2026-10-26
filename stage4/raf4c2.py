"""Stage 4c-2 building blocks: fixed-hash catalysis assignment and a sparse maxRAF detector.

Catalysis (spec docs/stage4c2_spec.md, section 1): sequence C catalyses the ligation A + B -> AB
if hash(C, A, B, salt) < p_cat(L_C), with p_cat(L) = p0 * (1 - exp(-(L - L_cat + 1) / lam)) for
L >= L_cat, else 0. A molecule is (chem, L, seq): chemistry flag, length, bit-packed sequence
(monomer i at bits [i*BITS, (i+1)*BITS), as in chem4c1). The hash is splitmix64 chained over
the nine integers, so the assignment is random but reproducible: the same molecules have the
same chemistry in every site and run. `salt` is a global constant; the shuffled control
(validation test 3) changes it every step.

maxRAF (Hordijk & Steel 2004): repeat { closure of the food set under the active reactions;
drop reactions with a reactant outside the closure or no catalyst inside it } until nothing
changes. The result is the maximal RAF; it is non-empty iff a RAF exists. Reactions are
ligations A + B -> AB; `reversible=True` also lets AB -> A + B extend the closure (the
convention of the stage-1e polymer model, used for the threshold test).
"""
import math
import torch

MASK64 = (1 << 64) - 1
GOLDEN = 0x9E3779B97F4A7C15
M1 = 0xBF58476D1CE4E5B9
M2 = 0x94D049BB133111EB


def _s64(x):
    """Python int (mod 2^64) -> signed int64 value."""
    x &= MASK64
    return x - (1 << 64) if x >= (1 << 63) else x


def _srl(z, k):
    """Logical right shift of int64 tensors."""
    return torch.bitwise_right_shift(z, k) & ((1 << (64 - k)) - 1)


def _mix(z):
    z = z + _s64(GOLDEN)
    z = (z ^ _srl(z, 30)) * _s64(M1)
    z = (z ^ _srl(z, 27)) * _s64(M2)
    return z ^ _srl(z, 31)


def hash_u01(parts, salt=0):
    """Uniform [0,1) float64 from a list of broadcastable int64 tensors (splitmix64 chain)."""
    z = None
    for v in parts:
        v = torch.as_tensor(v, dtype=torch.int64)
        z = _mix(v + _s64(salt)) if z is None else _mix(z ^ v)
    return _srl(z, 11).double() / float(1 << 53)


def p_cat(L, p0, L_cat, lam):
    L = torch.as_tensor(L, dtype=torch.float64)
    p = p0 * (1 - torch.exp(-(L - L_cat + 1) / lam))
    return torch.where(L >= L_cat, p, torch.zeros_like(p))


def catalyses(C, A, B, p0, L_cat, lam, salt=0):
    """C, A, B: tuples (chem, L, seq) of broadcastable int64 tensors. Returns bool tensor."""
    h = hash_u01([C[0], C[1], C[2], A[0], A[1], A[2], B[0], B[1], B[2]], salt)
    return h < p_cat(C[1], p0, L_cat, lam)


def max_raf(n_species, food, ra, rb, rab, cat_r, cat_c, reversible=False):
    """Sparse maxRAF. Species ids 0..n-1; food bool [n]; reactions ra + rb -> rab (long [R]);
    catalysis pairs (reaction cat_r, catalyst species cat_c). Returns bool [R]: reactions in
    the maxRAF (all False if there is no RAF)."""
    R = ra.numel()
    dev = ra.device
    active = torch.zeros(R, dtype=torch.bool, device=dev)
    if R == 0:
        return active
    active[cat_r] = True
    while True:
        avail = food.clone()
        while True:
            fire = active & avail[ra] & avail[rb]
            new = avail.clone()
            new[rab[fire]] = True
            if reversible:
                back = active & avail[rab]
                new[ra[back]] = True
                new[rb[back]] = True
            if torch.equal(new, avail):
                break
            avail = new
        has_cat = torch.zeros(R, dtype=torch.bool, device=dev)
        ok = avail[cat_c]
        has_cat[cat_r[ok]] = True
        reach = avail[ra] & avail[rb]
        if reversible:
            reach = reach | avail[rab]
        new_active = active & has_cat & reach
        if torch.equal(new_active, active):
            return active
        active = new_active


def raf_components(ra, rb, rab, in_raf):
    """Number of connected components (shared species) among the maxRAF reactions."""
    idx = in_raf.nonzero().squeeze(1).tolist()
    if not idx:
        return 0
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in idx:
        s = [int(ra[i]), int(rb[i]), int(rab[i])]
        for x in s[1:]:
            parent[find(x)] = find(s[0])
    return len({find(int(rab[i])) for i in idx})


def full_polymer_network(n, A=2, chem=0):
    """All strings of length 1..n over A letters (one chemistry): species tensors and the
    ligation reactions prefix + suffix -> string (as in the stage-1e model)."""
    bits = max(1, (A - 1).bit_length())
    Ls, seqs = [], []
    for L in range(1, n + 1):
        for v in range(A ** L):
            s, x = 0, v
            for i in range(L):
                s |= (x % A) << (i * bits)
                x //= A
            Ls.append(L)
            seqs.append(s)
    L = torch.tensor(Ls, dtype=torch.int64)
    seq = torch.tensor(seqs, dtype=torch.int64)
    index = {(int(l), int(s)): i for i, (l, s) in enumerate(zip(Ls, seqs))}
    ra, rb, rab = [], [], []
    for i, (l, s) in enumerate(zip(Ls, seqs)):
        for j in range(1, l):
            left = s & ((1 << (j * bits)) - 1)
            right = s >> (j * bits)
            ra.append(index[(j, left)])
            rb.append(index[(l - j, right)])
            rab.append(i)
    t = lambda z: torch.tensor(z, dtype=torch.int64)
    return torch.full_like(L, chem), L, seq, t(ra), t(rb), t(rab)
