"""Multi-word bit sequences for Gokyona Environment. A chain of length L (alphabet 2, one bit per letter) is stored
in W int64 words, letter i at bit i % 64 of word i // 64 (word 0 = letters 0-63). Tensors of sequences have a
trailing dimension W. With W = 1 every function reproduces the single-int64 operation of chem4c4 exactly
(so the multi-word engine with W = 1 is bit-identical to chem4c4)."""
import torch

W = 1                       # words per chain (set by gk_chem.configure)
_ALL = -1                   # int64 with all 64 bits set


def set_words(w):
    global W
    W = int(w)


def lmax():
    return 64 * W


def mask64(nbits):
    """Single-word mask of the low nbits bits (nbits clamped to 0..64)."""
    one = torch.ones_like(nbits)
    nb = nbits.clamp(min=0)
    return torch.where(nb >= 64, torch.full_like(nb, _ALL), torch.bitwise_left_shift(one, nb.clamp(max=63)) - 1)


def mask(nbits):
    """[...] -> [..., W]: words with the low nbits letters set."""
    w = torch.arange(W, device=nbits.device)
    return mask64(nbits.unsqueeze(-1) - 64 * w)


def srl64(x, r):
    """Logical right shift of int64 x by r (0..63), elementwise."""
    r = r.clamp(min=0, max=63)
    return torch.bitwise_right_shift(x, r) & mask64(64 - r)


def shr(seq, k):
    """Logical right shift of multi-word sequences by k letters (k [...] >= 0). Bits shifted in are 0."""
    if W == 1:
        return torch.where(k[..., None] >= 64, torch.zeros_like(seq), srl64(seq, k[..., None]))
    q, r = torch.div(k, 64, rounding_mode="floor"), k % 64
    w = torch.arange(W, device=seq.device)
    i0 = q[..., None] + w                                  # source word for the low part
    i1 = i0 + 1
    g = lambda i: torch.where(i < W, torch.gather(seq, -1, i.clamp(max=W - 1)), torch.zeros_like(seq))
    lo = srl64(g(i0), r[..., None].expand_as(i0))
    hi = torch.where(r[..., None] > 0, torch.bitwise_left_shift(g(i1), (64 - r[..., None]).clamp(max=63)), torch.zeros_like(seq))
    return lo | hi


def shl(seq, k):
    """Left shift by k letters (bits beyond W words are dropped)."""
    if W == 1:
        return torch.where(k[..., None] >= 64, torch.zeros_like(seq), torch.bitwise_left_shift(seq, k[..., None].clamp(max=63)))
    q, r = torch.div(k, 64, rounding_mode="floor"), k % 64
    w = torch.arange(W, device=seq.device)
    i0 = w - q[..., None]                                  # source word for the high part
    i1 = i0 - 1
    g = lambda i: torch.where((i >= 0) & (i < W), torch.gather(seq, -1, i.clamp(min=0, max=W - 1)), torch.zeros_like(seq))
    hi = torch.bitwise_left_shift(g(i0), r[..., None].expand_as(i0))
    lo = torch.where(r[..., None] > 0, srl64(g(i1), (64 - r[..., None]).clamp(max=63)), torch.zeros_like(seq))
    return hi | lo


def from_low(x):
    """int64 [...] (a value of <= 64 letters) -> [..., W] with x in word 0."""
    out = torch.zeros(x.shape + (W,), dtype=torch.int64, device=x.device)
    out[..., 0] = x
    return out


def pack(bits):
    """bool/0-1 [..., n] with n <= 64 W letters -> [..., W] words."""
    n = bits.shape[-1]
    pad = W * 64 - n
    b = bits.long()
    if pad > 0:
        b = torch.cat([b, torch.zeros(b.shape[:-1] + (pad,), dtype=torch.long, device=b.device)], -1)
    b = b.view(b.shape[:-1] + (W, 64))
    wbit = torch.bitwise_left_shift(torch.ones(64, dtype=torch.long, device=b.device), torch.arange(64, device=b.device))
    return (b * wbit).sum(-1)


def unpack(seq, n):
    """[..., W] -> [..., n] 0/1 letters (n <= 64 W)."""
    j = torch.arange(n, device=seq.device)
    word = torch.gather(seq, -1, (j // 64).expand(seq.shape[:-1] + (n,)))
    return (torch.bitwise_right_shift(word, j % 64) & 1)


def popcount(x):
    """Bit count per int64 element (vectorised)."""
    x = x - (srl64(x, torch.ones_like(x)) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + (srl64(x, torch.full_like(x, 2)) & 0x3333333333333333)
    x = (x + srl64(x, torch.full_like(x, 4))) & 0x0F0F0F0F0F0F0F0F
    return srl64(x * 0x0101010101010101, torch.full_like(x, 56))


def to_str(words, n, cap=None):
    """Host helper: words (list/tensor of W ints) -> '0101...' of the first min(n, cap) letters."""
    ws = [int(v) for v in words]
    m = n if cap is None else min(n, cap)
    return "".join(str((ws[i // 64] >> (i % 64)) & 1) for i in range(m))
