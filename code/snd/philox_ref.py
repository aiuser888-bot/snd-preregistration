"""Protocol 3.1 (theta generation) written out with Python integers only -- no numpy, no platform dependence.

Implements: protocol 3.1 (normative algorithm), 12.5.1 step 5 (E1: independent bit-exact recomputation).
Verified by: tests/test_snd.py::test_philox_cipher_known_answers, test_theta_matches_pure_python_reference, test_known_answer_vectors_file_is_current

This is the normative algorithm of protocol 3.1 in executable form.  `design.py` calls numpy for speed; `tests/test_snd.py` checks
that numpy and this module agree bit for bit, so a change in numpy's stream (numpy does not promise the stream of Generator methods)
is detected instead of silently changing theta.  Independent re-implementations (protocol 12.5.1 step 5, E1) can be checked against
the known-answer vectors in `tests/vectors/theta_kat.json`.

Algorithm (all arithmetic on unsigned 64-bit words; Philox4x64-10 of Random123):
  * key = (k0, k1): two 64-bit words, little-endian words 0..7 and 8..15 of the first 128 bits of SHA-256(seed material).
  * stream for (block b, redraw m): counter c = (m, b, 0, 0).  The j-th block of the stream is Philox4x64-10(c + 1 + j, key), where
    c is a 256-bit little-endian integer of four 64-bit words and the +1 is part of the definition (the counter is advanced before
    the first block is used).  The stream words are v0, v1, v2, v3 of block 0, then of block 1, and so on.
  * next_u64 returns the next stream word.  next_u32 returns the low 32 bits of a fresh word and keeps the high 32 bits for the
    following call; a pending high half is consumed by the next next_u32 call and by nothing else.
  * alpha = (2 pi) * ((w >> 11) * 2^-53) for the first word w.      (IEEE-754 binary64, round to nearest even)
  * permutation: pi = (0, 1, ..., L-1); for i = L-1 down to 1: j = bounded(i) in {0..i}; swap pi[i], pi[j].
    bounded(n) for n >= 1: mask = 2^bitlength(n) - 1; repeat v = next_u32() & mask until v <= n (n < 2^32 is assumed).
  * theta_i = alpha + ((2 pi) * pi[i]) / L, evaluated as binary64 operations in this order.
"""
import math
import struct

MASK64 = (1 << 64) - 1
M0, M1 = 0xD2E7470EE14C6C93, 0xCA5A826395121157          # Philox4x64 multipliers
W0, W1 = 0x9E3779B97F4A7C15, 0xBB67AE8584CAA73B          # key increments (golden ratio, sqrt 3 - 1)
TWO_PI = 2.0 * math.pi                                   # binary64 nearest to 2 pi: 6.283185307179586


def philox4x64_10(ctr, key):
    """Implements: protocol 3.1 (Philox4x64-10 block function).  ctr: 4 words, key: 2 words -> 4 words."""
    c0, c1, c2, c3 = ctr
    k0, k1 = key
    for _ in range(10):
        p0, p1 = M0 * c0, M1 * c2
        c0, c1, c2, c3 = (p1 >> 64) ^ c1 ^ k0, p1 & MASK64, (p0 >> 64) ^ c3 ^ k1, p0 & MASK64
        k0, k1 = (k0 + W0) & MASK64, (k1 + W1) & MASK64
    return c0, c1, c2, c3


class PhiloxStream:
    """Implements: protocol 3.1 (stream of 64-bit words for a key and an initial counter; 32-bit halves; masked rejection sampling)."""
    def __init__(self, key, counter):
        self.key = tuple(int(k) & MASK64 for k in key)
        self.ctr = sum(int(c) << (64 * i) for i, c in enumerate(counter)) & ((1 << 256) - 1)   # 256-bit little-endian integer
        self.buf, self.pos = (), 4
        self.half = None

    def next_u64(self):
        if self.pos == 4:
            self.ctr = (self.ctr + 1) & ((1 << 256) - 1)
            self.buf = philox4x64_10([(self.ctr >> (64 * i)) & MASK64 for i in range(4)], self.key)
            self.pos = 0
        v = self.buf[self.pos]
        self.pos += 1
        return v

    def next_u32(self):
        if self.half is not None:
            v, self.half = self.half, None
            return v
        w = self.next_u64()
        self.half = w >> 32
        return w & 0xFFFFFFFF

    def next_double(self):
        return (self.next_u64() >> 11) * (1.0 / 9007199254740992.0)

    def bounded(self, n):
        """Uniform integer in {0..n} by masked rejection on 32-bit halves (no modulo bias)."""
        if n == 0:
            return 0
        assert n < (1 << 32)
        mask = (1 << n.bit_length()) - 1
        while True:
            v = self.next_u32() & mask
            if v <= n:
                return v


def key_from_digest(digest: bytes):
    """Implements: protocol 3.1 (Philox key = first 128 bits of SHA-256, as two little-endian 64-bit words; platform independent)."""
    return int.from_bytes(digest[0:8], "little"), int.from_bytes(digest[8:16], "little")


def block_theta(key, block: int, m: int, L: int):
    """Implements: protocol 3.1 (i)-(iii) and counter (m, b, 0, 0).  Returns (theta[L] as Python floats, alpha, perm[L]) for block `block`, re-draw `m` (m = 0: observed assignment)."""
    s = PhiloxStream(key, (m, block, 0, 0))
    alpha = TWO_PI * s.next_double()
    perm = list(range(L))
    for i in range(L - 1, 0, -1):
        j = s.bounded(i)
        perm[i], perm[j] = perm[j], perm[i]
    theta = [alpha + (TWO_PI * perm[i]) / L for i in range(L)]
    return theta, alpha, perm


def theta_bytes(theta) -> bytes:
    """Implements: protocol 3.1 / 12.5.1 step 4 (canonical serialisation of angles: each theta as little-endian IEEE-754 binary64)."""
    return struct.pack("<%dd" % len(theta), *theta)
