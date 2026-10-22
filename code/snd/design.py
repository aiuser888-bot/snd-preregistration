"""Section 3.1: block angle generation.  Counter-based Philox4x64-10; the QRNG supplies only the seed (drawn once, after T_F; 12.5).

Implements: protocol 3.1 (block angles, seed material, key, counter), 3.2 (design check), 12.5.1 step 4-5 (theta digest, bit-exact recomputation, E1).
Verified by: tests/test_snd.py::test_theta_matches_pure_python_reference

This module is the fast (numpy) implementation.  `philox_ref.py` is the same algorithm written with Python integers only; the two are
compared bit for bit by the tests and by `stream_conformance()`, because numpy does not promise that the stream of its Generator
methods stays the same across versions.  Normative definition: protocol 3.1 and the docstring of `philox_ref.py`.
"""
import functools
import hashlib
import numpy as np

from . import philox_ref

MIN_RAW_BYTES = 64   # 3.1 / 12.5.1 step 3: raw QRNG output >= 512 bits.  The entropy requirement is separate (leak.seed_entropy_ok).
_U64 = (1 << 64) - 1


def _version() -> str:
    from . import __version__
    return __version__


def seed_material(raw: bytes, device_id: int, version: str = None) -> bytes:
    """Implements: protocol 3.1 (seed material = R_d || device id || protocol version).  R_d: raw QRNG bits, >= 512.
    The byte layout is raw || b"|dev=" || decimal device id || b"|ver=" || version (ASCII)."""
    if len(raw) < MIN_RAW_BYTES:
        raise ValueError("raw QRNG output must be >= 512 bits")
    return raw + b"|dev=%d|ver=%s" % (device_id, (version if version is not None else _version()).encode())


def seed_to_key(seed_bytes: bytes) -> np.ndarray:
    """Implements: protocol 3.1 (Philox key = first 128 bits of SHA-256 of the seed material; lambda = 128 in Lemma 5).
    Key words: bytes 0..7 and 8..15 of the digest, each read as a little-endian unsigned 64-bit integer, on every platform."""
    if len(seed_bytes) < MIN_RAW_BYTES:
        raise ValueError("seed material must be >= 512 bits")
    d = hashlib.sha256(seed_bytes).digest()
    return np.array(philox_ref.key_from_digest(d), dtype=np.uint64)


def make_block_rng(key: np.ndarray, block: int, m: int = 0) -> np.random.Generator:
    """Implements: protocol 3.1 (counter = (m, b, 0, 0); b is the absolute block number), 4 (m >= 1 are the re-draws).
    A new generator per (block, m): no state is shared between blocks, re-draws, threads or processes."""
    if not (0 <= int(block) <= _U64 and 0 <= int(m) <= _U64):
        raise ValueError("block and m must be unsigned 64-bit integers")
    try:
        words = [int(k) for k in key]
    except TypeError:
        words = []
    if len(words) != 2 or not all(0 <= k <= _U64 for k in words):
        raise ValueError("key must be two unsigned 64-bit words (design.seed_to_key)")
    key = np.array(words, dtype=np.uint64)                       # always a uint64 array: numpy reads a tuple or list key as seed entropy, i.e. as a different key
    counter = np.array([m, block, 0, 0], dtype=np.uint64)
    return np.random.Generator(np.random.Philox(key=key, counter=counter))


def block_angles(rng: np.random.Generator, L: int):
    """Implements: protocol 3.1 (i) alpha ~ U[0,2pi), (ii) equispaced set, (iii) Fisher-Yates with rejection sampling.
    alpha = 2pi * (w >> 11) * 2^-53 for the first stream word w; theta_i = alpha + (2pi * pi_i) / L (binary64, this order).
    Returns (theta[L], alpha, perm)."""
    alpha = rng.uniform(0.0, 2 * np.pi)
    base = alpha + 2 * np.pi * np.arange(L) / L
    perm = rng.permutation(L)
    return base[perm], alpha, perm


def theta_rounding_bound() -> float:
    """Implements: protocol 3.1 item 7 (binary64 angle versus the ideal angle).  Upper bound, in radians, on |theta - (alpha + 2 pi pi_i / L)| for the binary64 value
    theta = fl(alpha + fl(fl(c pi_i) / L)), c = fl(2 pi), where alpha is the binary64 rotation of the block and pi is exact.  Rounding model fl(x) = x (1 + d), |d| <= 2^-53:
    (c pi_i / L)(1 + d1)(1 + d2) differs from 2 pi pi_i / L by at most (2 pi - c) + 2 pi (2 u + u^2), u = 2^-53 (the quotient is below 2 pi); the final addition
    (result below 4 pi) adds at most 4 pi u.  Exact rational arithmetic with pi to 200 bits (error 2^-199)."""
    from fractions import Fraction
    from . import highprec
    two_pi = Fraction(2 * highprec.pi_fixed(200), 1 << 200)
    u = Fraction(1, 1 << 53)
    eps_c = abs(two_pi - Fraction(2 * float(np.pi)))
    return float(eps_c + two_pi * (2 * u + u * u) + 2 * two_pi * u)


def theta_digest(key: np.ndarray, blocks, L: int, m: int = 0) -> str:
    """Implements: protocol 12.5.1 step 4 (hash of all m = 0 angles written to the audit chain).
    SHA-256 over the concatenation, for the listed absolute block numbers in ascending order, of theta_{b,0..L-1} as little-endian
    IEEE-754 binary64 (philox_ref.theta_bytes).  Hex digest."""
    blocks = list(blocks)
    if len(set(blocks)) != len(blocks):
        raise ValueError("block numbers must be distinct")
    h = hashlib.sha256()
    for b in sorted(blocks):
        h.update(block_angles(make_block_rng(key, b, m), L)[0].astype("<f8").tobytes())
    return h.hexdigest()


def design_check(theta: np.ndarray, tol: float = 1e-9) -> bool:
    """Implements: protocol 3.2 (design check): |sum e^{i theta}| <= tol*L and |sum e^{2 i theta}| <= tol*L."""
    L = theta.size
    return abs(np.exp(1j * theta).sum()) <= tol * L and abs(np.exp(2j * theta).sum()) <= tol * L


_CONFORMANCE_KEY = (0x0123456789ABCDEF, 0xFEDCBA9876543210)
_CONFORMANCE_CASES = ((0, 0, 8), (1, 1, 16), (7, 0, 100), (2 ** 20 + 3, 1999, 33), (2 ** 40, 5, 257), (3, 2, 1))


@functools.lru_cache(maxsize=None)
def stream_conformance() -> bool:
    """Implements: protocol 12.5.1 step 5 (E1: bit-exact recomputation).  Checks that the installed numpy produces exactly the
    angles of the normative algorithm (philox_ref) for a fixed set of (key, b, m, L); raises RuntimeError otherwise.
    Cheap (a few milliseconds); the result is cached per process because it depends only on the installed numpy."""
    key = np.array(_CONFORMANCE_KEY, dtype=np.uint64)
    for b, m, L in _CONFORMANCE_CASES:
        th, a, p = block_angles(make_block_rng(key, b, m), L)
        th2, a2, p2 = philox_ref.block_theta(_CONFORMANCE_KEY, b, m, L)
        if not (a == a2 and list(p) == p2 and th.astype("<f8").tobytes() == philox_ref.theta_bytes(th2)):
            raise RuntimeError("numpy %s does not reproduce the protocol 3.1 stream (block=%d, m=%d, L=%d); "
                               "use snd.philox_ref or a numpy version that passes tests/test_snd.py" % (np.__version__, b, m, L))
    return True
