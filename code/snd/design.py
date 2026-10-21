"""Section 3.1: block angle generation. Counter-based Philox; QRNG supplies only the seed (drawn once, after T_F; 12.5).

numpy's Generator.permutation uses Fisher-Yates with bounded integers drawn by
rejection (Lemire / masked), i.e. no modulo bias, as required by protocol 3.1.
"""
import hashlib
import numpy as np


MIN_RAW_BYTES = 64   # V2.0.10 (3.1): raw QRNG output >= 512 bits. The entropy requirement is separate (leak.seed_entropy_ok).


def seed_material(raw: bytes, device_id: int, version: str = "2.0.12") -> bytes:
    """Protocol 3.1: seed material = R_d || device id || protocol version (R_d: raw QRNG bits, >= 512)."""
    if len(raw) < MIN_RAW_BYTES:
        raise ValueError("raw QRNG output must be >= 512 bits")
    return raw + b"|dev=%d|ver=%s" % (device_id, version.encode())


def seed_to_key(seed_bytes: bytes) -> np.ndarray:
    """Seed material (>= 512 bits) -> 2 x uint64 Philox key = first 128 bits of SHA-256 (lambda = 128 in Lemma 5)."""
    if len(seed_bytes) < MIN_RAW_BYTES:
        raise ValueError("seed material must be >= 512 bits")
    d = hashlib.sha256(seed_bytes).digest()
    return np.frombuffer(d[:16], dtype=np.uint64).copy()


def make_block_rng(key: np.ndarray, block: int, m: int = 0) -> np.random.Generator:
    """Generator for (block, m): m = 0 is the observed assignment, m >= 1 are re-draws."""
    counter = np.array([m, block, 0, 0], dtype=np.uint64)
    return np.random.Generator(np.random.Philox(key=key, counter=counter))


def block_angles(rng: np.random.Generator, L: int):
    """alpha ~ U[0,2pi); equispaced set alpha + 2 pi j / L; uniform random permutation.
    Returns (theta[L], alpha, perm)."""
    alpha = rng.uniform(0.0, 2 * np.pi)
    base = alpha + 2 * np.pi * np.arange(L) / L
    perm = rng.permutation(L)
    return base[perm], alpha, perm


def design_check(theta: np.ndarray, tol: float = 1e-9) -> bool:
    """Section 3.2: |sum e^{i theta}| <= tol*L and |sum e^{2 i theta}| <= tol*L."""
    L = theta.size
    return abs(np.exp(1j * theta).sum()) <= tol * L and abs(np.exp(2j * theta).sum()) <= tol * L
