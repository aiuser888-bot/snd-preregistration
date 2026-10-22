"""Section 2: two-way centred statistic.  Only x-side data (f) is used for centring.

Implements: protocol 2 (two-way centring, U_b, V_b, U, V, Lambda; sealed f~ and its digest), 12.5.1 step 2 (H_f~ is part of H_F).
Verified by: tests/test_snd.py::test_lemma1_variance, test_ucb95_and_scores, test_ftilde_digest_and_seal
"""
import hashlib
import numpy as np
from .numerics import u_stat, v_stat


def center_two_way(f: np.ndarray) -> np.ndarray:
    """Implements: protocol 2 (bidirectional centring: f_{b,j} - mean_b - g_j, then row mean removed again).  f: (K, L)."""
    r = f - f.mean(axis=1, keepdims=True)
    t = r - r.mean(axis=0, keepdims=True)
    return t - t.mean(axis=1, keepdims=True)


def block_V(ft: np.ndarray) -> np.ndarray:
    """Implements: protocol 2 (V_b = L / (2 (L - 1)) * sum_i f~^2, per block).  Used by the simulations; the analysis value of V is numerics.v_stat."""
    L = ft.shape[1]
    return L / (2.0 * (L - 1)) * (ft ** 2).sum(axis=1)


def block_scores(ft: np.ndarray, theta: np.ndarray):
    """Implements: protocol 2 (U = sum_b sum_i f~ e^{i theta}, V = sum_b V_b), summed in the fixed order of protocol 7(D) (numerics.u_stat, v_stat)."""
    return u_stat(ft, theta), v_stat(ft)


def lam(U: complex, V: float) -> float:
    """Implements: protocol 2 (Lambda = |U|^2 / V)."""
    return abs(U) ** 2 / V


def ftilde_digest(ft: np.ndarray, blocks=None) -> str:
    """Implements: protocol 2 (sealed f~: canonical serialisation and H_f~), 12.5.1 step 2.  SHA-256 over the concatenation, for the listed absolute block
    numbers in ascending order, of f~_{b,0..L-1} as little-endian IEEE-754 binary64 (the serialisation of protocol 3.1 point 8).  `blocks` gives the absolute
    block number of each row of `ft` (default 0..K-1).  Values must be finite.  Hex digest.  After T_F every quantity of the analysis is a function of this
    array and of theta; the array is not recomputed from x."""
    ft = np.asarray(ft)
    if ft.ndim != 2 or not np.isfinite(ft).all():
        raise ValueError("f~ must be a finite two-dimensional array")
    blocks = list(range(ft.shape[0])) if blocks is None else [int(b) for b in blocks]
    if len(blocks) != ft.shape[0] or len(set(blocks)) != len(blocks):
        raise ValueError("one distinct block number per row is required")
    h = hashlib.sha256()
    for row in sorted(range(len(blocks)), key=lambda r: blocks[r]):
        h.update(np.ascontiguousarray(ft[row], dtype="<f8").tobytes())
    return h.hexdigest()
