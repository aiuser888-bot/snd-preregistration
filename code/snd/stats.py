"""Section 2: two-way centred statistic. Only x-side data (f) is used for centring."""
import numpy as np


def center_two_way(f: np.ndarray) -> np.ndarray:
    """f: (K, L). Row-centre, subtract column profile g_j, row-centre again."""
    r = f - f.mean(axis=1, keepdims=True)
    t = r - r.mean(axis=0, keepdims=True)
    return t - t.mean(axis=1, keepdims=True)


def block_V(ft: np.ndarray) -> np.ndarray:
    L = ft.shape[1]
    return L / (2.0 * (L - 1)) * (ft ** 2).sum(axis=1)


def block_scores(ft: np.ndarray, theta: np.ndarray):
    """U = sum_b sum_i ft e^{i theta};  V = sum_b L/(2(L-1)) sum_i ft^2."""
    U = (ft * np.exp(1j * theta)).sum()
    return U, block_V(ft).sum()


def lam(U: complex, V: float) -> float:
    return abs(U) ** 2 / V
