"""Section 4: randomization test. V is fixed under re-draws, so Lambda_m >= Lambda_obs <=> |U_m| >= |U_obs|."""
import numpy as np
from .design import make_block_rng, block_angles


def draw_theta(key, K, L, m):
    return np.stack([block_angles(make_block_rng(key, b, m), L)[0] for b in range(K)])


def randomization_test(ft, key, M=1999, theta_obs=None):
    """ft: (K,L) centred f. theta_obs defaults to the m=0 draw. Returns dict with p_hat, tail count, U_obs, V."""
    K, L = ft.shape
    if theta_obs is None:
        theta_obs = draw_theta(key, K, L, 0)
    U_obs = (ft * np.exp(1j * theta_obs)).sum()
    V = (L / (2.0 * (L - 1)) * (ft ** 2).sum(axis=1)).sum()
    cnt = 0
    for m in range(1, M + 1):
        U_m = (ft * np.exp(1j * draw_theta(key, K, L, m))).sum()
        cnt += abs(U_m) >= abs(U_obs)
    return dict(p_hat=(1 + cnt) / (M + 1), tail=int(cnt), M=M, U_obs=U_obs, V=V,
                Lambda_obs=abs(U_obs) ** 2 / V, reject=(1 + cnt) / (M + 1) <= 0.05)
