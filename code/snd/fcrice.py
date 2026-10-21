"""Section 5.1: Feldman-Cousins ordering for Rice(nu,1); one-sided UCB u(x)."""
import numpy as np
from scipy.special import ive

TABLE_CHECKPOINTS = {0.1: 1.650, 1.0: 2.565, 1.794: 3.533, 2.448: 4.245, 3.0: 4.825, 5.0: 6.870, 10.0: 11.915, 15.0: 16.93}
XMAX_VALID = 18.0


def build_table(dx=0.004, xmax=22.0, dn=0.005, numax=20.0):
    xs = np.arange(dx, xmax, dx)
    nus = np.arange(0, numax + dn / 2, dn)
    LP = np.empty((len(nus), len(xs)))
    for i, nu in enumerate(nus):
        LP[i] = np.log(xs) - (xs - nu) ** 2 / 2 + np.log(ive(0, xs * nu))
    LPm = LP.max(axis=0)
    acc = np.zeros(LP.shape, bool)
    mass = np.empty(len(nus))
    for i in range(len(nus)):
        R = LP[i] - LPm
        p = np.exp(LP[i]) * dx
        o = np.argsort(-R)
        k = np.searchsorted(np.cumsum(p[o]), 0.95)
        acc[i] = R >= R[o[k]]
        mass[i] = p[acc[i]].sum()
    ucb = np.array([nus[acc[:, j]].max() for j in range(len(xs))])
    return xs, ucb, mass


def load_table(path):
    d = np.load(path)
    return d["xs"], d["ucb"]


def u_of_x(x, xs, ucb):
    x = np.asarray(x, float)
    return np.where(x <= XMAX_VALID, np.interp(x, xs, ucb), x + 1.96)


def ucb95(U, V, xs, ucb):
    """UCB_95 = sigma_hat * u(x_hat), sigma_hat = V^{-1/2}, x_hat = |U|/sqrt(V)."""
    xh = abs(U) / np.sqrt(V)
    return float(u_of_x(xh, xs, ucb)) / np.sqrt(V)
