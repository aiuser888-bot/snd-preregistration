"""Section 7: reduced-scale acceptance simulations (dimensionless; see protocol for scaling caveats).

Implements: protocol 7(B), 7(B′) (coverage simulations under A3 / A3′).
Verified by: tests/test_snd.py::test_coverage_iid_and_serial
"""
import numpy as np
from scipy import stats
from scipy.signal import lfilter
from .stats import center_two_way, block_V
from .fcrice import u_of_x

CRIT = stats.chi2.ppf(0.95, 2)


def wilson(k, n, z=1.96):
    """Implements: protocol 7(B) (Wilson interval used to report simulated coverage)."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def _angles(rng, K, L):
    base = 2 * np.pi * np.arange(L) / L
    a = rng.uniform(0, 2 * np.pi, (K, 1))
    return np.take_along_axis(a + base, np.argsort(rng.random((K, L)), axis=1), axis=1)  # uniform permutations


def simulate(L, K, rho, nu0, R, seed, xs, ucb, p0=0.2):
    """Implements: protocol 7(B), 7(B′).  Gaussian-copula Bernoulli(p0) H1 (A3′); f in {1.25-0.25, -0.25}; sigma_f^2 = 0.25 (p0=0.2).
    nu0 = r / sigma_0 with sigma_0 = V_iid^{-1/2} (iid normalisation).
    Returns dict(coverage, reject, kappa_emp, wilson)."""
    rng = np.random.default_rng([20261100, seed])
    EV0 = K * L * 0.25 / 2 * (1 - 1 / K)
    sig0 = 1 / np.sqrt(EV0)
    r_true = nu0 * sig0
    thr = stats.norm.ppf(1 - p0)
    cov = rej = 0
    Vs = []
    for _ in range(R):
        if abs(rho) > 0:
            e = rng.normal(size=(K, L))
            e[:, 0] /= np.sqrt(1 - rho ** 2)                    # stationary start: Var(Z_0)=1
            Z = lfilter([np.sqrt(1 - rho ** 2)], [1, -rho], e, axis=1)
        else:
            Z = rng.normal(size=(K, L))
        th = _angles(rng, K, L)
        phi = rng.uniform(0, 2 * np.pi)
        pp = p0 * (1 + nu0 * sig0 * np.cos(th - phi))
        f = 1.25 * (stats.norm.cdf(Z) > 1 - pp) - 0.25
        ft = center_two_way(f)
        U = (ft * np.exp(1j * th)).sum()
        V = block_V(ft).sum()
        xh = abs(U) / np.sqrt(V)
        rej += xh * xh > CRIT
        cov += float(u_of_x(xh, xs, ucb)) / np.sqrt(V) >= r_true
        Vs.append(V)
    return dict(coverage=cov / R, reject=rej / R, kappa_emp=float(np.mean(Vs)) / EV0,
                wilson_cov=wilson(cov, R), R=R, L=L, K=K, rho=rho, nu0=nu0)
