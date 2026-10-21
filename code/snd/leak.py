"""Section 12 (V2.0.08; seed channel / Lemma 5 added in V2.0.10): leakage sensitivity.  Angle-mechanism leakage -> bias beta -> type-I curve alpha(beta),
worst-phase coverage under bias, and the environment-projected leakage statistic U_mon (exact randomization test)."""
import numpy as np
from scipy import stats

CHI2_2_95 = float(stats.chi2.ppf(0.95, 2))  # 5.991464547...


def type1_of_beta(beta, level=0.05):
    """Lemma 3: Lambda ~ ncx2(2, beta^2) under bias beta=|E[U|G]|/sqrt(V); type-I = P(Lambda > chi2_{2,1-level})."""
    return stats.ncx2.sf(stats.chi2.ppf(1 - level, 2), 2, np.asarray(beta, float) ** 2)


def beta_of_B(ft, q):
    """Exact conditional bias: E[U|G] = sum f~ q (f~ is x-measurable). Returns (beta=|B|/sqrt(V), omega in [0,1], q_rms)."""
    K, L = ft.shape
    B = (ft * q).sum()
    V = L / (2.0 * (L - 1)) * (ft ** 2).sum()
    qrms = np.sqrt((np.abs(q) ** 2).mean())
    dmax = beta_max(K * L, qrms, L)
    d = abs(B) / np.sqrt(V)
    return d, (d / dmax if dmax > 0 else 0.0), qrms


def beta_max(N, q_rms, L):
    """Cauchy-Schwarz: |E[U|G]|/sqrt(V) <= sqrt(2(L-1)N/L) * q_rms  (equality iff q is proportional to f~)."""
    return np.sqrt(2.0 * (L - 1) * N / L) * q_rms


def q_rms_for_beta(beta, N, L, omega=1.0):
    """Inverse of beta = omega * beta_max: the rms angle tilt that produces normalized bias beta."""
    return beta / (omega * np.sqrt(2.0 * (L - 1) * N / L))


def coverage_worst_phase(nu, beta, xs, ucb):
    """Coverage of UCB>=r under bias: xhat ~ Rice(|nu e^{i phi} + beta e^{i psi}|); worst phase gives |nu-beta|."""
    xmin = xs[np.argmax(ucb >= nu)]
    return float(stats.ncx2.sf(xmin ** 2, 2, (nu - beta) ** 2))


def coverage_min_over_nu(beta, xs, ucb, nus=None):
    nus = np.arange(0.0, 15.0001, 0.005) if nus is None else nus   # table valid for xhat<=18 (nu<=15)
    c = np.array([coverage_worst_phase(n, beta, xs, ucb) for n in nus])
    i = int(c.argmin())
    return float(c[i]), float(nus[i])


def u_mon(ft, theta, cls, ncls):
    """Environment-projected statistic: sum_i f~_i * qhat_{c(i)},  qhat_c = class mean of e^{i theta}.
    = sum_c F_c Q_c / n_c with F_c = sum_{i in c} f~_i (x/e-measurable constants), Q_c = sum_{i in c} e^{i theta_i}."""
    c = cls.ravel()
    n = np.bincount(c, minlength=ncls).astype(float)
    F = np.bincount(c, weights=ft.ravel(), minlength=ncls)
    z = np.exp(1j * theta).ravel()
    Q = np.bincount(c, weights=z.real, minlength=ncls) + 1j * np.bincount(c, weights=z.imag, minlength=ncls)
    ok = n > 0
    return (F[ok] * Q[ok] / n[ok]).sum()


# ---- Section 12.5 / Lemma 5 (V2.0.10): the seed is the only random input of theta ------------------------------------
def chi2_divergence(p_g, p0=None):
    """D = E_0[(l-1)^2] = sum_s p0(s) (p_g(s)/p0(s) - 1)^2 on a finite seed space (p0 uniform if None)."""
    p_g = np.asarray(p_g, float)
    p0 = np.full(p_g.shape, 1.0 / p_g.size) if p0 is None else np.asarray(p0, float)
    return float((p0 * (p_g / p0 - 1.0) ** 2).sum())


def beta_ub_from_D(D):
    """Lemma 5(b): beta = |E_G U|/sqrt(V) <= sqrt(D)  (equality iff l-1 is proportional to Re(conj(w) U), w = B/|B|)."""
    return float(np.sqrt(D))


def size_ub_from_D(D, level=0.05):
    """Lemma 5(c): P_G(reject) <= level + TV <= level + sqrt(D)/2   (the test is exact at `level` under the design law)."""
    return float(level + 0.5 * np.sqrt(D))


def devices_D(D1, n_dev):
    """Independent seeds given G: chi2 divergence of the product law, (1+D1)^n - 1."""
    return float(np.expm1(n_dev * np.log1p(D1)))


def hash_D_bound(h2_bits, lam=128):
    """Lemma 5(e): E_h D_out = (2^lam-1) CP(R|G) <= 2^(lam-H2(R|G)) for h drawn from a 2-universal family (SHA-256 is *modelled* so)."""
    return float(2.0 ** (lam - h2_bits))


def seed_entropy_ok(raw_bits, h_min_per_bit, lam=128, margin=32):
    """Registered requirement (12.5.1): conservative min-entropy estimate x length >= lam + margin (min-entropy <= collision entropy)."""
    return raw_bits * h_min_per_bit >= lam + margin


def best_of_k_size(k, level=0.05):
    """Lemma 5(d): using the best of k independent draws (given G) raises the size to 1-(1-level)^k."""
    return float(1.0 - (1.0 - level) ** k)


def retry_size_ub(f, level=0.05):
    """Lemma 5(d'): re-drawing the seed after a theta-only check fails (first passing draw) gives law P_0(.|pass); size <= level/(1-f)."""
    return float(level / (1.0 - f))
