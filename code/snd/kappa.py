"""Lemma 2: kappa and the coverage prediction (corrected form, v15)."""
import numpy as np
from scipy import stats
from scipy.stats import multivariate_normal as mvn


def rho_x_bernoulli_copula(r, p=0.2):
    """Correlation of two Bernoulli(p) variables from a Gaussian copula with latent correlation r."""
    if abs(r) < 1e-12:
        return 0.0
    z = stats.norm.ppf(1 - p)
    both = 1 - 2 * stats.norm.cdf(z) + float(np.squeeze(mvn(mean=[0, 0], cov=[[1, r], [r, 1]]).cdf([z, z])))
    return float((both - p * p) / (p * (1 - p)))


def kappa_exact(L, rho_latent, p=0.2):
    """kappa = (L-1-2 S_L)/(L-1), S_L = sum (1-k/L) rho_k, AR(1) latent => latent corr rho^k."""
    ks = np.arange(1, L)
    rk = np.array([rho_x_bernoulli_copula(np.sign(rho_latent) ** k * abs(rho_latent) ** k, p) for k in ks])
    S = float(np.sum((1 - ks / L) * rk))
    return (L - 1 - 2 * S) / (L - 1), S


def coverage_prediction(nu, kappa, xs, ucb):
    """Lemma 2(iv): nu = r/sigma_hat (actual V). Coverage ~ ncx2.sf(x_min(nu)^2; 2, (nu/kappa)^2)."""
    xmin = xs[np.argmax(ucb >= nu)]
    return float(stats.ncx2.sf(xmin ** 2, 2, (nu / kappa) ** 2))
