"""Section 10: the threshold comparisons of the decision rule, made with numerically certified values.

Implements: protocol 10 (exclusion test UCB_95 < r_min; reproduction criterion (ii)), 5.2 (UCB_95), 7(D) (reference value, precision schedule).
Verified by: tests/test_snd.py::test_decision_comparisons_certified

The p_hat comparisons of protocol 4 and 12.3 are in randtest.py.  The two comparisons here are of the same kind: a quantity computed from U is
compared with a constant.  The floating-point value decides when it is farther from the constant than its error bound; otherwise the quantity
is recomputed from the enclosure of the reference value (highprec.py) with the precision schedule of protocol 4, and the comparison is made
with a relative margin MARGIN that covers the floating-point evaluation of the table lookup and the square roots.  A comparison that
the enclosure does not decide at the last precision is reported as 'uncertain'.
"""
import math
from fractions import Fraction

import numpy as np

from . import highprec
from .fcrice import u_of_x
from .numerics import u_stat, v_stat
from .randtest import TOL_U

R_MIN = 1e-4                                       # protocol 10: exclusion threshold r_min
MARGIN = 1e-12                                     # relative margin for floating-point evaluation inside the reference branch
CHI2_LO = Fraction("5.991464547107981986870447")   # 2 ln 20 = 5.99146454710798198687044715...: the chi-square quantile of protocol 10 (ii)
CHI2_HI = Fraction("5.991464547107981986870448")
CHI2_FLOAT = 5.991464547107982


def exclusion(ft, theta_obs, xs, ucb, r_min=R_MIN, tol_u=TOL_U, p_start=highprec.P_START, p_cap=highprec.P_CAP):
    """Implements: protocol 10 (row 1: UCB_95 < r_min), 5.2.  UCB_95 = u(x_hat) / sqrt(V) < r_min  <=>  u(x_hat) < r_min sqrt(V), x_hat = |U| / sqrt(V); u is non-decreasing
    (checked on the stored table by the tests).  Returns dict(decision in {'exclude', 'not_exclude', 'uncertain'}, x_hat (floating value), ucb95 (floating value),
    reference_used (bool), precision)."""
    V = v_stat(ft)
    sv = math.sqrt(V)
    T = r_min * sv
    xh = float(abs(u_stat(ft, theta_obs))) / sv
    u = lambda x: float(u_of_x(max(x, 0.0), xs, ucb))
    out = dict(x_hat=xh, ucb95=u(xh) / sv, reference_used=False, precision=0)
    d = 2.0 * tol_u                                # x_hat = |U| / sqrt(V) is within tol_u of its reference value; 2 tol_u keeps the 7(D) factor 2
    if u(xh + d) < T * (1 - MARGIN):
        return dict(out, decision="exclude")
    if u(xh - d) >= T * (1 + MARGIN):
        return dict(out, decision="not_exclude")
    Vq, P = Fraction(V), p_start
    out["reference_used"] = True
    while P <= p_cap:
        lo2, hi2 = highprec.modulus_sq_bounds(highprec.reference_U(ft, theta_obs, P))
        x_lo, x_hi = math.sqrt(float(lo2 / Vq)) * (1 - MARGIN), math.sqrt(float(hi2 / Vq)) * (1 + MARGIN)
        out["precision"] = P
        if u(x_hi) < T * (1 - MARGIN):
            return dict(out, decision="exclude")
        if u(x_lo) >= T * (1 + MARGIN):
            return dict(out, decision="not_exclude")
        P *= 2
    return dict(out, decision="uncertain")


def reproduction_ii(ft_o, theta_o, ft_r, theta_r, tol_u=TOL_U, p_start=highprec.P_START, p_cap=highprec.P_CAP):
    """Implements: protocol 10 (reproduction criterion (ii)): |mu_rep - mu_orig|^2 / (sigma_rep^2 + sigma_orig^2) <= chi2_{2,0.95} = 2 ln 20, mu = U / V, sigma^2 = 1 / V.
    Written without division: D = |U_r V_o - U_o V_r|^2 / (V_o V_r (V_o + V_r)).  Returns dict(decision in {'pass', 'fail', 'uncertain'}, D (floating value),
    reference_used, precision).  The floating value decides when D is farther from the constant than 1e-4 (the error of D is below 10 tol_u); otherwise the
    enclosures of U_orig and U_rep give an enclosure of D."""
    Vo, Vr = v_stat(ft_o), v_stat(ft_r)
    Uo, Ur = u_stat(ft_o, theta_o), u_stat(ft_r, theta_r)
    D = float(abs(Ur * Vo - Uo * Vr) ** 2 / (Vo * Vr * (Vo + Vr)))
    out = dict(D=D, reference_used=False, precision=0)
    if D < CHI2_FLOAT - max(1e-4, 100 * tol_u):
        return dict(out, decision="pass")
    if D > CHI2_FLOAT + max(1e-4, 100 * tol_u):
        return dict(out, decision="fail")
    qo, qr = Fraction(Vo), Fraction(Vr)
    den = qo * qr * (qo + qr)
    out["reference_used"] = True
    P = p_start
    while P <= p_cap:
        re_o, im_o, e_o = highprec.reference_U(ft_o, theta_o, P)
        re_r, im_r, e_r = highprec.reference_U(ft_r, theta_r, P)
        zr, zi, ez = re_r * qo - re_o * qr, im_r * qo - im_o * qr, e_r * qo + e_o * qr      # box for Z = U_r V_o - U_o V_r
        lo2, hi2 = highprec.modulus_sq_bounds((zr, zi, ez))
        out["precision"] = P
        if hi2 / den <= CHI2_LO:
            return dict(out, decision="pass")
        if lo2 / den > CHI2_HI:
            return dict(out, decision="fail")
        P *= 2
    return dict(out, decision="uncertain")
