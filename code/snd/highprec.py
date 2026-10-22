"""Sections 4 and 7(D): reference value of the statistic in integer arithmetic, and the exact comparisons that resolve near ties.

Implements: protocol 7(D) (reference value of U and U_mon; error bound; precision schedule), 4 (second-tier resolution of near ties), 12.2 (U_mon).
Verified by: tests/test_snd.py::test_highprec_sincos_against_independent_references, test_highprec_reference_sum_and_comparison, test_near_tie_resolution

Everything here uses Python integers and fractions only: no floating-point hardware, no library cos/sin.  The inputs are binary64 numbers
(the sealed f~ of protocol 2 and the angles theta of protocol 3.1), taken as the exact binary values they are.

Reference value.  U = sum_i f~_i e^{i theta_i} as an exact real number.  `reference_U` returns an enclosure: the real and imaginary
parts lie within `err` of the returned rationals.  The enclosure comes from
  * `sincos_fixed(theta, P)`: cos and sin as integers in units of 2^-P, each within one unit (in fact within 0.5 + 1e-10 units) of the exact value.
    Method: theta is split as k*(pi/2) + r with |r| <= pi/4 using pi from Machin's formula carried with more bits than the
    magnitude of k requires, then Taylor series in integer arithmetic with GUARD extra bits, truncated when the term is below one unit; every
    floor-rounding is counted: the total working error is below 50 units of 2^-(P+GUARD), so far below one unit of 2^-P;
  * f~ taken as exact integers (mantissa * 2^e), summed exactly against the fixed-point values.
For a partition into classes (U_mon, protocol 12.2) the class sums are formed first and combined with the exact class totals F_c and sizes n_c.

Resolution.  Two enclosures are compared through the squares of their moduli (rational numbers, no square roots): |a| >= |b| is certain
when the lower bound of |a|^2 is at least the upper bound of |b|^2, |a| < |b| is certain when the upper bound of |a|^2 is below the lower
bound of |b|^2; otherwise the precision P is doubled (P_START = 128 up to P_CAP = 1024) and the comparison is repeated.
"""
import math
from fractions import Fraction

import numpy as np

P_START = 128            # first precision (bits) of the second tier, protocol 4
P_CAP = 1024             # last precision; a comparison still open there is reported as open
GUARD = 40               # extra bits carried inside sincos_fixed
_PI_CACHE = {}


def _arccot(x, bits):
    """arctan(1/x) * 2**bits for an integer x >= 2; every floor-rounding is below one unit, 3 per term at most: total error below 2**9 units."""
    x2 = x * x
    term = (1 << bits) // x
    total, k, sign = term, 1, -1
    while term:
        term //= x2
        total += sign * (term // (2 * k + 1))
        sign, k = -sign, k + 1
    return total


def pi_fixed(bits):
    """Implements: protocol 7(D) (the constant pi used by the reference value).  floor-style fixed-point pi: |result - pi * 2**bits| <= 2 units.
    Machin: pi = 16 arctan(1/5) - 4 arctan(1/239), evaluated with 16 guard bits."""
    if bits not in _PI_CACHE:
        g = bits + 16
        _PI_CACHE[bits] = (16 * _arccot(5, g) - 4 * _arccot(239, g)) >> 16
    return _PI_CACHE[bits]


def sincos_fixed(theta, P):
    """Implements: protocol 7(D) (reference value: cos and sin of a binary64 angle, exact up to 2^-P).
    Returns (c, s): integers with |c / 2**P - cos(theta)| <= 2**-P and |s / 2**P - sin(theta)| <= 2**-P.  |theta| must be below 2**24."""
    theta = float(theta)
    if not math.isfinite(theta) or abs(theta) >= 2.0 ** 24:
        raise ValueError("angle must be finite with |theta| < 2**24")
    W = P + GUARD
    X = 32                                                  # extra bits for the reduction (|k| < 2**24 costs at most 26 bits of the reduced angle)
    n, d = theta.as_integer_ratio()                         # theta = n / d exactly, d a power of two
    t = (n << (W + X)) // d                                 # floor(theta * 2**(W+X)), error < 1 unit
    half = pi_fixed(W + X) >> 1                             # pi/2, error <= 3 units
    k = (t + (half >> 1)) // half                           # nearest multiple of pi/2
    r = (t - k * half) >> X                                 # reduced angle in units of 2**-W, |r| <= (pi/4) 2**W + 2 units; error < 2 units
    sign, a = (-1, -r) if r < 0 else (1, r)
    a2 = (a * a) >> W
    cs, sn = 1 << W, a                                      # series for cos a = sum (-1)^j a^(2j)/(2j)!,  sin a = sum (-1)^j a^(2j+1)/(2j+1)!
    tc, ts, j, neg = 1 << W, a, 1, True
    while tc or ts:
        tc = ((tc * a2) >> W) // ((2 * j - 1) * (2 * j))
        ts = ((ts * a2) >> W) // ((2 * j) * (2 * j + 1))
        cs += -tc if neg else tc
        sn += -ts if neg else ts
        neg, j = not neg, j + 1
    sn *= sign
    q = k % 4                                               # rotate back: (cos, sin)(theta) from (cos, sin)(r)
    c, s = ((cs, sn), (-sn, cs), (-cs, -sn), (sn, -cs))[q]
    half_unit = 1 << (GUARD - 1)
    return (c + half_unit) >> GUARD, (s + half_unit) >> GUARD


def exact_weights(w):
    """Implements: protocol 7(D) (inputs are taken as the exact binary values).  Integers M_i and a common exponent e0 with w_i = M_i * 2**e0 exactly
    (float array or sequence; zeros give M_i = 0; non-finite values are rejected)."""
    parts = []
    for v in w:
        v = float(v)
        if not math.isfinite(v):
            raise ValueError("weights must be finite")
        if v == 0.0:
            parts.append((0, 0))
            continue
        m, e = math.frexp(v)                                # v = m * 2**e, 0.5 <= |m| < 1, m has at most 53 significant bits
        parts.append((int(m * (1 << 53)), e - 53))
    nz = [e for mm, e in parts if mm]
    e0 = min(nz) if nz else 0
    return [mm << (e - e0) if mm else 0 for mm, e in parts], e0


def reference_U(w, theta, P, classes=None):
    """Implements: protocol 7(D) (reference value of U), 12.2 (U_mon).  Enclosure (re, im, err) of U = sum_i w_i e^{i theta_i}, all exact Fractions:
    |Re U - re| <= err and |Im U - im| <= err.  w and theta: arrays of the same shape (flattened in C order), binary64.
    classes = (cls, ncls) with integer class labels of the same shape gives U_mon = sum_c F_c (n_c^-1 sum_{i in c} e^{i theta_i}) with the exact
    totals F_c = sum_{i in c} w_i; empty classes are left out.  Without classes every unit is its own class (n_c = 1, F_c = w_i)."""
    wl = np.asarray(w, float).ravel().tolist()
    tl = np.asarray(theta, float).ravel().tolist()
    if len(wl) != len(tl):
        raise ValueError("w and theta must have the same number of elements")
    ints, e0 = exact_weights(wl)
    scale = Fraction(2) ** (e0 - P)
    if classes is None:
        re = im = mass = 0
        for W_, t in zip(ints, tl):
            if W_:
                c, s = sincos_fixed(t, P)
                re += W_ * c
                im += W_ * s
                mass += abs(W_)
        return Fraction(re) * scale, Fraction(im) * scale, Fraction(mass) * scale
    cls, ncls = classes
    cl = np.asarray(cls).ravel().tolist()
    if len(cl) != len(tl):
        raise ValueError("classes must have the same number of elements as theta")
    S = [0] * ncls                                          # F_c in units of 2**e0
    R = [0] * ncls
    I = [0] * ncls
    n = [0] * ncls
    for c_, W_, t in zip(cl, ints, tl):
        c, s = sincos_fixed(t, P)
        S[c_] += W_
        R[c_] += c
        I[c_] += s
        n[c_] += 1
    re = sum((Fraction(S[c] * R[c], n[c]) for c in range(ncls) if n[c]), Fraction(0))
    im = sum((Fraction(S[c] * I[c], n[c]) for c in range(ncls) if n[c]), Fraction(0))
    mass = sum(abs(S[c]) for c in range(ncls) if n[c])
    return re * scale, im * scale, Fraction(mass) * scale


def modulus_sq_bounds(ref):
    """Implements: protocol 4 (comparison of |U| values by enclosures).  (lo, hi) Fractions with lo <= |U|^2 <= hi for every point of the enclosure box."""
    re, im, err = ref
    lo = max(abs(re) - err, Fraction(0)) ** 2 + max(abs(im) - err, Fraction(0)) ** 2
    hi = (abs(re) + err) ** 2 + (abs(im) + err) ** 2
    return lo, hi


def compare_moduli(a, b):
    """Implements: protocol 4 (|U_m| >= |U_obs| decided from enclosures).  +1 if |a| >= |b| is certain, -1 if |a| < |b| is certain, 0 if the enclosures do not decide."""
    alo, ahi = modulus_sq_bounds(a)
    blo, bhi = modulus_sq_bounds(b)
    if alo >= bhi:
        return 1
    if ahi < blo:
        return -1
    return 0


def classify(w, theta_obs, theta_of, ms, classes=None, p_start=P_START, p_cap=P_CAP):
    """Implements: protocol 4 (second-tier resolution), 7(D) (precision schedule).  Decide, for every re-draw m in `ms`, whether |U_m| >= |U_obs| in
    the reference arithmetic.  theta_of(m) returns the angles of re-draw m.  A re-draw whose angles are bit for bit those of the observed
    assignment is an exact tie and counts as >= (the definition of p_hat counts ties).  Precision starts at p_start and doubles up to p_cap.
    Returns dict(status = {m: 'ge' | 'lt' | 'open'}, precision = largest precision used, evaluations = number of reference values computed)."""
    theta_obs = np.asarray(theta_obs, float)
    status, pending, evals = {}, [], 0
    for m in ms:
        th = np.asarray(theta_of(m), float)
        if np.array_equal(th, theta_obs):
            status[m] = "ge"
        else:
            pending.append(m)
    P, used = p_start, 0
    while pending and P <= p_cap:
        ref_obs = reference_U(w, theta_obs, P, classes)
        evals += 1
        still = []
        for m in pending:
            ref_m = reference_U(w, theta_of(m), P, classes)
            evals += 1
            c = compare_moduli(ref_m, ref_obs)
            if c:
                status[m] = "ge" if c > 0 else "lt"
            else:
                still.append(m)
        pending, used = still, P
        P *= 2
    for m in pending:
        status[m] = "open"
    return dict(status=status, precision=used, evaluations=evals)
