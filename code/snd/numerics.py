"""Section 7(D): fixed-order summation and the checks that bound its error.

Implements: protocol 7(D) (numerical correctness of U and V; error bound; library qualification), 2 (statistic definition: U_b, V_b, U, V), 12.2 (U_mon in a fixed order).
Verified by: tests/test_snd.py::test_numerics_fixed_order_and_sum_checks, test_error_bound_and_library_qualification, test_near_tie_resolution

What is and is not bit-reproducible
  * theta (protocol 3.1) is: integer arithmetic plus IEEE-754 +, *, / only (philox_ref.py).
  * tree_sum below is: every step is one elementwise IEEE-754 addition in an order fixed by the algorithm, so the sum depends on
    the input values only -- not on numpy's internal blocking, SIMD width, thread count or chunking.
  * cos/sin are library calls and may differ in the last bit between platforms (they are not part of the bit-exact claim);
    protocol 7(D) bounds their effect on U (u_error_bound), checks the library against the integer reference (libm_conformance), and
    protocol 4 resolves re-draws that are closer to the observed value than that bound (highprec.py).
"""
import functools
import hashlib
import math
from fractions import Fraction
import numpy as np

from . import highprec

LIBM_ABS_TOL = 2.0 ** -50     # 7(D): the cos/sin library is assumed to be within 4 ulp, i.e. absolute error <= 2^-50 on values in [-1, 1]


def tree_sum(x, axis=-1):
    """Implements: protocol 7(D) (fixed reduction order; pairwise algorithm).  Pairwise sum in an explicit order.  Pad with +0.0 at the end to length n = 2^k, then repeat
    x <- x[..., :n/2] + x[..., n/2:] until one element is left.  Works for float64 and complex128 (componentwise)."""
    x = np.moveaxis(np.asarray(x), axis, -1)
    n = x.shape[-1]
    if n == 0:
        return np.zeros(x.shape[:-1], dtype=x.dtype)
    m = 1 << (n - 1).bit_length()
    if m != n:
        x = np.concatenate([x, np.zeros(x.shape[:-1] + (m - n,), dtype=x.dtype)], axis=-1)
    while x.shape[-1] > 1:
        h = x.shape[-1] // 2
        x = x[..., :h] + x[..., h:]
    return x[..., 0]


def kahan_sum(a):
    """Implements: protocol 7(D) (Kahan compensated summation).  Sequential in the given order (float sequence)."""
    s = c = 0.0
    for v in a:
        y = v - c
        t = s + y
        c = (t - s) - y
        s = t
    return s


def dd_sum(a):
    """Implements: protocol 7(D) (double-double accumulation as the reference).  Knuth two-sum into a hi/lo pair, sequential in the given order."""
    hi = lo = 0.0
    for v in a:
        s = hi + v
        bb = s - hi
        lo += (hi - (s - bb)) + (v - bb)
        hi, lo = s + lo, lo - ((s + lo) - s)          # fast two-sum renormalisation (|s| >= |lo|)
    return hi + lo


def exact_sum(a):
    """Implements: protocol 7(D) (second reference).  Correctly rounded sum of the exact real sum of the inputs (math.fsum)."""
    return math.fsum(a)


def u_terms(ft, theta):
    """Implements: protocol 2 (the terms f~ e^{i theta} of U).  Real and imaginary parts, elementwise (K, L)."""
    return ft * np.cos(theta), ft * np.sin(theta)


def u_stat(ft, theta):
    """Implements: protocol 2 (U), 7(D) (fixed order).  U = sum_b sum_i f~_{b,i} e^{i theta_{b,i}}.  Fixed order: tree over positions i within each block, then tree over
    blocks b = 0..K-1 (blocks in ascending index)."""
    re, im = u_terms(ft, theta)
    return complex(tree_sum(tree_sum(re, axis=1)), tree_sum(tree_sum(im, axis=1)))


def v_stat(ft):
    """Implements: protocol 2 (V), 7(D) (fixed order).  V = sum_b L/(2(L-1)) sum_i f~^2.  Same fixed order as u_stat."""
    L = ft.shape[1]
    return float(tree_sum(L / (2.0 * (L - 1)) * tree_sum(ft * ft, axis=1)))


def sum_check(ft, theta):
    """Implements: protocol 7(D).  Compute U and V by two independent algorithms in a fixed order (Kahan; pairwise tree) and compare both with a
    double-double accumulation and with the exactly rounded sum.  Returns the largest |dU|/sqrt(V) and |dV|/V over the pairs
    (tree, Kahan) x (double-double, exact), and V.  Acceptance: dU_rel <= 1e-6 and dV_rel <= 1e-9."""
    re, im = u_terms(ft, theta)
    L = ft.shape[1]
    w = L / (2.0 * (L - 1))
    vb = (w * tree_sum(ft * ft, axis=1))                        # per-block V_b, then summed in four different ways
    flat = lambda a: a.reshape(-1).tolist()
    res = {}
    for name, f in (("tree", None), ("kahan", kahan_sum), ("dd", dd_sum), ("exact", exact_sum)):
        if f is None:
            res[name] = (complex(tree_sum(tree_sum(re, axis=1)), tree_sum(tree_sum(im, axis=1))), float(tree_sum(vb)))
        else:
            res[name] = (complex(f(flat(re)), f(flat(im))), f(vb.tolist()))
    V = res["exact"][1]
    dU = max(abs(res[a][0] - res[b][0]) for a in ("tree", "kahan") for b in ("dd", "exact")) / math.sqrt(V)
    dV = max(abs(res[a][1] - res[b][1]) for a in ("tree", "kahan") for b in ("dd", "exact")) / V
    return dict(dU_rel=dU, dV_rel=dV, V=V, ok=bool(dU <= 1e-6 and dV <= 1e-9))


def class_plan(cls, ncls):
    """Implements: protocol 7(D) (fixed order of U_mon), 12.2 (classes).  For the non-empty classes in ascending class number, the flat indices (C order, i.e.
    ascending (block, position)) of their members.  The order inside a class is the order of summation."""
    c = np.asarray(cls).ravel()
    if c.size and (c.min() < 0 or c.max() >= ncls):
        raise ValueError("class labels must be in 0..ncls-1")
    order = np.argsort(c, kind="stable")
    counts = np.bincount(c, minlength=ncls)
    ends = np.cumsum(counts)
    return [order[e - n:e] for e, n in zip(ends, counts) if n]


def u_mon_stat(ft, theta, plan):
    """Implements: protocol 12.2 (U_mon = sum_c F_c Qbar_c), 7(D) (fixed order).  F_c and Q_c = sum_{i in c} e^{i theta_i} are summed by the pairwise tree over the
    members of the class in plan order; U_mon is the pairwise tree over the classes (ascending) of (F_c / n_c) Q_c, real and imaginary parts separately."""
    f, cs, sn = ft.ravel(), np.cos(theta).ravel(), np.sin(theta).ravel()
    tr, ti = [], []
    for idx in plan:
        w = float(tree_sum(f[idx])) / len(idx)
        tr.append(w * float(tree_sum(cs[idx])))
        ti.append(w * float(tree_sum(sn[idx])))
    return complex(tree_sum(np.array(tr)), tree_sum(np.array(ti)))


def error_bound_formula(K, L, sum_abs, V, plan_dims=None):
    """Implements: protocol 7(D) (worst-case bound on |dU| / sqrt(V) of the floating-point evaluation).  K blocks of L units, sum_abs = sum|f~|, V.
    Assumes the cos/sin library error is at most LIBM_ABS_TOL absolute (checked by libm_conformance) and a binary summation tree of depth d adds at most
    d 2^-53 sum|terms| per component (first order; the higher-order terms are covered by the factor g = 1 / (1 - (d + 4) 2^-53), which is 1 + 4e-15 at d = 35).
    U (plan_dims None): g sqrt(2) (LIBM_ABS_TOL + 2^-53 + d 2^-53) sum|f~| / sqrt(V), d = ceil(log2 L) + ceil(log2 K).
    U_mon (plan_dims = (largest class size, number of classes)): g (sqrt(2) (LIBM_ABS_TOL + d_n 2^-53) + d_n 2^-53 + 2 2^-53 + sqrt(2) d_c 2^-53), times
    sum|f~| / sqrt(V), d_n = ceil(log2 largest class), d_c = ceil(log2 number of classes) (class totals, division and product, and the two trees), d = d_n + d_c in g."""
    u = 2.0 ** -53
    if plan_dims is None:
        d = math.ceil(math.log2(L)) + (math.ceil(math.log2(K)) if K > 1 else 0)
        rel = math.sqrt(2) * (LIBM_ABS_TOL + u + d * u) / (1.0 - (d + 4) * u)
    else:
        nmax, ncl = plan_dims
        dn = math.ceil(math.log2(nmax)) if nmax > 1 else 0
        dc = math.ceil(math.log2(ncl)) if ncl > 1 else 0
        rel = (math.sqrt(2) * (LIBM_ABS_TOL + dn * u) + dn * u + 2 * u + math.sqrt(2) * dc * u) / (1.0 - (dn + dc + 4) * u)
    return float(rel * sum_abs / math.sqrt(V))


def u_error_bound(ft, plan=None):
    """Implements: protocol 7(D) (the bound of error_bound_formula for the data at hand: sum|f~|, V and the class structure taken from ft and `plan`)."""
    K, L = ft.shape
    dims = None if plan is None else (max(len(i) for i in plan), len(plan))
    return error_bound_formula(K, L, float(np.abs(ft).sum()), v_stat(ft), dims)


def _libm_sample():
    """Fixed sample of angles: the neighbourhood of every multiple of pi/2 up to 4 pi, tiny values, and 3000 pseudo-random angles in [0, 4 pi) from SHA-256 of a counter."""
    pts = [0.0, 5e-324, 1e-300, 1e-8, 1e-3]
    for k in range(0, 9):
        a = k * math.pi / 2
        for _ in range(3):
            pts += [a, float(np.nextafter(a, 0.0)), float(np.nextafter(a, 100.0))]
            a = float(np.nextafter(a, 100.0))
    pts += [4 * math.pi - 1e-9, float(np.nextafter(4 * math.pi, 0.0))]
    for i in range(3000):
        pts.append(4 * math.pi * int.from_bytes(hashlib.sha256(b"libm-sample-%d" % i).digest()[:8], "big") / 2.0 ** 64)
    return np.array(pts, dtype=float)


@functools.lru_cache(maxsize=None)
def libm_conformance():
    """Implements: protocol 7(D) (qualification of the cos/sin library of the analysis environment).  Compares numpy's cos and sin with the integer reference
    (highprec.sincos_fixed, 80 bits) on a fixed sample (_libm_sample) and requires absolute error <= LIBM_ABS_TOL (2^-50); raises RuntimeError otherwise.
    A sample check, not a proof for every angle.  Cached per process."""
    ang = _libm_sample()
    c, s = np.cos(ang), np.sin(ang)
    worst = 0.0
    for a, ci, si in zip(ang.tolist(), c.tolist(), s.tolist()):
        rc, rs = highprec.sincos_fixed(a, 80)
        worst = max(worst, float(abs(Fraction(ci) - Fraction(rc, 1 << 80))), float(abs(Fraction(si) - Fraction(rs, 1 << 80))))
    if worst > LIBM_ABS_TOL:
        raise RuntimeError("numpy %s: cos/sin differ from the reference by %.3g > 2^-50 on the qualification sample (protocol 7(D))" % (np.__version__, worst))
    return worst
