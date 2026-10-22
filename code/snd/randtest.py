"""Section 4: randomization test.  V is fixed under re-draws, so Lambda_m >= Lambda_obs <=> |U_m| >= |U_obs|.

Implements: protocol 4 (randomization test, p_hat, M_perm = 1999, near ties and the second-tier resolution), 3.1 (re-draws: same key, counter (b, m)), 7(D) (numerical uncertainty), 12.2 (p_hat_mon).
Verified by: tests/test_snd.py::test_randtest_size_under_H0, test_randtest_independent_of_workers_and_order, test_near_tie_resolution

Reproducibility.  Nothing here touches a global random state: the angles of block b at re-draw m come from a fresh generator that
depends only on (key, b, m) (design.make_block_rng).  A re-draw m can therefore be computed in any order, in any process, on any number
of workers; tail counts are integers and add exactly.  Bit-exact quantities: theta (protocol 3.1) and the fixed-order sums
(numerics.tree_sum).  Not bit-exact across platforms: cos/sin inside U_m.

Numerical uncertainty (protocol 4, 7(D)).  p_hat is defined by the reference values of |U| (exact real arithmetic on the sealed f~ and the binary64
angles; highprec.py).  The floating-point statistic is an approximation with error at most TOL_U*sqrt(V) (numerics.u_error_bound is the worst-case
bound for the data at hand; numerics.libm_conformance checks the cos/sin library).  A re-draw whose floating value is within
delta = 2*TOL_U*sqrt(V) of |U_obs| is a near tie; every other re-draw has the same order relative to |U_obs| as its reference value.  If the
near ties cannot change the decision (p_hat_lo and p_hat_hi on the same side of 0.05) the decision is that side; if they can
(p_hat_lo <= 0.05 < p_hat_hi) the near ties and |U_obs| are recomputed in the reference arithmetic (highprec.classify) -- always, with no choice left to the analyst.
"""
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import highprec
from .design import make_block_rng, block_angles, stream_conformance
from .numerics import u_stat, v_stat, u_mon_stat, class_plan, u_error_bound, libm_conformance

TOL_U = 1e-6        # 7(D): |dU| / sqrt(V) <= 1e-6 is the accepted numerical error of U
ALPHA = 0.05        # 4: reject H0 iff p_hat <= 0.05


def draw_theta(key, K, L, m):
    """Implements: protocol 3.1 / 4.  Angles of blocks 0..K-1 at re-draw m (m = 0: observed assignment), shape (K, L)."""
    return np.stack([block_angles(make_block_rng(key, b, m), L)[0] for b in range(K)])


def _stat(ft, theta, plan):
    return u_stat(ft, theta) if plan is None else u_mon_stat(ft, theta, plan)


def abs_U_of_draw(ft, key, m, plan=None):
    """Implements: protocol 4 (|U_m|), 7(D).  |U_m| for re-draw m (m = 0: observed); with `plan` (numerics.class_plan) the statistic is U_mon (12.2).
    Pure function of its arguments; used by the reproducibility tests."""
    return float(abs(_stat(ft, draw_theta(key, ft.shape[0], ft.shape[1], m), plan)))


def tail_counts(ft, key, ms, absU_obs, delta, plan=None):
    """Implements: protocol 4 (tail count) and 7(D) (near ties).  For re-draws m in `ms` return (c, c_lo, c_hi, ties): integer counts of
    |U_m| >= |U_obs|,  |U_m| >= |U_obs| + delta,  |U_m| >= |U_obs| - delta,  and the sorted list of the near ties
    (those m with |U_obs| - delta <= |U_m| < |U_obs| + delta).  Pure function of its arguments (no global state);
    counts over disjoint sets of m add, lists merge."""
    K, L = ft.shape
    c = c_lo = c_hi = 0
    ties = []
    for m in ms:
        a = abs(_stat(ft, draw_theta(key, K, L, m), plan))
        c += a >= absU_obs
        c_lo += a >= absU_obs + delta
        c_hi += a >= absU_obs - delta
        if absU_obs - delta <= a < absU_obs + delta:
            ties.append(int(m))
    return int(c), int(c_lo), int(c_hi), sorted(ties)


def _chunk(args):                                  # top-level so that worker processes can import it
    ft, key, ms, absU_obs, delta, plan = args
    stream_conformance()
    libm_conformance()
    return tail_counts(ft, key, ms, absU_obs, delta, plan)


def split_ms(M, parts):
    """Implements: protocol 4 (re-draws 1..M; any partition gives the same tail count).  Contiguous, near-equal chunks of 1..M."""
    edges = [1 + (M * i) // parts for i in range(parts + 1)]
    return [range(edges[i], edges[i + 1]) for i in range(parts) if edges[i + 1] > edges[i]]


def randomization_test(ft, key, M=1999, theta_obs=None, workers=1, chunks=None, tol_u=TOL_U, classes=None, resolve=True,
                       p_start=highprec.P_START, p_cap=highprec.P_CAP):
    """Implements: protocol 4, 7(D), 12.2 (p_hat_mon with `classes`).  ft: (K, L) sealed f~.  theta_obs defaults to the m = 0 draw.
    workers > 1 evaluates the re-draws in that many processes (start method 'spawn': the calling program must be a script file with an
    `if __name__ == "__main__"` guard, as for any multiprocessing use); `chunks` (list of iterables of m, a partition of 1..M)
    overrides the default split -- the result does not depend on either.
    classes = (cls, ncls) replaces U by U_mon (protocol 12.2) built from the class labels `cls` (same shape as ft); V stays the V of ft.
    resolve=False skips the second tier (used only to look at the floating-point report).  p_start / p_cap: precision schedule of the second tier.
    Returns dict: p_hat = (1 + #{m: |U_m| >= |U_obs|}) / (M + 1), tail, M, U_obs (floating value), V, Lambda_obs (None for U_mon), reject (decision == 'reject'),
    near_ties (#{m: |U_obs| - delta <= |U_m| < |U_obs| + delta}, delta = 2 tol_u sqrt(V)), near_tie_draws, p_hat_lo / p_hat_hi (ties counted as below /
    at or above), uncertain, decision in {reject, not_reject, uncertain}, error_bound (numerics.u_error_bound), p_hat_float (before the second tier) and
    resolution (None, or dict(precision, evaluations, n_ties, n_ge, n_lt, n_open) of the second tier).  After the second tier p_hat, tail, p_hat_lo, p_hat_hi
    and decision are those of the reference values; a comparison left open at p_cap counts toward the tail in p_hat (and in p_hat_hi only)."""
    stream_conformance()
    libm_conformance()
    if M < 1:
        raise ValueError("M must be at least 1")
    K, L = ft.shape
    plan = None if classes is None else class_plan(*classes)
    bound = u_error_bound(ft, plan)
    if bound > TOL_U:
        raise ValueError("7(D): worst-case numerical error %.3g sqrt(V) of U exceeds the accepted %.0e" % (bound, TOL_U))
    if theta_obs is None:
        theta_obs = draw_theta(key, K, L, 0)
    U_obs = _stat(ft, theta_obs, plan)
    V = v_stat(ft)
    absU, delta = abs(U_obs), 2.0 * tol_u * np.sqrt(V)
    if chunks is None:
        chunks = split_ms(M, max(1, workers))
    if sorted(m for c in chunks for m in c) != list(range(1, M + 1)):
        raise ValueError("chunks must partition 1..M")
    jobs = [(ft, key, list(c), absU, delta, plan) for c in chunks]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
            parts = list(ex.map(_chunk, jobs))
    else:
        parts = [tail_counts(*j) for j in jobs]
    cnt, c_lo, c_hi = (sum(p[i] for p in parts) for i in range(3))
    ties = sorted(m for p in parts for m in p[3])
    p_hat, p_lo, p_hi = ((1 + c) / (M + 1) for c in (cnt, c_lo, c_hi))
    uncertain = (p_lo <= ALPHA) != (p_hi <= ALPHA)
    out = dict(p_hat=p_hat, tail=cnt, M=M, U_obs=U_obs, V=V, Lambda_obs=None if plan is not None else absU ** 2 / V,
               near_ties=c_hi - c_lo, near_tie_draws=ties, p_hat_lo=p_lo, p_hat_hi=p_hi, error_bound=bound, p_hat_float=p_hat, resolution=None)
    resolution = None
    if uncertain and resolve:
        res = highprec.classify(ft, theta_obs, lambda m: draw_theta(key, K, L, m), ties, classes, p_start, p_cap)
        st = res["status"]
        n_ge, n_lt, n_open = (sum(1 for v in st.values() if v == s) for s in ("ge", "lt", "open"))
        cnt = c_lo + n_ge + n_open                                # open comparisons count toward the tail (ties count, p_hat is conservative)
        p_hat, p_lo, p_hi = (1 + cnt) / (M + 1), (1 + c_lo + n_ge) / (M + 1), (1 + c_lo + n_ge + n_open) / (M + 1)
        uncertain = (p_lo <= ALPHA) != (p_hi <= ALPHA)
        resolution = dict(precision=res["precision"], evaluations=res["evaluations"], n_ties=len(ties), n_ge=n_ge, n_lt=n_lt, n_open=n_open)
    decision = "uncertain" if uncertain else ("reject" if p_hat <= ALPHA else "not_reject")
    if resolution is None:                                        # floating-point report only: keep p_hat <= 0.05 as the plain comparison
        out.update(reject=p_hat <= ALPHA)
    else:
        out.update(reject=decision == "reject")
    out.update(p_hat=p_hat, tail=cnt, p_hat_lo=p_lo, p_hat_hi=p_hi, uncertain=bool(uncertain), decision=decision, resolution=resolution)
    return out
