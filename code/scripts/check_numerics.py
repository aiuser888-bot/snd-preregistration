"""Protocol 7(D): numerical correctness of the sums behind U and V, measured at reduced scale.  Writes results/numerics_check.json.

Implements: protocol 7(D).
Verified by: tests/test_snd.py::test_numerics_fixed_order_and_sum_checks, test_error_bound_and_library_qualification

For each (K, L): centred f~ from a Bernoulli(0.2) sample, theta from the 3.1 generator (m = 0); U and V are summed in a fixed order by the pairwise
tree and by Kahan and compared with a double-double accumulation and with the exactly rounded sum (numerics.sum_check).  Also reported:
the effect of the cos/sin implementation on U (numpy vs the C library, same angles), the largest cos/sin difference in ulp, and the whole floating-point
U (numpy cos/sin, fixed-order tree) against the integer reference value (highprec.reference_U, 80 bits: error below 1e-18 sqrt(V) here) together with the
worst-case bound of numerics.u_error_bound.  The last block times the integer reference and extrapolates its cost to the registered scale.
Acceptance (7(D)): |dU|/sqrt(V) <= 1e-6 and |dV|/V <= 1e-9, measured error against the reference <= the bound <= 1e-6.   Exit status 1 if any case fails.
"""
import json, math, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from snd.design import seed_to_key, seed_material
from snd import highprec
from snd.numerics import sum_check, u_stat, v_stat, exact_sum, u_error_bound, error_bound_formula, libm_conformance
from snd.randtest import draw_theta
from snd.stats import center_two_way

ROOT = os.path.join(os.path.dirname(__file__), "..", "..", "results")


def one(K, L, seed):
    rng = np.random.default_rng([20261300, seed])
    ft = center_two_way(1.25 * (rng.random((K, L)) < 0.2) - 0.25)
    th = draw_theta(seed_to_key(seed_material(bytes(range(64)), 1)), K, L, 0)
    r = sum_check(ft, th)
    V = v_stat(ft)
    cs = np.cos(th), np.sin(th)
    lib = (np.array([math.cos(a) for a in th.reshape(-1).tolist()]).reshape(th.shape),
           np.array([math.sin(a) for a in th.reshape(-1).tolist()]).reshape(th.shape))
    U_np = complex(exact_sum((ft * cs[0]).reshape(-1).tolist()), exact_sum((ft * cs[1]).reshape(-1).tolist()))
    U_lib = complex(exact_sum((ft * lib[0]).reshape(-1).tolist()), exact_sum((ft * lib[1]).reshape(-1).tolist()))
    ulp = max(float((np.abs(cs[i] - lib[i]) / np.spacing(np.abs(lib[i]))).max()) for i in (0, 1))
    re, im, err = highprec.reference_U(ft, th, 80)                                      # exact value of U for these f~ and theta, to within err
    U_fl = u_stat(ft, th)
    d_ref = math.hypot(U_fl.real - float(re), U_fl.imag - float(im)) / math.sqrt(V)
    return dict(K=K, L=L, N=K * L, dU_rel_sum=r["dU_rel"], dV_rel_sum=r["dV_rel"], ok_sum=r["ok"],
                dU_rel_trig=abs(U_np - U_lib) / math.sqrt(V), max_trig_ulp=ulp, abs_U_over_sqrtV=abs(U_np) / math.sqrt(V),
                dU_rel_vs_reference=d_ref, reference_enclosure_rel=float(err) / math.sqrt(V), error_bound=u_error_bound(ft))


if __name__ == "__main__":
    t0 = time.time()
    rows = [one(K, L, s) for s, (K, L) in enumerate(((50, 100), (100, 1000), (100, 10000), (1000, 1000)))]
    for r in rows:
        print("K=%5d L=%6d N=%8d  dU/sqrtV sums %.2e  dV/V sums %.2e  | cos/sin: dU/sqrtV %.2e, max %.1f ulp | vs reference: %.2e (enclosure %.1e, bound %.2e)  %s"
              % (r["K"], r["L"], r["N"], r["dU_rel_sum"], r["dV_rel_sum"], r["dU_rel_trig"], r["max_trig_ulp"], r["dU_rel_vs_reference"],
                 r["reference_enclosure_rel"], r["error_bound"], "OK" if r["ok_sum"] and r["dU_rel_vs_reference"] <= r["error_bound"] <= 1e-6 else "FAIL"))
    # Worst case at the registered scale (N = 1.2e10, V_target = 1.5e9).  |f~| <= 6: f~ = (f - mean_j f) - (g_j - mean_j g) with |f| <= 1 (A2),
    # so |f - mean_j f| <= 2 and |g_j - mean_j g| <= 4.  Each term f~ (cos, sin)(theta) has components with relative error at most
    # u ulp (<= u 2^-52 of the value) from the library plus 2^-53 from the product, all in the same direction: the error vector of one term has
    # length <= sqrt(2) (u 2^-52 + 2^-53) |f~|, so |dU| <= sqrt(2) (u 2^-52 + 2^-53) sum|f~| <= sqrt(2) (u 2^-52 + 2^-53) 6 N.
    # A binary summation tree of depth d adds at most d 2^-53 sum|terms| per component, i.e. sqrt(2) d 2^-53 sum|f~| in length.
    u53, N, V, fmax = 2.0 ** -53, 1.2e10, 1.5e9, 6.0
    bound = {u: math.sqrt(2) * (u * 2.0 ** -52 + u53) * fmax * N / math.sqrt(V) for u in (1, 4)}
    print("worst-case |dU|/sqrt(V) at N=1.2e10, V=1.5e9, |f~|<=6, cos/sin error of 1 / 4 ulp on every term: %.2e / %.2e (limit 1e-6)" % (bound[1], bound[4]))
    depth = math.ceil(math.log2(10000)) + math.ceil(math.log2(N / 10000))             # tree depth: positions within a block, then blocks
    tree_bound = math.sqrt(2) * depth * u53 * fmax * N / math.sqrt(V)                   # worst-case summation error of the pairwise tree
    print("worst-case |dU|/sqrt(V) of the pairwise tree itself (depth %d): %.2e (limit 1e-6)" % (depth, tree_bound))
    full = error_bound_formula(N / 10000, 10000, fmax * N, V)                           # the bound of numerics.u_error_bound with sum|f~| <= 6 N
    print("numerics.error_bound_formula at the registered scale (sum|f~| <= 6 N): %.2e (limit 1e-6)" % full)
    # cost of one reference evaluation (integer arithmetic, 128 bits): time per term, extrapolated to one evaluation of U at N = 1.2e10
    rng = np.random.default_rng(1)
    ang = rng.uniform(0, 4 * math.pi, 20000).tolist()
    def best(fn, repeats=5):                                     # the fastest of several runs: the least disturbed by other load on the machine
        times = []
        for _ in range(repeats):
            t1 = time.perf_counter(); fn(); times.append(time.perf_counter() - t1)
        return min(times) / len(ang) * 1e6
    us = best(lambda: [highprec.sincos_fixed(a, 128) for a in ang])
    us_term = best(lambda: highprec.reference_U(np.ones(20000), np.array(ang), 128))
    cost = dict(sincos_us=us, reference_U_us_per_term=us_term, extrapolated_core_days_per_evaluation_N_1p2e10=us_term * 1.2e10 / 1e6 / 86400,
                note="fastest of 5 runs on this machine, one core, pure Python; extrapolation, not measured at the registered scale")
    print("reference value: %.1f us per cos/sin pair, %.1f us per term of U (128 bits) -> %.1f core-days per evaluation at N = 1.2e10 (extrapolated)"
          % (us, us_term, cost["extrapolated_core_days_per_evaluation_N_1p2e10"]))
    lib_worst = libm_conformance()
    print("library qualification (numerics.libm_conformance): largest absolute cos/sin error on the sample %.2e (limit 2^-50 = %.2e)" % (lib_worst, 2.0 ** -50))
    ok = (all(r["ok_sum"] and r["dU_rel_trig"] <= 1e-6 and r["max_trig_ulp"] <= 4 and r["dU_rel_vs_reference"] <= r["error_bound"] <= 1e-6 for r in rows)
          and bound[4] + tree_bound <= 1e-6 and full <= 1e-6)
    json.dump(dict(rows=rows, worst_case_bound=dict(N=N, V=V, f_abs_max=fmax, ulp1=bound[1], ulp4=bound[4], tree_depth=depth, tree_sum=tree_bound, formula_total=full),
                   libm_sample_worst_abs_error=lib_worst, timing=cost, seconds=time.time() - t0, ok=ok),
              open(os.path.join(ROOT, "numerics_check.json"), "w"), indent=1)
    print("numerics check:", "OK" if ok else "FAIL", "(%.1fs)" % (time.time() - t0))
    sys.exit(0 if ok else 1)
