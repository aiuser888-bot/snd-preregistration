import os, sys, time, traceback, numpy as np
from snd.design import seed_to_key, make_block_rng, block_angles, design_check
from snd.stats import center_two_way, block_V
from snd.randtest import randomization_test
from snd.fcrice import TABLE_CHECKPOINTS, load_table
from snd.kappa import kappa_exact, coverage_prediction
from snd.accept import simulate
from snd.leak import (type1_of_beta, beta_of_B, beta_max, q_rms_for_beta, coverage_worst_phase,
                      coverage_min_over_nu, u_mon, CHI2_2_95, chi2_divergence, beta_ub_from_D, size_ub_from_D,
                      devices_D, hash_D_bound, seed_entropy_ok, best_of_k_size, retry_size_ub)
from snd.design import seed_material, MIN_RAW_BYTES
import itertools

HERE = os.path.dirname(__file__)
TABLE = os.path.join(HERE, "..", "..", "results", "fc_rice_table.npz")
KEY = seed_to_key(b"0123456789abcdef-test-seed-" + bytes(range(40)))


def test_seed_too_short():
    try:
        seed_to_key(b"short")
    except ValueError:
        return
    raise AssertionError("short seed accepted")


def test_reproducible_and_design_check():
    a = block_angles(make_block_rng(KEY, 3, 5), 64)[0]
    b = block_angles(make_block_rng(KEY, 3, 5), 64)[0]
    assert np.array_equal(a, b) and design_check(a)
    assert not np.array_equal(a, block_angles(make_block_rng(KEY, 3, 6), 64)[0])


def test_permutation_uniform_L4():
    counts = {}
    for m in range(1, 4801):
        p = tuple(block_angles(make_block_rng(KEY, 0, m), 4)[2])
        counts[p] = counts.get(p, 0) + 1
    assert len(counts) == 24
    chi2 = sum((c - 200) ** 2 / 200 for c in counts.values())
    assert chi2 < 52.6  # chi2_{23} 0.9995 quantile


def test_lemma1_variance():
    """Lemma 1: Var(Re U_b) = Var(Im U_b) = V_b, Cov = 0 for fixed centred f (equispaced + uniform permutation)."""
    rng = np.random.default_rng(1)
    L = 200
    f = 1.25 * (rng.random((1, L)) < 0.2) - 0.25
    ft = f - f.mean(axis=1, keepdims=True)      # any fixed vector with zero sum (Lemma 1 premise)
    Vb = block_V(ft)[0]
    U = np.array([(ft[0] * np.exp(1j * block_angles(make_block_rng(KEY, 0, m), L)[0])).sum() for m in range(1, 6001)])
    assert abs(U.real.var() / Vb - 1) < 0.06 and abs(U.imag.var() / Vb - 1) < 0.06
    assert abs(np.corrcoef(U.real, U.imag)[0, 1]) < 0.06


def test_randtest_size_under_H0():
    """Type-I error of the randomization test (M=99 -> exact level 0.05). Code check, not a proof (see protocol 7(A))."""
    rng = np.random.default_rng(2)
    K, L, R = 12, 40, 200
    rej = 0
    for r in range(R):
        ft = center_two_way(1.25 * (rng.random((K, L)) < 0.2) - 0.25)
        key = seed_to_key(b"size-test-%08d-" % r + bytes(range(50)))
        rej += randomization_test(ft, key, M=99)["reject"]
    se = np.sqrt(0.05 * 0.95 / R)
    assert rej / R <= 0.05 + 3.5 * se


def test_table_checkpoints():
    xs, ucb = load_table(TABLE)
    for x, v in TABLE_CHECKPOINTS.items():
        assert abs(np.interp(x, xs, ucb) - v) <= 0.01


def test_kappa_values():
    k0, _ = kappa_exact(100, 0.0)
    assert k0 == 1.0
    k9, _ = kappa_exact(100, 0.9)
    assert abs(k9 - 0.9007) < 0.0015
    km, _ = kappa_exact(100, -0.6)
    assert abs(km - 1.0015) < 0.0015
    for rho in (-0.9, -0.6, 0.6, 0.9):
        assert kappa_exact(100, rho)[0] <= 1 + 1 / 99 + 1e-9


def test_coverage_iid_and_serial():
    xs, ucb = load_table(TABLE)
    for rho in (0.0, 0.9):
        res = simulate(L=40, K=80, rho=rho, nu0=3.0, R=300, seed=7 + int(10 * rho), xs=xs, ucb=ucb)
        lo, hi = res["wilson_cov"]
        assert lo >= 0.93, (rho, res)           # nominal 0.95; true ~0.97 (conservative)
        kt, _ = kappa_exact(40, rho)
        assert abs(res["kappa_emp"] - kt) < 0.03
    # prediction function sanity (iid): between Rice minimum 0.9709 and 1
    assert 0.97 <= coverage_prediction(2.65, 1.0, xs, ucb) <= 1.0


def test_leak_curves():
    xs, ucb = load_table(TABLE)
    assert abs(CHI2_2_95 - 5.991464547) < 1e-6
    assert abs(float(type1_of_beta(0.0)) - 0.05) < 1e-12
    for b, ref in ((0.25, 0.0547), (0.5, 0.0693), (1.0, 0.1327), (2.0, 0.4154)):      # V2.0.08 table
        assert abs(float(type1_of_beta(b)) - ref) < 6e-5, (b, float(type1_of_beta(b)))
    for b in (0.1, 0.25, 0.5):                                                          # quadratic approximation
        assert abs(float(type1_of_beta(b)) - (0.05 + 0.075 * b * b)) < 6e-4
    c0, nu0 = coverage_min_over_nu(0.0, xs, ucb)
    assert abs(c0 - 0.9709) < 0.003 and abs(nu0 - 2.65) < 0.2        # Rice minimum, section 5.1
    assert abs(coverage_min_over_nu(1.0, xs, ucb)[0] - 0.8297) < 0.005
    assert coverage_worst_phase(3.0, 0.0, xs, ucb) > 0.96
    assert abs(q_rms_for_beta(0.25, 1.1e10, 10**4) - 1.69e-6) < 2e-8
    assert abs(beta_max(1.1e10, q_rms_for_beta(0.25, 1.1e10, 10**4), 10**4) - 0.25) < 1e-9


def test_leak_cauchy_schwarz_and_exact_bias():
    rng = np.random.default_rng(11)
    K, L = 30, 40
    for _ in range(200):
        ft = center_two_way(rng.normal(size=(K, L)))
        q = (rng.normal(size=(K, L)) + 1j * rng.normal(size=(K, L))) * rng.uniform(0, 0.1)
        b, om, qr = beta_of_B(ft, q)
        assert b <= beta_max(K * L, qr, L) * (1 + 1e-12) and 0 <= om <= 1 + 1e-12
        # equality when q is proportional to f~
        b2, om2, qr2 = beta_of_B(ft, 0.01 * ft.astype(complex))
        assert abs(om2 - 1) < 1e-9 and abs(b2 - beta_max(K * L, qr2, L)) < 1e-9


def test_u_mon_exact_size_and_unbiasedness():
    # (a) A1 holds: randomization p-value of |U_mon| has size <= 0.05 (exact up to Monte-Carlo error)
    rng = np.random.default_rng(5)
    K, L, M, R = 40, 40, 99, 300
    j = np.arange(L); base = 2 * np.pi * j / L
    cls = np.broadcast_to(2 * (np.arange(K)[:, None] % 2) + (j[None, :] >= L // 2), (K, L)).astype(int)
    def draw():
        a = rng.uniform(0, 2 * np.pi, (K, 1)); return np.take_along_axis(a + base, np.argsort(rng.random((K, L)), axis=1), axis=1)
    rej = 0
    for _ in range(R):
        ft = center_two_way(1.25 * (rng.random((K, L)) < 0.2 + 0.1 * (cls % 2) * (j[None, :] / L)) - 0.25)
        uo = abs(u_mon(ft, draw(), cls, 4))
        cnt = sum(abs(u_mon(ft, draw(), cls, 4)) >= uo for _ in range(M))
        rej += (1 + cnt) / (M + 1) <= 0.05
    assert 0.02 <= rej / R <= 0.085, rej / R           # 0.05 +- 3.5 SE (SE 0.0126 at R=300)
    # (b) E[U_mon | G] = sum_c F_c qbar_c exactly: leak theta->0 on class 1 with prob 0.3 => qbar_1 = 0.3
    ft = center_two_way(rng.normal(size=(K, L))); acc = 0; Rr = 400
    for _ in range(Rr):
        th = draw(); th = np.where((cls == 1) & (rng.random((K, L)) < 0.3), 0.0, th); acc += u_mon(ft, th, cls, 4)
    n = np.bincount(cls.ravel(), minlength=4); F = np.bincount(cls.ravel(), weights=ft.ravel(), minlength=4)
    pred = 0.3 * F[1]
    se = np.sqrt((F ** 2 / n).sum() / 2 / Rr) * 1.0
    assert abs(acc / Rr - pred) < 4 * se + 1e-9, (acc / Rr, pred, se)


def test_seed_material_rules():
    # 3.1 (V2.0.10): raw QRNG output >= 512 bits; key = first 128 bits of SHA-256 of raw || device || version
    try:
        seed_material(bytes(63), 1)
    except ValueError:
        pass
    else:
        raise AssertionError("63-byte raw material accepted")
    m1, m2 = seed_material(bytes(range(64)), 1), seed_material(bytes(range(64)), 2)
    assert len(m1) >= MIN_RAW_BYTES and not np.array_equal(seed_to_key(m1), seed_to_key(m2))
    assert seed_to_key(m1).shape == (2,) and seed_to_key(m1).dtype == np.uint64      # m = 128 bits
    try:
        seed_to_key(bytes(63))
    except ValueError:
        return
    raise AssertionError("short seed material accepted")


def _ub_table(f, alphas, L):
    return np.array([[(f * np.exp(1j * (a + 2 * np.pi * np.array(p) / L))).sum()
                      for p in itertools.permutations(range(L))] for a in alphas]).ravel()


def test_lemma5_design_moments_exact():
    # 12.5.3(a): over the design law (uniform rotation grid x all L! permutations) E U = 0, E U^2 = 0, E|U|^2 = 2V,
    # hence E[Re(conj(w) U)^2] = V for every unit w.  Exact enumeration, L = 4, 5, 6.
    rng = np.random.default_rng(11)
    for L in (4, 5, 6):
        f = rng.normal(size=L); f -= f.mean()
        U = _ub_table(f, 2 * np.pi * np.arange(7) / 7, L)
        V = L / (2.0 * (L - 1)) * (f ** 2).sum()
        assert abs(U.mean()) < 1e-12 and abs((U ** 2).mean()) < 1e-12
        assert abs((abs(U) ** 2).mean() / (2 * V) - 1) < 1e-12
        for w in (1.0, np.exp(0.7j), np.exp(-2.1j)):
            assert abs(((np.conj(w) * U).real ** 2).mean() / V - 1) < 1e-12


def test_lemma5_beta_le_sqrtD_tight_and_tv():
    # 12.5.3(b)(c): finite seed space (K=2 blocks, L=4, 6 rotation grid points -> 20736 seeds), exact.
    rng = np.random.default_rng(12)
    L, A = 4, 6
    al = 2 * np.pi * np.arange(A) / A
    fs = [rng.normal(size=L) for _ in range(2)]
    fs = [f - f.mean() for f in fs]
    V = sum(L / (2.0 * (L - 1)) * (f ** 2).sum() for f in fs)
    t0, t1 = _ub_table(fs[0], al, L), _ub_table(fs[1], al, L)
    U = (t0[:, None] + t1[None, :]).ravel()
    assert abs((abs(U) ** 2).mean() / (2 * V) - 1) < 1e-12
    for eps in (0.05, 0.2, 0.4):                                   # aligned tilt l-1 = eps Re(conj(w)U)/sqrt(V): equality
        l = 1 + eps * (np.conj(np.exp(0.9j)) * U).real / np.sqrt(V)
        assert l.min() >= 0
        D = chi2_divergence(l / l.size)
        beta = abs((l * U).mean()) / np.sqrt(V)
        assert abs(D - eps ** 2) < 1e-12 and abs(beta - beta_ub_from_D(D)) < 1e-12
        tv = 0.5 * np.abs(l - 1).mean()
        assert tv <= 0.5 * np.sqrt(D) + 1e-12
        assert size_ub_from_D(D) >= 0.05 + tv - 1e-12
    for _ in range(200):                                           # arbitrary tilts: bound holds, TV bound holds, event bound holds
        l = np.exp(rng.normal(size=U.size) * rng.uniform(0.01, 1.5)); l /= l.mean()
        D = chi2_divergence(l / l.size)
        assert abs((l * U).mean()) / np.sqrt(V) <= beta_ub_from_D(D) + 1e-12
        tv = 0.5 * np.abs(l - 1).mean()
        assert tv <= 0.5 * np.sqrt(D) + 1e-12
        ev = rng.random(U.size) < rng.uniform(0.05, 0.95)
        assert (l * ev).mean() - ev.mean() <= tv + 1e-12


def test_lemma5_hash_collision_equality_gf2():
    # 12.5.3(e): h(x) = A x over GF(2), A uniform lam x n (lam = m below): E_A D_out = (2^lam - 1) CP(R) exactly.  Exhaustive n=6, lam=3.
    n, m = 6, 3
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n)) & 1
    pr = np.prod(np.where(bits == 1, 0.7, 0.3), axis=1)
    CP = (pr ** 2).sum()
    par = np.array([[bin(x & r).count("1") & 1 for r in range(2 ** n)] for x in range(2 ** n)])
    tot = 0.0; R = np.arange(2 ** n)
    for r1 in range(2 ** n):
        for r2 in range(2 ** n):
            z = (par[:, r1] + 2 * par[:, r2])[None, :] + 4 * par[:, R].T
            pm = np.stack([((z == k) * pr[None, :]).sum(1) for k in range(2 ** m)], axis=1)
            tot += (2 ** m * (pm ** 2).sum(1) - 1).sum()
    assert abs(tot / 2 ** (3 * n) - (2 ** m - 1) * CP) < 1e-12
    assert tot / 2 ** (3 * n) <= hash_D_bound(-np.log2(CP), m) + 1e-12


def test_lemma5_numbers_and_best_of_k():
    D1 = hash_D_bound(128 + 32)                                    # registered requirement: H2 >= lambda + 32
    assert D1 == 2.0 ** -32 and seed_entropy_ok(512, 0.3125) and not seed_entropy_ok(512, 0.30)
    D4 = devices_D(D1, 4)
    assert abs(D4 / (4 * D1) - 1) < 1e-8
    assert abs(beta_ub_from_D(D1) - 2 ** -16) < 1e-18 and abs(size_ub_from_D(D1) - (0.05 + 2 ** -17)) < 1e-15
    assert [round(best_of_k_size(k), 4) for k in (1, 2, 3, 20)] == [0.05, 0.0975, 0.1426, 0.6415]
    # best-of-k size by simulation of the exact randomization test (K=4, L=8, M=199, reject iff #<=9)
    rng = np.random.default_rng(13)
    K, L, M, R = 4, 8, 199, 3000
    f = center_two_way(1.25 * (rng.random((K, L)) < 0.3) - 0.25)
    base = 2 * np.pi * np.arange(L) / L
    def U_abs(n):
        a = rng.uniform(0, 2 * np.pi, (n, K, 1)); th = np.take_along_axis(a + base, np.argsort(rng.random((n, K, L)), axis=2), axis=2)
        return np.abs((f * np.exp(1j * th)).sum(axis=(1, 2)))
    rej = np.zeros((R, 3), bool)
    for k in range(3):
        for c in range(0, R, 500):
            u0 = U_abs(500); um = np.stack([U_abs(500) for _ in range(M)], axis=1)
            rej[c:c + 500, k] = (um >= u0[:, None]).sum(1) <= 9
    p1, p2, p3 = rej[:, 0].mean(), (rej[:, 0] | rej[:, 1]).mean(), rej.any(axis=1).mean()
    for p, k in ((p1, 1), (p2, 2), (p3, 3)):
        se = np.sqrt(best_of_k_size(k) * (1 - best_of_k_size(k)) / R)
        assert abs(p - best_of_k_size(k)) < 4 * se, (k, p)


def test_lemma5_retry_on_theta_only_failure():
    # 12.5.3(d'): finite toy.  1000 seeds, 'pass' set of size 900 (f = 0.1); A subset of pass with P_0(A) = 0.05 (worst case).
    # Sequential redraw (<=3 draws), first passing draw used: exact law P_0(s)(1+f+f^2) on pass; conditional size = 0.05/(1-f).
    n, f = 1000, 0.1
    pass_ = np.arange(n) < 900; A = np.arange(n) < 50
    p_use = np.where(pass_, (1 + f + f ** 2) / n, 0.0)
    p_analysis = p_use.sum()
    assert abs(p_analysis - (1 - f ** 3)) < 1e-12
    assert abs((p_use * A).sum() / p_analysis - retry_size_ub(f)) < 1e-12
    rng = np.random.default_rng(14); R = 200000; hit = ana = 0
    for _ in range(R):
        for _k in range(3):
            s = rng.integers(n)
            if pass_[s]:
                ana += 1; hit += A[s]; break
    assert abs(ana / R - (1 - f ** 3)) < 4 * np.sqrt(f ** 3 / R) + 1e-9
    assert abs(hit / ana - retry_size_ub(f)) < 4 * np.sqrt(retry_size_ub(f) * (1 - retry_size_ub(f)) / ana)
    assert abs(retry_size_ub(1e-3) - 0.05005005) < 1e-8


def test_theta_is_pure_function_of_seed_and_counter():
    # E1 (12.5.1 step 5): theta = Theta(seed material, device, version, (block, m')) -- no other input.  Recompute via the raw
    # numpy Philox API and with perturbed global RNG state: bit-identical; a different device id or version gives a different key.
    raw = bytes(range(64)); mat = seed_material(raw, 3)
    key = seed_to_key(mat); L = 32
    for b in (0, 5, 123):
        for m in (0, 1, 7):
            th = block_angles(make_block_rng(key, b, m), L)[0]
            g = np.random.Generator(np.random.Philox(key=key, counter=np.array([m, b, 0, 0], dtype=np.uint64)))
            alpha = g.uniform(0.0, 2 * np.pi); perm = g.permutation(L)
            assert np.array_equal(th, (alpha + 2 * np.pi * np.arange(L) / L)[perm])
    np.random.seed(1); a = block_angles(make_block_rng(key, 9, 2), L)[0]
    np.random.seed(2); _ = np.random.random(1000); b2 = block_angles(make_block_rng(key, 9, 2), L)[0]
    assert np.array_equal(a, b2)
    assert not np.array_equal(seed_to_key(seed_material(raw, 4)), key)
    assert not np.array_equal(seed_to_key(seed_material(raw, 3, "2.0.11")), key)


if __name__ == "__main__":      # dependency-free runner (pytest also works)
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            t = time.time()
            try:
                fn(); print("PASS %-32s %.1fs" % (name, time.time() - t))
            except Exception:
                fails += 1; print("FAIL", name); traceback.print_exc()
    sys.exit(1 if fails else 0)
