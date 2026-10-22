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
    import snd
    assert seed_material(raw, 3) == seed_material(raw, 3, snd.__version__)          # the default version tag IS the package version
    assert not np.array_equal(seed_to_key(seed_material(raw, 3, "0.0.0-other")), key)   # a different version string gives a different key


# ---------------------------------------------------------------------------------------------------------------------------
# V2.0.13: normative theta algorithm (protocol 3.1), fixed-order numerics (7(D)), reproducibility of the re-draws (4)
# ---------------------------------------------------------------------------------------------------------------------------
import json, math, random as _pyrandom, hashlib, runpy, re
from snd import philox_ref as P
from snd import design as D
from snd import numerics as NU
from snd.randtest import draw_theta, tail_counts, split_ms, abs_U_of_draw

VEC = os.path.join(HERE, "vectors", "theta_kat.json")
MAXW = (1 << 64) - 1


def _bits(x):
    return np.asarray(x, dtype="<f8").tobytes()


def test_philox_cipher_known_answers():
    """Philox4x64-10 in pure Python against: the ISO C++ std::philox4x64 value; numpy's published Philox test data; numpy's Philox for random
    keys and counters, including counter carries.  These are the independent checks of the cipher itself."""
    assert P.philox4x64_10((2499, 0, 0, 0), (20111115, 0))[3] == 3409172418970261260       # C++ engine: output #10000 (counter increments after a block)
    data_dir = os.path.join(os.path.dirname(np.random.__file__), "tests", "data")
    seen = 0
    for fn in ("philox-testset-1.csv", "philox-testset-2.csv"):
        path = os.path.join(data_dir, fn)
        if not os.path.exists(path):
            continue
        rows = open(path).read().strip().splitlines()
        seed = int(rows[0].split(",")[1], 16); vals = [int(r.split(",")[1], 16) for r in rows[1:]]
        st = np.random.Philox(seed).state["state"]
        s = P.PhiloxStream([int(k) for k in st["key"]], [int(c) for c in st["counter"]])
        assert [s.next_u64() for _ in range(len(vals))] == vals
        seen += 1
    rnd = _pyrandom.Random(99)
    cases = [((0, 0), (0, 0, 0, 0)), ((MAXW, MAXW), (MAXW, 0, 0, 0)), ((1, 2), (MAXW, MAXW, MAXW, MAXW)), ((5, 6), (MAXW, MAXW, 0, 7))]
    cases += [((rnd.getrandbits(64), rnd.getrandbits(64)), tuple(rnd.getrandbits(64) for _ in range(4))) for _ in range(40)]
    for key, ctr in cases:
        raw = [int(v) for v in np.random.Philox(key=np.array(key, dtype=np.uint64), counter=np.array(ctr, dtype=np.uint64)).random_raw(12)]
        s = P.PhiloxStream(key, ctr)
        assert raw == [s.next_u64() for _ in range(12)], (key, ctr)
    print("   numpy test-data files used:", seen, end=" ")


def test_theta_matches_pure_python_reference():
    """E1 (12.5.1 step 5): the numpy implementation and the integer-only reference give bit-identical alpha, permutation and angles."""
    rnd = _pyrandom.Random(2026)
    n = 0
    for t in range(40):
        raw = bytes(rnd.getrandbits(8) for _ in range(64))
        mat = seed_material(raw, rnd.randrange(1, 5))
        key = seed_to_key(mat)
        import struct
        assert key.dtype == np.uint64 and tuple(int(k) for k in key) == struct.unpack("<QQ", hashlib.sha256(mat).digest()[:16])
        assert [int(k) for k in key] == list(P.key_from_digest(hashlib.sha256(mat).digest()))
        for L in (1, 2, 3, 7, 10, 64, 100, 257, 1000):
            b = rnd.choice([0, 1, 5, 99, 2 ** 20, 2 ** 40]); m = rnd.choice([0, 1, 2, 1999])
            th, a, p = block_angles(make_block_rng(key, b, m), L)
            th2, a2, p2 = P.block_theta([int(k) for k in key], b, m, L)
            assert a == a2 and list(p) == p2 and _bits(th) == P.theta_bytes(th2), (t, L, b, m)
            n += 1
    key = seed_to_key(seed_material(bytes(range(64)), 1))
    th, a, p = block_angles(make_block_rng(key, 3, 7), 10000)
    th2, a2, p2 = P.block_theta([int(k) for k in key], 3, 7, 10000)
    assert a == a2 and list(p) == p2 and _bits(th) == P.theta_bytes(th2)
    assert D.stream_conformance() is True
    ref = hashlib.sha256()
    for b in (0, 3, 9):
        ref.update(P.theta_bytes(P.block_theta([int(k) for k in key], b, 0, 50)[0]))
    assert D.theta_digest(key, [9, 0, 3], 50) == ref.hexdigest()                    # order of the list does not matter: ascending block number
    assert D.theta_digest(key, [0, 3, 9], 50, m=1) != ref.hexdigest()
    print("   cases:", n + 1, end=" ")


def test_known_answer_vectors_file_is_current():
    """tests/vectors/theta_kat.json is exactly what scripts/make_kat.py generates, and the numpy implementation reproduces every entry."""
    gen = runpy.run_path(os.path.join(HERE, "..", "scripts", "make_kat.py"))["generate"]()
    text = json.dumps(gen, indent=1, sort_keys=True) + "\n"
    assert open(VEC).read() == text, "run: python scripts/make_kat.py"
    ext = gen["external_cpp_standard"]
    assert ext["computed_hex"] == ext["value_hex"] and int(ext["value_hex"], 16) == ext["value_dec"] == 3409172418970261260
    w = lambda x: int(x, 16)
    for e in gen["cipher"]:
        assert P.philox4x64_10([w(c) for c in e["counter"]], [w(k) for k in e["key"]]) == tuple(w(o) for o in e["out"])
    for e in gen["stream"]:
        key, ctr = [w(k) for k in e["key"]], [w(c) for c in e["counter"]]
        raw = np.random.Philox(key=np.array(key, dtype=np.uint64), counter=np.array(ctr, dtype=np.uint64)).random_raw(8)
        assert [int(v) for v in raw] == [w(x) for x in e["words"]]
    for e in gen["seed_to_key"]:
        mat = bytes.fromhex(e["seed_material_hex"])
        assert hashlib.sha256(mat).hexdigest() == e["sha256"] and [int(k) for k in seed_to_key(mat)] == [w(k) for k in e["key"]]
    for e in gen["block_theta"]:
        key = np.array([w(k) for k in e["key"]], dtype=np.uint64)
        th, a, p = block_angles(make_block_rng(key, e["block"], e["m"]), e["L"])
        assert int.from_bytes(np.float64(a).tobytes(), "little") == w(e["alpha_bits"]) and list(p) == e["perm"]
        assert [int.from_bytes(t.tobytes(), "little") for t in th.astype("<f8")] == [w(x) for x in e["theta_bits"]]
    for e in gen["theta_digest"]:
        key = np.array([w(k) for k in e["key"]], dtype=np.uint64)
        assert D.theta_digest(key, e["blocks"], e["L"], e["m"]) == e["sha256"]


def _ref_tree(v):
    n = 1 << max(0, (len(v) - 1).bit_length())
    v = list(v) + [0.0] * (n - len(v))
    while len(v) > 1:
        h = len(v) // 2
        v = [v[i] + v[i + h] for i in range(h)]
    return v[0]


def test_numerics_fixed_order_and_sum_checks():
    """7(D): the tree sum equals an independent scalar implementation of the documented order bit for bit and does not depend on memory
    layout; Kahan, tree, double-double and exact sums agree within the 7(D) tolerances; cos/sin library error is within the assumed 4 ulp."""
    rng = np.random.default_rng(5)
    for n in (1, 2, 3, 5, 8, 13, 1000):
        x = rng.normal(size=n) * 10.0 ** rng.integers(-3, 4, size=n)
        assert NU.tree_sum(x) == _ref_tree(x.tolist())
        z = x + 1j * x[::-1]
        assert NU.tree_sum(z) == complex(_ref_tree(z.real.tolist()), _ref_tree(z.imag.tolist()))
    A = rng.normal(size=(7, 37))
    assert np.array_equal(NU.tree_sum(A, axis=1), NU.tree_sum(np.asfortranarray(A), axis=1))
    assert np.array_equal(NU.tree_sum(A, axis=0), NU.tree_sum(A.T.copy(), axis=1))
    assert NU.tree_sum(np.zeros((3, 0)), axis=1).shape == (3,)
    hard = [1e16, 1.0, -1e16, 1.0] * 50 + [3.0]                                       # exact sum 103
    assert NU.exact_sum(hard) == 103.0 and NU.dd_sum(hard) == 103.0
    K, L = 20, 500
    ft = center_two_way(1.25 * (rng.random((K, L)) < 0.2) - 0.25)
    th = draw_theta(KEY, K, L, 0)
    r = NU.sum_check(ft, th)
    assert r["ok"] and r["dU_rel"] < 1e-12 and r["dV_rel"] < 1e-12, r
    assert abs(NU.u_stat(ft, th) - complex(NU.exact_sum((ft * np.cos(th)).reshape(-1).tolist()),
                                           NU.exact_sum((ft * np.sin(th)).reshape(-1).tolist()))) < 1e-9
    ang = rng.uniform(0.0, 4 * np.pi, 100000)
    for fn, ref in ((np.cos, math.cos), (np.sin, math.sin)):
        r0 = np.array([ref(a) for a in ang.tolist()])
        ulp = np.abs(fn(ang) - r0) / np.spacing(np.abs(r0))
        assert ulp.max() <= 4, (fn.__name__, ulp.max())


def test_randtest_independent_of_workers_and_order():
    """Item 4 + 3.1: the test result does not depend on the number of workers, on the partition or order of the re-draws, or on any global
    random state; a re-draw evaluated in another process gives the same bits."""
    rng = np.random.default_rng(11)
    K, L, M = 5, 24, 60
    ft = center_two_way(1.25 * (rng.random((K, L)) < 0.2) - 0.25)
    key = seed_to_key(b"order-test-" + bytes(range(60)))
    keys = ("p_hat", "tail", "M", "U_obs", "V", "Lambda_obs", "reject", "near_ties", "near_tie_draws", "p_hat_lo", "p_hat_hi", "uncertain", "decision",
            "error_bound", "p_hat_float", "resolution")
    base = randomization_test(ft, key, M=M)
    np.random.seed(123); _pyrandom.seed(7); _ = np.random.random(500)               # poisoned global state
    r2 = randomization_test(ft, key, M=M, chunks=[range(31, 61), range(1, 11), range(11, 31)])
    r3 = randomization_test(ft, key, M=M, workers=2)
    r4 = randomization_test(ft, key, M=M, chunks=[[m] for m in range(M, 0, -1)])
    for r in (r2, r3, r4):
        assert all(r[k] == base[k] for k in keys), [k for k in keys if r[k] != base[k]]
    assert base["p_hat"] == (1 + base["tail"]) / (M + 1) and base["near_ties"] == 0
    assert tail_counts(ft, key, list(range(1, M + 1)), abs(base["U_obs"]), 0.0)[0] == base["tail"]
    assert sum(len(c) for c in split_ms(M, 7)) == M and sorted(m for c in split_ms(M, 7) for m in c) == list(range(1, M + 1))
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    local = [abs_U_of_draw(ft, key, m) for m in (1, 17, 60)]
    with ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context("spawn")) as ex:
        remote = list(ex.map(abs_U_of_draw, [ft] * 3, [key] * 3, (1, 17, 60)))
    assert local == remote                                                          # identical floats, not merely close
    # 7(D) uncertainty report: a huge tolerance makes every re-draw a near tie -> decision uncertain; zero tolerance -> p_lo = p_hat = p_hi
    big = randomization_test(ft, key, M=M, tol_u=1e9, resolve=False)
    assert big["uncertain"] and big["decision"] == "uncertain" and big["near_ties"] == M
    zero = randomization_test(ft, key, M=M, tol_u=0.0)
    assert zero["p_hat_lo"] == zero["p_hat"] == zero["p_hat_hi"] and not zero["uncertain"]
    for bad in (lambda: randomization_test(ft, key, M=M, chunks=[range(1, M)]),                 # not a partition of 1..M
                lambda: randomization_test(ft, key, M=M, chunks=[range(1, M + 1), range(1, 2)]),    # overlapping
                lambda: randomization_test(ft, key, M=0),
                lambda: make_block_rng(np.array([1, 2, 3], dtype=np.uint64), 0, 0),                 # wrong key size
                lambda: make_block_rng(key, -1, 0), lambda: make_block_rng(key, 0, 2 ** 64),
                lambda: D.theta_digest(key, [1, 1], 8)):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("invalid input accepted")



def test_ucb95_and_scores():
    """5.2: UCB_95 = sigma_hat * u(x_hat) (table below 18, x + 1.96 above, never below the table); 2: Lambda = |U|^2 / V and the fixed-order scores."""
    from snd.fcrice import u_of_x, ucb95, XMAX_VALID
    from snd.stats import block_scores, lam
    xs, ucb = load_table(TABLE)
    for x, v in TABLE_CHECKPOINTS.items():
        assert abs(float(u_of_x(x, xs, ucb)) - v) <= 0.01
    assert XMAX_VALID == 18.0 and float(u_of_x(18.0, xs, ucb)) == float(np.interp(18.0, xs, ucb))
    assert float(u_of_x(18.5, xs, ucb)) == 18.5 + 1.96 and float(u_of_x(30.0, xs, ucb)) == 31.96
    assert np.interp(18.0, xs, ucb) <= 18.0 + 1.96 + 1e-9                           # the tail formula is not below the table at the switch
    grid = np.linspace(0.05, 18.0, 400)
    assert np.all(np.diff(u_of_x(grid, xs, ucb)) >= -1e-9)                           # u(x) is non-decreasing
    U, V = 3 + 4j, 4.0                                                               # x_hat = 5 / 2
    assert abs(ucb95(U, V, xs, ucb) - float(np.interp(2.5, xs, ucb)) / 2.0) < 1e-12
    c = 7.0
    assert abs(ucb95(c * U, c * c * V, xs, ucb) - ucb95(U, V, xs, ucb) / c) < 1e-12       # scale: x_hat invariant, sigma_hat -> sigma_hat / c
    rng = np.random.default_rng(3)
    ft = center_two_way(1.25 * (rng.random((9, 30)) < 0.2) - 0.25)
    th = draw_theta(KEY, 9, 30, 0)
    Us, Vs = block_scores(ft, th)
    assert abs(Us - (ft * np.exp(1j * th)).sum()) < 1e-12 and abs(Vs - block_V(ft).sum()) < 1e-12
    assert abs(lam(Us, Vs) - abs(Us) ** 2 / Vs) < 1e-15


def test_trace_labels_parse():
    """scripts/check_trace.py: the reference grammar parses as documented, and existing / non-existing protocol labels are told apart."""
    ct = runpy.run_path(os.path.join(HERE, "..", "scripts", "check_trace.py"))
    root = os.path.join(HERE, "..", "..", "protocol")
    import snd
    path = os.path.join(root, "SND_protocol_V%s.md" % snd.__version__)
    P = ct["read_protocol"](path)
    assert ct["parse_refs"]("protocol 7(D), Lemma 5(c)(f), protocol 12.5.1 step 3, Lemma 5(d'), Lemma 2 (iv)") == [
        ("sec", "7", "D", None), ("sec", "12.5.1", None, 3), ("lem", 5, ["c", "f"], None), ("lem", 5, ["d′"], None), ("lem", 2, ["iv"], None)]
    assert ct["validate"](ct["parse_refs"]("protocol 3.1, 4, 5.1, 7(B′), 12.5.3, Lemma 1, Lemma 2(vi), Lemma 3(a), Lemma 4(g), Lemma 5(d'), protocol 12.5.1 step 5"), P) == []
    for bad in ("protocol 99", "protocol 7(Z)", "protocol 12.5.1 step 9", "Lemma 9", "Lemma 5(z)", "Lemma 1(a)", "protocol 3.9"):
        assert ct["validate"](ct["parse_refs"](bad), P), bad
    assert P["seven"] == {"A", "B", "B′", "C", "D", "E", "F"} and P["steps_1251"] == 5 and sorted(P["lemma_items"]) == [1, 2, 3, 4, 5]
    assert {"a", "b", "c", "d", "d′", "e", "f", "g"} <= P["lemma_items"][5]


def test_fc_rice_table_rebuilds_identically():
    """5.1: the published table is exactly what build_table() computes (arrays identical), with the acceptance mass in [0.9500, 0.9505]
    (the file's own SHA-256 is locked in SHA256SUMS.txt; a rebuilt .npz differs byte-wise only in zip timestamps)."""
    from snd.fcrice import build_table
    xs, ucb, mass = build_table()
    xs0, ucb0 = load_table(TABLE)
    assert np.array_equal(xs, xs0) and np.array_equal(ucb, ucb0)
    assert len(xs) == 5499 and abs(xs[-1] - 21.996) < 1e-9 and abs(xs[0] - 0.004) < 1e-12
    assert 0.9500 <= mass.min() and mass.max() <= 0.9505, (mass.min(), mass.max())


def test_version_checks_detect_a_stale_label():
    """check_trace.version_checks: the repository is consistent, and each kind of stale label (version, test count, copy, extra protocol file,
    second 'current version' marker) is reported."""
    import shutil, tempfile
    import snd
    ct = runpy.run_path(os.path.join(HERE, "..", "scripts", "check_trace.py"))
    root0 = ct["ROOT"]; v = snd.__version__
    ntests = len(re.findall(r"^def (test_\w+)", open(os.path.join(HERE, "test_snd.py"), encoding="utf8").read(), re.M))
    assert all(ok for ok, _ in ct["version_checks"](root0, v, ntests)), [d for ok, d in ct["version_checks"](root0, v, ntests) if not ok]
    files = ["CITATION.cff", ".zenodo.json", "metadata/zenodo.json", "protocol/SND_protocol_V%s.md" % v, "VERSIONING.md", "README.md", "CHANGELOG.md",
             "docs/CHANGELOG.md", "docs/audits/SND_audit_V1.0.08_to_V%s.md" % v, "docs/SND_evolution_for_everyone.md", "docs/SND_evolution_summary.md",
             "OSF_REGISTRATION_GUIDE.md"]

    def copy_to(dst):
        for rel in files:
            os.makedirs(os.path.dirname(os.path.join(dst, rel)), exist_ok=True)
            shutil.copy(os.path.join(root0, rel), os.path.join(dst, rel))

    def edit(dst, rel, a, b):
        path = os.path.join(dst, rel); t = open(path, encoding="utf8").read()
        assert a in t, (rel, a)
        open(path, "w", encoding="utf8", newline="").write(t.replace(a, b, 1))

    mutations = [("CITATION.cff", 'version: "%s"' % v, 'version: "0.0.0"'),
                 ("README.md", "# %d PASS" % ntests, "# %d PASS" % (ntests - 1)),
                 ("docs/CHANGELOG.md", "### V%s" % v, "### Vx"),
                 ("VERSIONING.md", "**V%s**" % v, "**V0.0.0**"),
                 (".zenodo.json", '"version": "%s"' % v, '"version": "0.0.0"'),
                 ("CHANGELOG.md", "\uff08", "\uff08\u76ee\u524d\u7248\u672c\uff09\uff08")]
    for rel, a, b in mutations:
        with tempfile.TemporaryDirectory() as d:
            copy_to(d)
            assert all(ok for ok, _ in ct["version_checks"](d, v, ntests))
            edit(d, rel, a, b)
            assert not all(ok for ok, _ in ct["version_checks"](d, v, ntests)), (rel, a)
    with tempfile.TemporaryDirectory() as d:                      # a second protocol file next to the current one
        copy_to(d)
        shutil.copy(os.path.join(d, "protocol", "SND_protocol_V%s.md" % v), os.path.join(d, "protocol", "SND_protocol_V0.0.0.md"))
        assert not all(ok for ok, _ in ct["version_checks"](d, v, ntests))
    # personal-identifier check: the word list is never in the repository; here it is a made-up list
    dummy = ["zq-sec" "ret-one", "Zq-Sec" "ret-Two"]                      # written in pieces: the scan also reads this file
    assert all(ok for ok, _ in ct["identity_checks"](root0, dummy))
    for word in dummy:                                            # a listed word in any file is reported, whatever its case
        with tempfile.TemporaryDirectory() as d:
            copy_to(d)
            assert all(ok for ok, _ in ct["identity_checks"](d, dummy))
            with open(os.path.join(d, "README.md"), "a", encoding="utf8") as f:
                f.write("\nauthor: " + word.swapcase() + "\n")
            bad = ct["identity_checks"](d, dummy)
            assert not all(ok for ok, _ in bad) and "README.md" in bad[0][1], word
    old_env = os.environ.get("SND_IDENTITY_WORDS")
    try:                                                          # no list available -> reported as skipped, not as passed silently
        os.environ["SND_IDENTITY_WORDS"] = os.path.join(root0, "no-such-file")
        skipped = ct["identity_checks"](root0)
        assert len(skipped) == 1 and skipped[0][0] and skipped[0][1].startswith("SKIPPED")
    finally:
        if old_env is None:
            os.environ.pop("SND_IDENTITY_WORDS", None)
        else:
            os.environ["SND_IDENTITY_WORDS"] = old_env
        pass


def test_checksum_tool_detects_every_kind_of_stale_list():
    """scripts/make_checksums.py: a list written after the last edit verifies; a changed, added, removed or reordered file, a duplicate
    line or a malformed line is reported; caches, .git/, .github/ and .gitignore are not part of the list."""
    import tempfile
    mc = runpy.run_path(os.path.join(HERE, "..", "scripts", "make_checksums.py"))
    with tempfile.TemporaryDirectory() as root:
        def put(rel, data=b"x"):
            path = os.path.join(root, *rel.split("/")); os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "wb").write(data)
        for rel in ("b.md", "a.md", "docs/z.md", "docs/a.md", "code/x.py"):
            put(rel, rel.encode())
        for rel in (".github/workflows/w.yml", ".gitignore", "code/__pycache__/x.pyc", ".git/HEAD", "code/y.pyc"):
            put(rel)
        assert mc["main"](["--write", "--root", root]) == 0
        lines = open(os.path.join(root, "SHA256SUMS.txt")).read().splitlines()
        assert [l.split("  ")[1] for l in lines] == ["./a.md", "./b.md", "./code/x.py", "./docs/a.md", "./docs/z.md"]
        assert mc["check"](root) == []
        put("docs/a.md", b"changed")
        assert mc["check"](root) == ["hash differs: docs/a.md"]
        mc["main"](["--write", "--root", root]); assert mc["check"](root) == []
        put("new.md"); assert mc["check"](root) == ["not listed: new.md"]
        os.remove(os.path.join(root, "new.md")); os.remove(os.path.join(root, "b.md"))
        assert mc["check"](root) == ["listed but missing: b.md"]
        mc["main"](["--write", "--root", root])
        path = os.path.join(root, "SHA256SUMS.txt"); good = open(path).read().splitlines()
        open(path, "w").write("\n".join(reversed(good)) + "\n")
        assert any("not sorted" in p for p in mc["check"](root))
        open(path, "w").write("\n".join(good + good[:1]) + "\n")
        assert any("listed twice" in p for p in mc["check"](root))
        open(path, "w").write("\n".join(good + ["garbage"]) + "\n")
        assert any("not '<sha256>" in p for p in mc["check"](root))
        os.remove(path)
        assert mc["check"](root) == ["SHA256SUMS.txt is missing"]


# ---------------------------------------------------------------------------------------------------------------------------
# V2.0.13: integer reference value, near ties and the second-tier resolution (protocol 4, 7(D), 10, 12.2), sealed f~ (protocol 2)
# ---------------------------------------------------------------------------------------------------------------------------
import shutil, subprocess, functools
from decimal import Decimal, localcontext, getcontext
from fractions import Fraction
from snd import highprec as HP
from snd import randtest as RT
from snd import decide as DE
from snd.stats import ftilde_digest
from snd.fcrice import u_of_x

REF_VEC = os.path.join(HERE, "vectors", "reference_kat.json")
getcontext().prec = 500          # comparisons below subtract numbers with up to 330 significant digits


@functools.lru_cache(maxsize=None)
def _dec_pi(prec):
    """pi by Machin's formula in decimal arithmetic (independent of highprec.pi_fixed, which works in binary integers)."""
    with localcontext() as ctx:
        ctx.prec = prec + 12

        def arccot(x):
            x = Decimal(x); x2 = x * x; term = 1 / x; s = term; n = 1; sign = -1
            while term > Decimal(10) ** (-(prec + 10)):
                term /= x2
                s += sign * term / (2 * n + 1)
                sign, n = -sign, n + 1
            return s
        return 16 * arccot(5) - 4 * arccot(239)


@functools.lru_cache(maxsize=None)
def _dec_sincos(theta, prec=90):
    """cos and sin of the exact binary value of `theta` in decimal arithmetic with `prec` digits (independent of highprec.sincos_fixed)."""
    with localcontext() as ctx:
        ctx.prec = prec + 12
        x, pi = Decimal(theta), _dec_pi(prec + 12)
        k = (x / (pi / 2)).to_integral_value()
        r = x - k * pi / 2
        r2, c, s, tc, ts, j = r * r, Decimal(1), r, Decimal(1), r, 1
        eps = Decimal(10) ** (-(prec + 6))
        while abs(tc) > eps or abs(ts) > eps:
            tc = -tc * r2 / ((2 * j - 1) * (2 * j))
            ts = -ts * r2 / ((2 * j) * (2 * j + 1))
            c, s, j = c + tc, s + ts, j + 1
        return ((c, s), (-s, c), (-c, -s), (s, -c))[int(k) % 4]


def _dec_U(w, th, classes=None, prec=60):
    """Exact-in-decimal U (or U_mon) of binary64 arrays: sum w_i (cos + i sin)(theta_i), with the class means of protocol 12.2 when `classes` is given."""
    wl, tl = np.asarray(w, float).ravel().tolist(), np.asarray(th, float).ravel().tolist()
    with localcontext() as ctx:
        ctx.prec = prec + 12
        if classes is None:
            cs = [_dec_sincos(t, prec) for t in tl]
            return sum((Decimal(a) * c for a, (c, s) in zip(wl, cs)), Decimal(0)), sum((Decimal(a) * s for a, (c, s) in zip(wl, cs)), Decimal(0))
        cl = np.asarray(classes[0]).ravel().tolist()
        re = im = Decimal(0)
        for c in sorted(set(cl)):
            idx = [i for i, v in enumerate(cl) if v == c]
            F = sum((Decimal(wl[i]) for i in idx), Decimal(0))
            re += F * sum((_dec_sincos(tl[i], prec)[0] for i in idx), Decimal(0)) / len(idx)
            im += F * sum((_dec_sincos(tl[i], prec)[1] for i in idx), Decimal(0)) / len(idx)
        return re, im


def _dec_abs2(w, th, classes=None):
    re, im = _dec_U(w, th, classes)
    return re * re + im * im


def test_highprec_sincos_against_independent_references():
    """7(D): the integer cos/sin is within one unit of 2^-P (in fact within 0.55; measured 0.50) of independent decimal-arithmetic values, of GNU bc if installed, for
    edge angles, random angles in [0, 4 pi) and beyond, and several precisions; pi likewise; invalid angles are rejected."""
    for bits in (64, 128, 256, 1024):
        assert abs(Decimal(HP.pi_fixed(bits)) - _dec_pi(400) * Decimal(2) ** bits) <= 2
    rnd = _pyrandom.Random(31)
    angles = [0.0, -0.0, 5e-324, 1e-300, 1e-10, 0.5, 1.0, math.pi / 2, math.pi, 2 * math.pi, 3 * math.pi / 2, 4 * math.pi, 12.566370614359172, -1.0, -math.pi,
              1e5, 2.0 ** 23, float(np.nextafter(math.pi / 2, 0.0)), float(np.nextafter(math.pi / 2, 9.0))]
    angles += [rnd.uniform(0, 4 * math.pi) for _ in range(120)] + [rnd.uniform(-1000, 1000) for _ in range(20)]
    worst = 0.0
    for t in angles:
        c, s = HP.sincos_fixed(t, 128)
        dc, ds = _dec_sincos(t)
        worst = max(worst, float(abs(Decimal(c) - dc * Decimal(2) ** 128)), float(abs(Decimal(s) - ds * Decimal(2) ** 128)))
        assert abs(c * c + s * s - (1 << 256)) <= 4 * (1 << 128)
    assert worst <= 0.55, worst
    assert HP.sincos_fixed(0.0, 128) == (1 << 128, 0)
    for P_, ts in ((20, angles[:30]), (300, angles[:6]), (1024, [0.5, 1e5, 3.0])):
        for t in ts:
            c, s = HP.sincos_fixed(t, P_)
            dc, ds = _dec_sincos(t, 330 if P_ == 1024 else 100)
            assert abs(Decimal(c) - dc * Decimal(2) ** P_) <= 1 and abs(Decimal(s) - ds * Decimal(2) ** P_) <= 1, (t, P_)
    if shutil.which("bc"):
        sel = angles[:60]
        inp = "scale=70\n" + "".join("s(%s)\nc(%s)\n" % (format(Decimal(a), "f"), format(Decimal(a), "f")) for a in sel)
        out = subprocess.run(["bc", "-l"], input=inp, capture_output=True, text=True).stdout.replace("\\\n", "").split()
        assert len(out) == 2 * len(sel)
        for i, a in enumerate(sel):
            c, s = HP.sincos_fixed(a, 128)
            assert abs(Decimal(s) - Decimal(out[2 * i]) * Decimal(2) ** 128) <= 0.55 and abs(Decimal(c) - Decimal(out[2 * i + 1]) * Decimal(2) ** 128) <= 0.55
        print("   bc used:", len(sel), "angles", end=" ")
    for bad in (float("nan"), float("inf"), 2.0 ** 24, -1e30):
        try:
            HP.sincos_fixed(bad, 64)
        except ValueError:
            continue
        raise AssertionError("invalid angle accepted")
    print("   worst error vs decimal: %.3f units" % worst, end=" ")


def test_highprec_reference_sum_and_comparison():
    """7(D), 4: exact integer representation of floats; the enclosure of U and of U_mon contains the independent decimal value and is tight;
    the floating-point statistics lie inside it up to the proven bound; comparisons of enclosures are decided only when certain."""
    rnd = _pyrandom.Random(5)
    vals = [0.0, 5e-324, -3.5, 1e300, 1e-300, 1.0, 0.1, -2.0 ** -40] + [rnd.uniform(-6, 6) * 2.0 ** rnd.randrange(-30, 30) for _ in range(60)]
    ints, e0 = HP.exact_weights(vals)
    assert all(Fraction(v) == Fraction(i) * Fraction(2) ** e0 for v, i in zip(vals, ints))
    for bad in (float("nan"), float("inf")):
        try:
            HP.exact_weights([1.0, bad])
        except ValueError:
            continue
        raise AssertionError("non-finite weight accepted")
    rng = np.random.default_rng(8)
    K, L = 3, 7
    ft = center_two_way(1.25 * (rng.random((K, L)) < 0.3) - 0.25)
    th = draw_theta(KEY, K, L, 0)
    re, im, err = HP.reference_U(ft, th, 128)
    dre, dim = _dec_U(ft, th)
    assert abs(Decimal(re.numerator) / Decimal(re.denominator) - dre) <= Decimal(err.numerator) / Decimal(err.denominator)
    assert abs(Decimal(im.numerator) / Decimal(im.denominator) - dim) <= Decimal(err.numerator) / Decimal(err.denominator)
    assert err < Fraction(1, 2 ** 100) * float(np.abs(ft).sum()) * 2 and err > 0
    ub = NU.u_error_bound(ft)
    U_fl = NU.u_stat(ft, th)
    assert abs(U_fl - complex(float(re), float(im))) <= ub * math.sqrt(NU.v_stat(ft)) and ub < 1e-6
    cls = (np.arange(K * L) % 4).reshape(K, L)
    plan = NU.class_plan(cls, 5)                                                       # label 4 is empty
    assert [len(p) for p in plan] == [6, 5, 5, 5] and all((np.diff(p) > 0).all() for p in plan)
    mre, mim, merr = HP.reference_U(ft, th, 128, (cls, 5))
    dmre, dmim = _dec_U(ft, th, (cls, 5))
    assert abs(Decimal(mre.numerator) / Decimal(mre.denominator) - dmre) <= Decimal(merr.numerator) / Decimal(merr.denominator)
    assert abs(Decimal(mim.numerator) / Decimal(mim.denominator) - dmim) <= Decimal(merr.numerator) / Decimal(merr.denominator)
    Um = NU.u_mon_stat(ft, th, plan)
    assert abs(Um - complex(float(mre), float(mim))) <= NU.u_error_bound(ft, plan) * math.sqrt(NU.v_stat(ft))
    assert abs(Um - u_mon(ft, th, cls, 5)) < 1e-12                                       # the fixed-order version agrees with the bincount version
    for bad in (lambda: NU.class_plan(np.array([0, 3]), 3), lambda: HP.reference_U(ft, th[:, :3], 64)):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("invalid input accepted")
    F = Fraction
    a, b = (F(3), F(4), F(1, 10)), (F(5), F(0), F(1, 10))
    assert HP.compare_moduli(a, b) == 0                                                  # |a| = 5 = |b|: the enclosures overlap
    assert HP.compare_moduli((F(3), F(4) + F(1, 2), F(1, 10)), b) == 1 and HP.compare_moduli(b, (F(3), F(4) + F(1, 2), F(1, 10))) == -1
    assert HP.compare_moduli((F(0), F(0), F(1, 10)), (F(0), F(0), F(1, 10))) == 0
    assert HP.compare_moduli((F(1), F(0), F(0)), (F(1), F(0), F(0))) == 1                # exact equality counts as >= (a tie counts toward the tail)
    lo, hi = HP.modulus_sq_bounds((F(3), F(-4), F(1, 10)))
    assert lo < 25 < hi and lo == F(29, 10) ** 2 + F(39, 10) ** 2 and hi == F(31, 10) ** 2 + F(41, 10) ** 2
    ms = list(range(1, 13))
    thof = lambda m: draw_theta(KEY, K, L, m)
    r = HP.classify(ft, th, thof, ms)
    exact = _dec_abs2(ft, th)
    for m in ms:
        assert r["status"][m] == ("ge" if _dec_abs2(ft, thof(m)) >= exact else "lt"), m
    assert r["precision"] == 128 and r["evaluations"] == 1 + len(ms)
    r = HP.classify(ft, th, lambda m: th if m == 3 else thof(m), [3, 4], p_start=3, p_cap=3)       # identical angles: exact tie; low precision: open
    assert r["status"] == {3: "ge", 4: "open"}


def _toy(K, L, seed):
    rng = np.random.default_rng(seed)
    return center_two_way(1.25 * (rng.random((K, L)) < 0.3) - 0.25)


def test_near_tie_resolution():
    """4, 7(D), 12.2: (a) every re-draw a near tie -> the second tier gives the tail count of independent decimal arithmetic, for U and for U_mon, for any
    number of workers and partition; (b) a precision cap that is too low leaves comparisons open -> 'uncertain', p_hat counts them; (c) a genuine near tie
    at the protocol's tolerance, on either side of 0.05, is decided as the exact arithmetic says; (d) a floating-point statistic that is wrong within the
    tolerance is corrected; (e) near ties that cannot change the decision are reported but not resolved."""
    key = seed_to_key(b"near-tie-test-" + bytes(range(64)))
    K, L, M = 3, 12, 99
    ft = _toy(K, L, 21)
    th_of = lambda m: draw_theta(key, K, L, m)
    ex_obs = _dec_abs2(ft, th_of(0))
    ex_tail = sum(_dec_abs2(ft, th_of(m)) >= ex_obs for m in range(1, M + 1))
    base = randomization_test(ft, key, M=M)
    assert base["tail"] == ex_tail and base["near_ties"] == 0 and base["resolution"] is None and base["error_bound"] < 1e-6
    big = randomization_test(ft, key, M=M, tol_u=1e9)                                    # (a) every re-draw is a near tie, p_lo <= 0.05 < p_hi
    assert big["near_ties"] == M and big["near_tie_draws"] == list(range(1, M + 1)) and big["p_hat_float"] == base["p_hat"]
    rs = big["resolution"]
    assert rs["n_ties"] == M and rs["n_open"] == 0 and rs["n_ge"] + rs["n_lt"] == M and rs["precision"] == 128 and rs["n_ge"] == ex_tail
    assert big["tail"] == ex_tail and big["p_hat"] == (1 + ex_tail) / (M + 1) == big["p_hat_lo"] == big["p_hat_hi"] and not big["uncertain"]
    assert big["decision"] == base["decision"] and big["reject"] == (base["decision"] == "reject")
    floating = randomization_test(ft, key, M=M, tol_u=1e9, resolve=False)               # the floating-point report alone is 'uncertain'
    assert floating["uncertain"] and floating["decision"] == "uncertain" and floating["resolution"] is None
    keys = ("p_hat", "tail", "p_hat_lo", "p_hat_hi", "uncertain", "decision", "reject", "near_tie_draws", "resolution")
    for kw in (dict(workers=2), dict(chunks=[range(60, 100), range(1, 60)]), dict(chunks=[[m] for m in range(M, 0, -1)])):
        r = randomization_test(ft, key, M=M, tol_u=1e9, **kw)
        assert all(r[k] == big[k] for k in keys), kw
    cap = randomization_test(ft, key, M=M, tol_u=1e9, p_start=2, p_cap=2)                # (b) precision cap too low
    assert cap["resolution"]["n_open"] > 0 and cap["uncertain"] and cap["decision"] == "uncertain" and not cap["reject"]
    assert cap["p_hat"] == cap["p_hat_hi"] > 0.05 >= cap["p_hat_lo"] and cap["tail"] == cap["resolution"]["n_ge"] + cap["resolution"]["n_open"]
    cls = np.random.default_rng(9).integers(0, 4, (K, L))                               # U_mon: the gate statistic (classes unrelated to the centring, so F_c != 0)
    ex_mon = _dec_abs2(ft, th_of(0), (cls, 4))
    assert max(abs(sum(Decimal(v) for v in ft[cls == c].tolist())) for c in range(4)) > 0.1
    mon_tail = sum(_dec_abs2(ft, th_of(m), (cls, 4)) >= ex_mon for m in range(1, M + 1))
    mon = randomization_test(ft, key, M=M, tol_u=1e9, classes=(cls, 4))
    assert mon["Lambda_obs"] is None and mon["tail"] == mon_tail and mon["resolution"]["n_open"] == 0 and mon["error_bound"] < 1e-6
    assert mon["p_hat_float"] == randomization_test(ft, key, M=M, classes=(cls, 4))["p_hat"] == (1 + mon_tail) / (M + 1)
    # (c) a near tie that straddles 0.05: observed angles = the 5th largest re-draw with one angle moved by 1e-8 (|U| up or down)
    mag = sorted(((abs(NU.u_stat(ft, th_of(m))), m) for m in range(1, M + 1)), reverse=True)
    m5 = mag[4][1]
    th5 = th_of(m5)
    exact_all = {m: _dec_abs2(ft, th_of(m)) for m in range(1, M + 1)}
    noise = {}
    sides = {}
    for i in range(th5.size):
        for sgn in (1.0, -1.0):
            tho = th5.copy(); tho.flat[i] += sgn * 1e-8
            above = _dec_abs2(ft, tho) > exact_all[m5]
            sides.setdefault(above, (tho, i))
        if len(sides) == 2:
            break
    assert len(sides) == 2
    wrong_float = 0
    for above, (tho, i) in sides.items():
        ex_t = sum(v >= _dec_abs2(ft, tho) for v in exact_all.values())
        assert ex_t == (4 if above else 5)
        r = randomization_test(ft, key, M=M, theta_obs=tho)
        assert r["near_ties"] >= 1 and m5 in r["near_tie_draws"] and r["resolution"]["n_ties"] >= 1 and r["resolution"]["n_open"] == 0
        assert r["tail"] == ex_t and r["decision"] == ("reject" if above else "not_reject") and r["p_hat_lo"] == r["p_hat_hi"] == r["p_hat"]
        # (d) the same instance with a floating-point statistic that is wrong by up to 0.9 tol sqrt(V): the second tier corrects it
        tol, sv, orig = 0.02, math.sqrt(NU.v_stat(ft)), RT._stat
        for salt in (b"a", b"b", b"c"):
            def noisy(f, t, plan, salt=salt):
                u = orig(f, t, plan)
                x = int.from_bytes(hashlib.sha256(salt + np.ascontiguousarray(t).tobytes()).digest()[:8], "big") / 2.0 ** 64
                return u + (2 * x - 1) * 0.9 * tol * sv * u / abs(u)
            RT._stat = noisy
            try:
                rn = randomization_test(ft, key, M=M, theta_obs=tho, tol_u=tol)
            finally:
                RT._stat = orig
            assert rn["tail"] == ex_t and rn["decision"] == r["decision"] and rn["resolution"] is not None and rn["resolution"]["n_open"] == 0, (above, salt)
            wrong_float += rn["p_hat_float"] != r["p_hat"]
    assert wrong_float >= 1, "the corrupted floating-point statistic never changed the floating-point p_hat: the test is not informative"
    # (e) a near tie far from 0.05: reported, not resolved, decision is the common side of p_hat_lo and p_hat_hi
    m30 = mag[29][1]
    tho = th_of(m30).copy(); tho.flat[0] += 1e-8
    r = randomization_test(ft, key, M=M, theta_obs=tho)
    assert r["near_ties"] >= 1 and r["resolution"] is None and r["decision"] == "not_reject" and r["p_hat_lo"] <= r["p_hat"] <= r["p_hat_hi"]
    print("   corrupted-float p_hat differed %d times" % wrong_float, end=" ")


def test_reject_implies_reference_reject():
    """4 (validity including the numerical tier): on many small instances, with tolerances from tiny to huge (so that near ties are rare, common or
    universal) and with precision schedules that settle everything or leave comparisons open, a decision 'reject' implies that the p_hat of the exact
    (decimal) arithmetic is <= 0.05 and 'not_reject' that it is > 0.05; p_hat_lo <= exact p_hat <= p_hat_hi always; the reported tail never falls below
    the exact tail after the second tier.  Hence P(reject) <= P(p_hat_exact <= 0.05) <= 0.05."""
    rnd = _pyrandom.Random(4242)
    K, L, M = 3, 8, 39
    seen = dict(reject=0, not_reject=0, uncertain=0, second_tier=0, open=0)
    for case in range(60):
        ft = _toy(K, L, 100 + case)
        key = seed_to_key(b"reject-implies-" + bytes([case]) + bytes(range(64)))
        th_of = lambda m: draw_theta(key, K, L, m)
        ex_obs = _dec_abs2(ft, th_of(0))
        ex_cnt = sum(_dec_abs2(ft, th_of(m)) >= ex_obs for m in range(1, M + 1))
        ex_p = (1 + ex_cnt) / (M + 1)
        # an observed assignment that is a large re-draw with one angle nudged makes near ties across 0.05 likely; case % 3 == 0 uses the m = 0 assignment as it is
        theta_obs = None
        if case % 3:
            mag = sorted(((abs(NU.u_stat(ft, th_of(m))), m) for m in range(1, M + 1)), reverse=True)
            theta_obs = th_of(mag[rnd.randrange(0, 3)][1]).copy()               # one of the three largest re-draws: p_hat is at the 0.05 edge (tail 1 or 2 of 39)
            theta_obs.flat[rnd.randrange(theta_obs.size)] += rnd.choice([-1.0, 1.0]) * rnd.choice([1e-9, 1e-6, 1e-3])
            ex_obs = _dec_abs2(ft, theta_obs)
            ex_cnt = sum(_dec_abs2(ft, th_of(m)) >= ex_obs for m in range(1, M + 1))
            ex_p = (1 + ex_cnt) / (M + 1)
        tol = rnd.choice([1e-6, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1])
        sched = rnd.choice([dict(), dict(p_start=2, p_cap=2), dict(p_start=4, p_cap=16)])
        r = randomization_test(ft, key, M=M, theta_obs=theta_obs, tol_u=tol, **sched)
        assert r["p_hat_lo"] <= ex_p <= r["p_hat_hi"], (case, tol, sched, r["p_hat_lo"], ex_p, r["p_hat_hi"])
        if r["decision"] == "reject":
            assert ex_p <= 0.05, (case, tol, sched)
        if r["decision"] == "not_reject":
            assert ex_p > 0.05, (case, tol, sched)
        if r["resolution"] is not None:
            seen["second_tier"] += 1
            seen["open"] += r["resolution"]["n_open"] > 0
            assert r["p_hat"] >= ex_p, (case, tol, sched)
        seen[r["decision"]] += 1
    assert seen["reject"] > 0 and seen["not_reject"] > 0 and seen["second_tier"] > 5 and seen["open"] > 0, seen
    print("   %s" % seen, end=" ")


def test_error_bound_and_library_qualification():
    """7(D): the error bound reproduces the numbers of the text at the registered scale, is computed from the data at hand, makes the analysis stop when it
    exceeds the accepted tolerance; the cos/sin library of the environment is qualified on a fixed sample and a degraded library is rejected."""
    N, V = 1.2e10, 1.5e9
    b = NU.error_bound_formula(1.2e6, 1e4, 6 * N, V)
    assert 1.2e-8 < b < 1.3e-8 and b < 1e-6
    lib_part = math.sqrt(2) * (NU.LIBM_ABS_TOL + 2.0 ** -53) * 6 * N / math.sqrt(V)
    assert 2.5e-9 < lib_part < 2.7e-9
    ft = _toy(4, 50, 3)
    assert NU.u_error_bound(ft) == NU.error_bound_formula(4, 50, float(np.abs(ft).sum()), NU.v_stat(ft))
    cls = (np.arange(200) % 7).reshape(4, 50)
    assert NU.u_error_bound(ft, NU.class_plan(cls, 7)) == NU.error_bound_formula(4, 50, float(np.abs(ft).sum()), NU.v_stat(ft), (29, 7))
    assert len(NU._libm_sample()) >= 1000
    assert 0.0 <= NU.libm_conformance() <= NU.LIBM_ABS_TOL
    real_cos = np.cos
    try:
        np.cos = lambda x: real_cos(x) + 1e-12
        NU.libm_conformance.cache_clear()
        try:
            NU.libm_conformance()
        except RuntimeError:
            pass
        else:
            raise AssertionError("a degraded cos was accepted")
    finally:
        np.cos = real_cos
        NU.libm_conformance.cache_clear()
    assert NU.libm_conformance() <= NU.LIBM_ABS_TOL
    saved = NU.LIBM_ABS_TOL
    try:
        NU.LIBM_ABS_TOL = 1e-3
        try:
            randomization_test(ft, KEY, M=9)
        except ValueError:
            pass
        else:
            raise AssertionError("a bound above the accepted tolerance did not stop the analysis")
    finally:
        NU.LIBM_ABS_TOL = saved


def test_decision_comparisons_certified():
    """10, 5.2: UCB_95 < r_min and the reproduction criterion (ii) are decided by the floating-point value when it is far from the constant and by the
    enclosure of the reference value otherwise; undecided comparisons are reported as such.  The stored FC-Rice table is non-decreasing (which the
    exclusion test relies on)."""
    xs, ucb = load_table(TABLE)
    assert (np.diff(ucb) >= 0).all() and float(u_of_x(18.0, xs, ucb)) < float(u_of_x(18.0 + 1e-9, xs, ucb))
    K, L = 3, 12
    ft = _toy(K, L, 4)
    th = draw_theta(KEY, K, L, 0)
    xs_id = np.linspace(0.0, 100.0, 100001); ucb_id = xs_id.copy()                       # u(x) = x
    sv = math.sqrt(NU.v_stat(ft))
    xh = abs(NU.u_stat(ft, th)) / sv
    ex = float(DE.exclusion(ft, th, xs_id, ucb_id, r_min=2 * xh / sv)["x_hat"])
    assert ex == xh
    far1, far2 = DE.exclusion(ft, th, xs_id, ucb_id, r_min=2 * xh / sv), DE.exclusion(ft, th, xs_id, ucb_id, r_min=0.5 * xh / sv)
    assert (far1["decision"], far1["reference_used"], far2["decision"], far2["reference_used"]) == ("exclude", False, "not_exclude", False)
    near1, near2 = DE.exclusion(ft, th, xs_id, ucb_id, r_min=xh * (1 + 1e-9) / sv), DE.exclusion(ft, th, xs_id, ucb_id, r_min=xh * (1 - 1e-9) / sv)
    assert (near1["decision"], near1["reference_used"], near1["precision"]) == ("exclude", True, 128)
    assert (near2["decision"], near2["reference_used"]) == ("not_exclude", True)
    tie = DE.exclusion(ft, th, xs_id, ucb_id, r_min=xh * (1 + 3e-13) / sv)
    assert tie["decision"] == "uncertain" and tie["reference_used"] and tie["precision"] == 1024
    real = DE.exclusion(ft, th, xs, ucb)                                                   # the stored table
    assert real["decision"] in ("exclude", "not_exclude") and real["ucb95"] == float(u_of_x(xh, xs, ucb)) / sv
    # reproduction criterion (ii)
    tho = draw_theta(KEY, K, L, 0)
    sig = center_two_way(np.cos(tho))                                                       # a signal aligned with theta: Lambda is huge
    assert DE.reproduction_ii(sig, tho, sig, tho)["decision"] == "pass" and not DE.reproduction_ii(sig, tho, sig, tho)["reference_used"]
    assert DE.reproduction_ii(sig, tho, sig, tho + math.pi)["decision"] == "fail"
    ft_r, th_r = _toy(K, L, 5), draw_theta(KEY, K, L, 1)
    r0 = DE.reproduction_ii(ft, th, ft_r, th_r)
    D = r0["D"]
    saved = (DE.CHI2_FLOAT, DE.CHI2_LO, DE.CHI2_HI)
    try:
        for shift, want in ((1e-9, "pass"), (-1e-9, "fail")):
            c = D + shift
            DE.CHI2_FLOAT, DE.CHI2_LO, DE.CHI2_HI = c, Fraction(c), Fraction(c)
            r = DE.reproduction_ii(ft, th, ft_r, th_r)
            assert r["decision"] == want and r["reference_used"] and r["precision"] == 128, (shift, r)
        Vo, Vr = NU.v_stat(ft), NU.v_stat(ft_r)                                             # a constant inside the enclosure of D: undecided
        re_o, im_o, e_o = HP.reference_U(ft, th, 1024)
        re_r, im_r, e_r = HP.reference_U(ft_r, th_r, 1024)
        zr, zi = re_r * Fraction(Vo) - re_o * Fraction(Vr), im_r * Fraction(Vo) - im_o * Fraction(Vr)
        c = (zr * zr + zi * zi) / (Fraction(Vo) * Fraction(Vr) * (Fraction(Vo) + Fraction(Vr)))
        DE.CHI2_FLOAT, DE.CHI2_LO, DE.CHI2_HI = float(c), c, c
        r = DE.reproduction_ii(ft, th, ft_r, th_r, p_cap=256)
        assert r["decision"] == "uncertain" and r["reference_used"] and r["precision"] == 256
    finally:
        DE.CHI2_FLOAT, DE.CHI2_LO, DE.CHI2_HI = saved
    assert DE.CHI2_LO <= Fraction(2) * Fraction(Decimal(20).ln()) <= DE.CHI2_HI or abs(DE.CHI2_FLOAT - 2 * math.log(20)) < 1e-12
    assert float(DE.CHI2_LO) == float(DE.CHI2_HI) == DE.CHI2_FLOAT == 2 * math.log(20) or abs(DE.CHI2_FLOAT - 2 * math.log(20)) < 1e-14


def test_ftilde_digest_and_seal():
    """2, 12.5.1 step 2: H_f~ is the SHA-256 of the little-endian binary64 rows in ascending absolute block number; independent of the order of the rows
    as listed, sensitive to one bit, and defined only for finite values."""
    ft = _toy(4, 9, 6)
    blocks = [10, 3, 7, 4]
    h = hashlib.sha256()
    for b in (3, 4, 7, 10):
        h.update(b"".join(__import__("struct").pack("<d", x) for x in ft[blocks.index(b)].tolist()))
    assert ftilde_digest(ft, blocks) == h.hexdigest()
    perm = [2, 0, 3, 1]
    assert ftilde_digest(ft[perm], [blocks[i] for i in perm]) == h.hexdigest()
    assert ftilde_digest(ft) == ftilde_digest(ft, [0, 1, 2, 3]) != h.hexdigest()
    ft2 = ft.copy(); ft2[1, 4] = np.nextafter(ft2[1, 4], 1e9)
    assert ftilde_digest(ft2, blocks) != h.hexdigest()
    assert ftilde_digest(np.asfortranarray(ft), blocks) == h.hexdigest()                       # memory layout does not matter
    for bad in (lambda: ftilde_digest(np.where(np.arange(36).reshape(4, 9) == 5, np.nan, ft)), lambda: ftilde_digest(ft, [1, 1, 2, 3]),
                lambda: ftilde_digest(ft, [1, 2, 3]), lambda: ftilde_digest(ft[0])):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("invalid input accepted")
    assert np.array_equal(center_two_way(_toy(5, 8, 1) + 0.0), center_two_way(_toy(5, 8, 1) + 0.0))           # deterministic on one machine


def test_reference_vectors_file_is_current():
    """tests/vectors/reference_kat.json is what scripts/make_kat.py generates, and each entry agrees with the independent decimal-arithmetic values (pi, cos,
    sin, enclosures of U and U_mon) and with ftilde_digest."""
    gen = runpy.run_path(os.path.join(HERE, "..", "scripts", "make_kat.py"))["generate_reference"]()
    assert open(REF_VEC).read() == json.dumps(gen, indent=1, sort_keys=True) + "\n", "run: python scripts/make_kat.py"
    for e in gen["pi_fixed"]:
        assert abs(Decimal(e["value"]) - _dec_pi(400) * Decimal(2) ** e["bits"]) <= 2
    for e in gen["sincos"]:
        t = float.fromhex(e["theta"])
        dc, ds = _dec_sincos(t, 330 if e["P"] > 200 else 90)
        assert abs(Decimal(int(e["cos"])) - dc * Decimal(2) ** e["P"]) <= 1 and abs(Decimal(int(e["sin"])) - ds * Decimal(2) ** e["P"]) <= 1
    for e in gen["reference_U"]:
        K, L = e["shape"]
        w, th = np.array([float.fromhex(x) for x in e["w"]]).reshape(K, L), np.array([float.fromhex(x) for x in e["theta"]]).reshape(K, L)
        cl = None if e["classes"] is None else (np.array(e["classes"]).reshape(K, L), 3)
        dre, dim = _dec_U(w, th, cl)
        fr = lambda s: Fraction(int(s.split("/")[0]), int(s.split("/")[1]))
        err = fr(e["err"])
        assert abs(Decimal(fr(e["re"]).numerator) / Decimal(fr(e["re"]).denominator) - dre) <= Decimal(err.numerator) / Decimal(err.denominator)
        assert abs(Decimal(fr(e["im"]).numerator) / Decimal(fr(e["im"]).denominator) - dim) <= Decimal(err.numerator) / Decimal(err.denominator)
    for e in gen["ftilde_digest"]:
        K, L = e["shape"]
        assert ftilde_digest(np.array([float.fromhex(x) for x in e["ft"]]).reshape(K, L), e["blocks"]) == e["sha256"]


def test_theta_rounding_bound():
    """3.1 point 7: the binary64 angle differs from alpha + 2 pi pi_i / L (alpha the block's binary64 rotation, pi exact) by at most design.theta_rounding_bound()
    (about 3.04e-15 rad, derived from the rounding model).  Checked in exact arithmetic (pi to 200 bits) on angles of real blocks and on extreme inputs (rotation
    near 2 pi, permutation entries near L - 1, where the quotient and the sum are largest); so the statistic built on the binary64 angles differs from the one
    on ideal angles by at most that bound times sum|f~|."""
    pi2 = Fraction(2 * HP.pi_fixed(200), 1 << 200)
    bound = D.theta_rounding_bound()
    assert 3.0e-15 < bound < 3.1e-15, bound
    worst = Fraction(0)
    for (L, b, m) in ((3, 0, 0), (10, 1, 1), (100, 5, 3), (1000, 2 ** 20, 7), (10000, 3, 0), (10000, 9, 4)):
        th, alpha, perm = block_angles(make_block_rng(KEY, b, m), L)
        a = Fraction(alpha)
        for t, p in zip(th.tolist(), perm.tolist()):
            worst = max(worst, abs(Fraction(t) - (a + pi2 * p / L)))
    rnd = _pyrandom.Random(5)
    two_pi = P.TWO_PI                                                           # the expression of philox_ref.block_theta, on chosen inputs
    for _ in range(60000):
        L = rnd.choice([10000, 9999, 8192, 5000, 1000, 100, 33])
        alpha = two_pi * ((1 << 53) - 1 - rnd.getrandbits(rnd.choice([1, 3, 8, 20, 53]))) * 2.0 ** -53
        p = L - 1 - rnd.randrange(min(L, 40)) if rnd.random() < 0.7 else rnd.randrange(L)
        t = alpha + (two_pi * p) / L
        worst = max(worst, abs(Fraction(t) - (Fraction(alpha) + pi2 * p / L)))
    assert float(worst) <= bound, (float(worst), bound)
    print("   max |theta - ideal| = %.2e rad (bound %.3e)" % (float(worst), bound), end=" ")


def _rayleigh_top(rng, n, m, trials):
    """Largest m of n iid Rayleigh values per trial, in descending order (shape (trials, m)).  Uniform order statistics from the top, by Renyi's
    recursion u_(n) = U1^(1/n), u_(n-j) = u_(n-j+1) U_(j+1)^(1/(n-j)) (given the maximum the others are iid uniform below it); x = sqrt(-2 ln(1-u))."""
    lu = np.cumsum(np.log(rng.random((trials, m))) / (n - np.arange(m)), axis=1)
    return np.sqrt(-2.0 * np.log1p(-np.exp(lu)))


def test_near_tie_rate_formulas():
    """Protocol 7(D) (trigger probability, mean number of near ties).  Under H0 the n = M+1 values |U|/sqrt(V) (observed and re-draws) are exchangeable and, in
    the Rayleigh idealisation, iid with density f(x) = x exp(-x^2/2).  (a) Exact: the second tier is triggered iff the k-th largest re-draw (k = 0.05 n) lies
    in [x_obs - d, x_obs + d); checked against the definition with p_hat_lo, p_hat_hi.  (b) The k-th and (k+1)-th largest of n values are
    ~ Exp(n f(x_q)) apart, x_q = sqrt(2 ln 20), so P(trigger) = (2/n)(1 - exp(-n f(x_q) d)) to first order: simulated for n = 200, 2000 (the generator of
    the largest values is checked against full sorting).  (c) E[#near ties] = M P(window) = M 2 d int f^2 = M 2 d sqrt(pi)/4 for small d."""
    from scipy import integrate, stats as _st
    rng = np.random.default_rng(20261014)
    # (a) the trigger event equals the order-statistic condition, on the definition of p_hat_lo / p_hat_hi
    for n in (20, 60, 200):
        k = int(round(0.05 * n))
        assert abs(k - 0.05 * n) < 1e-9
        trig = 0
        for _ in range(4000):
            d = float(rng.choice([0.0, 1e-3, 0.02, 0.2]))
            v = np.sqrt(-2.0 * np.log(rng.random(n)))
            i = int(rng.integers(n))
            obs, others = v[i], np.delete(v, i)
            p_lo = (1 + np.count_nonzero(others >= obs + d)) / n
            p_hi = (1 + np.count_nonzero(others >= obs - d)) / n
            hit = (p_lo <= 0.05) != (p_hi <= 0.05)
            y = np.sort(others)[::-1][k - 1]                       # k-th largest re-draw
            assert hit == (obs - d <= y < obs + d), (n, d, hit, y, obs)
            trig += hit
        assert trig > 0
    # the top-order-statistic generator agrees with full sorting (two-sample KS on z_(k) and on the gap z_(k) - z_(k+1))
    n, k = 200, 10
    top = _rayleigh_top(rng, n, k + 1, 20000)
    full = np.sort(np.sqrt(-2.0 * np.log(rng.random((20000, n)))), axis=1)[:, ::-1][:, :k + 1]
    for a, b in ((top[:, k - 1], full[:, k - 1]), (top[:, k - 1] - top[:, k], full[:, k - 1] - full[:, k])):
        assert _st.ks_2samp(a, b).pvalue > 1e-3
    # (b) the gap law at the registered n and at n = 200; first-order trigger probability against simulation
    f = lambda x: x * np.exp(-x * x / 2.0)
    xq = np.sqrt(2.0 * np.log(20.0))
    assert abs(f(xq) - xq / 20.0) < 1e-15
    for n, k, trials in ((2000, 100, 40000), (200, 10, 200000)):
        lam = n * f(xq)
        z = _rayleigh_top(rng, n, k + 1, trials)
        g = z[:, k - 1] - z[:, k]
        for lam_d in (0.1, 1.0):
            d = lam_d / lam
            got = np.mean(g < d)
            want = 1.0 - np.exp(-lam_d)
            assert abs(got - want) < 5.0 * np.sqrt(want * (1 - want) / trials) + 0.01 * want, (n, lam_d, got, want)
    n, k, lam_d, trials = 200, 10, 0.01, 8000000
    d = lam_d / (n * f(xq))
    hits = 0
    for _ in range(8):
        z = _rayleigh_top(rng, n, k + 4, trials // 8)
        for r in range(1, k + 5):                                # rank r of the observation, uniform on 1..n, independent of the sorted values
            y = z[:, k] if r <= k else z[:, k - 1]
            hits += np.count_nonzero((z[:, r - 1] - d <= y) & (y < z[:, r - 1] + d))
    got = hits / n / trials                                       # each rank has probability 1/n; ranks above k+4 need five points inside the window
    want = 2.0 / n * (1.0 - np.exp(-lam_d))
    assert abs(got / want - 1.0) < 0.15, (got, want)
    reg = 2.0 / 2000 * (1.0 - np.exp(-2000 * f(xq) * 2e-6))
    assert 4.8e-7 < reg < 5.0e-7, reg                                # the figure quoted in the protocol (4.9e-7)
    # (c) the mean number of near ties
    Fx = lambda x: 1.0 - np.exp(-np.maximum(x, 0.0) ** 2 / 2.0)
    def window(d):
        return integrate.quad(lambda y: f(y) * (Fx(y + d) - Fx(y - d)), 0.0, 12.0, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    for d in (2e-6, 1e-2):
        assert abs(window(d) / (2.0 * d * np.sqrt(np.pi) / 4.0) - 1.0) < (1e-5 if d < 1e-5 else 1e-3)
    d, N = 1e-2, 10000000
    x, y = np.sqrt(-2.0 * np.log(rng.random(N))), np.sqrt(-2.0 * np.log(rng.random(N)))
    got = np.mean((x - d <= y) & (y < x + d))
    assert abs(got / window(d) - 1.0) < 0.02, (got, window(d))
    reg_mean = 1999 * 2.0 * 2e-6 * np.sqrt(np.pi) / 4.0
    assert 3.4e-3 < reg_mean < 3.6e-3, reg_mean                      # the figure quoted in the protocol (3.5e-3)
    print("   P(trigger) = %.2e, mean near ties = %.2e at the registered scale" % (reg, reg_mean), end=" ")


def test_trace_checks_detect_a_wrong_number():
    """check_trace.constants_checks: every constant and every measured number that protocol 4, 7(D) and 10 quote is compared with the code or with
    results/numerics_check.json; the repository is consistent, and altering a quoted number is reported."""
    import snd
    ct = runpy.run_path(os.path.join(HERE, "..", "scripts", "check_trace.py"))
    ptxt = open(os.path.join(HERE, "..", "..", "protocol", "SND_protocol_V%s.md" % snd.__version__), encoding="utf8").read()
    base = ct["constants_checks"](ptxt)
    assert all(ok for ok, _ in base), [d for ok, d in base if not ok]
    muts = [("4.9\\times10^{-7}", "4.8\\times10^{-7}"), ("3.5\\times10^{-3}", "3.6\\times10^{-3}"), ("3088 個角度", "3000 個角度"),
            ("7.4\\times10^{-16}", "7.5\\times10^{-16}"), ("每項約 10–16 微秒", "每項約 30–31 微秒"), ("2\\ln20=5.99146", "2\\ln20=5.99147"),
            ("1.0\\times10^{-8}", "1.1\\times10^{-8}"), ("3.04\\times10^{-15}", "2.04\\times10^{-15}"), ("1024 位元", "2048 位元"), ("\\delta=2\\times10^{-6}\\sqrt V", "\\delta=3\\times10^{-6}\\sqrt V"),
            ("`numerics.u_mon_stat`", "`numerics.u_mon_total`"), ("約 1.4–2.3 核", "約 4.5–4.7 核")]
    for a, b in muts:
        assert a in ptxt, a
        bad = ct["constants_checks"](ptxt.replace(a, b))
        assert not all(ok for ok, _ in bad), a


def test_c_implementation_matches_python():
    """crosscheck/philox_theta.c is a second implementation of protocol 3.1 written from the protocol text (C, unsigned __int128, no numpy, no Python code
    shared).  Its output is bit-identical to the known-answer vectors (blocks, theta digests via the same serialisation) and to the production path
    on random keys, blocks, re-draws and sizes.  Also: make_block_rng gives the same angles for a tuple, a list and a uint64 array key."""
    gcc = shutil.which("gcc")
    if gcc is None:
        print("   SKIPPED: no C compiler on this machine (the check is repeated by whoever has one)", end=" ")
        return
    import tempfile
    src = os.path.join(HERE, "..", "crosscheck", "philox_theta.c")
    with tempfile.TemporaryDirectory() as tmp:
        exe = os.path.join(tmp, "philox_theta")
        subprocess.run([gcc, "-O1", "-ffp-contract=off", "-o", exe, src], check=True)
        def run(key, b, m, L, raw=False):
            args = [exe, "%x" % key[0], "%x" % key[1], str(b), str(m), str(L)] + (["raw"] if raw else [])
            return subprocess.run(args, capture_output=True, check=True).stdout
        gen = json.load(open(os.path.join(HERE, "vectors", "theta_kat.json")))
        w = lambda x: int(x, 16)
        for e in gen["block_theta"]:
            lines = run([w(k) for k in e["key"]], e["block"], e["m"], e["L"]).decode().split("\n")
            rows = [l.split() for l in lines[1:] if l]
            assert int(lines[0], 16) == w(e["alpha_bits"]) and [int(r[0]) for r in rows] == e["perm"]
            assert [int(r[1], 16) for r in rows] == [w(x) for x in e["theta_bits"]]
        for e in gen["theta_digest"]:
            h = hashlib.sha256()
            for b in sorted(e["blocks"]):
                h.update(run([w(k) for k in e["key"]], b, e["m"], e["L"], raw=True))
            assert h.hexdigest() == e["sha256"]
        rnd = _pyrandom.Random(20260101)
        for _ in range(24):
            key = (rnd.getrandbits(64), rnd.getrandbits(64))
            b = rnd.choice([rnd.randrange(0, 10 ** 6), rnd.randrange(0, 2 ** 40), rnd.randrange(0, 2 ** 64)])
            m, L = rnd.randrange(0, 2000), rnd.choice([1, 2, 3, 7, 64, 65, 1000, 10000])
            mine = np.frombuffer(run(key, b, m, L, raw=True), dtype="<u8")
            for k in (key, list(key), np.array(key, dtype=np.uint64)):
                th = block_angles(make_block_rng(k, b, m), L)[0]
                assert np.array_equal(mine, th.astype("<f8").view("<u8")), (key, b, m, L, type(k))
    print("   KAT blocks %d, digests %d, random cases 24" % (len(gen["block_theta"]), len(gen["theta_digest"])), end=" ")


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
