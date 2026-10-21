# r21b (V2.0.08) Part 2: worst-phase coverage under bias (full pipeline).  Part 3: U_mon randomization test (size, power, UB coverage).
import sys, numpy as np
sys.path.insert(0, "code")
from scipy import stats
from snd.leak import type1_of_beta, beta_of_B, coverage_worst_phase, u_mon, CHI2_2_95
from snd.stats import center_two_way
from snd.fcrice import load_table
from snd.accept import wilson
xs, ucb = load_table("results/fc_rice_table.npz")
L, K = 100, 400
j = np.arange(L); base = 2*np.pi*j/L; g = (j/L-0.5)*0.4
part = sys.argv[1]
def gen(rng, dm, lk, q, t=0.0, phi=np.pi):
    c = rng.choice([-1.0, 1.0], (K, 1))
    d = np.broadcast_to(g, (K, L)) if dm == "common" else (c*g[None, :] if dm == "block" else np.zeros((K, L)))
    a = rng.uniform(0, 2*np.pi, (K, 1)); th = np.take_along_axis(a+base, np.argsort(rng.random((K, L)), axis=1), axis=1)
    qm = (np.broadcast_to(np.where(g > 0, q, 0.0), (K, L)) if lk == "static" else (np.where(d > 0, q, 0.0) if lk == "blockvar" else np.zeros((K, L))))
    th = np.where(rng.random((K, L)) < qm, 0.0, th)
    p = np.clip(0.2+d+t*np.cos(th-phi), .01, .99); f = 1.25*(rng.random((K, L)) < p)-.25
    return c, d, th, qm, f
if part == "cov":
    R = int(sys.argv[2]); ci = 0
    for q in (0.03, 0.05):
        for nu0 in (1.0, 2.0, 3.0, 4.0):
            rng = np.random.default_rng([20261030, 2, ci]); ci += 1; hit = 0; pred = 0.0; dl = []; pred0 = 0.0
            t = nu0/(5*np.sqrt(L/(2*(L-1))*K*L*0.25))   # r = 5 t,  nu = r*sqrt(V),  V ~ N*0.25*L/(2(L-1))
            for _ in range(R):
                c, d, th, qm, f = gen(rng, "block", "blockvar", q, t)
                ft = center_two_way(f); U = (ft*np.exp(1j*th)).sum(); V = L/(2*(L-1))*(ft**2).sum()
                r = 5*t; nu = r*np.sqrt(V); xh = abs(U)/np.sqrt(V)
                u = np.interp(xh, xs, ucb); hit += (u/np.sqrt(V) >= r)
                dd, w, qr = beta_of_B(ft, qm.astype(complex)); dl.append(dd)
                pred += coverage_worst_phase(nu, dd, xs, ucb); pred0 += coverage_worst_phase(nu, 0.0, xs, ucb)
            lo, hi = wilson(hit, R)
            print("COV C-class q=%.2f nu~%.0f delta=%.2f | obs %.3f [%.3f,%.3f]  pred(|nu-delta|) %.3f  nominal-no-bias %.3f" % (q, nu0, np.mean(dl), hit/R, lo, hi, pred/R, pred0/R), flush=True)
else:
    R = int(sys.argv[2]); M = 199; ci = int(sys.argv[3]) * 10
    def classes(c, mode, rng):
        pos = np.broadcast_to((g > 0).astype(int), (K, L))
        s = ((c[:, 0] > 0).astype(int))
        if mode == "perfect": sg = s
        elif mode == "noisy": sg = np.where(rng.random(K) < 0.25, 1-s, s)
        else: sg = rng.integers(0, 2, K)          # uninformative proxy
        return (2*sg[:, None]+pos).astype(int)
    cells = [("no leak, drift only", "block", "none", 0.0, 0.0, "perfect"), ("C q=.03 proxy perfect", "block", "blockvar", .03, 0.0, "perfect"),
             ("C q=.03 proxy noisy(25% flip)", "block", "blockvar", .03, 0.0, "noisy"), ("C q=.03 proxy uninformative", "block", "blockvar", .03, 0.0, "none"),
             ("H1 nu~3 no leak, proxy perfect", "block", "none", 0.0, 3/(5*np.sqrt(L/(2*(L-1))*K*L*0.25)), "perfect")]
    sel = [int(x) for x in sys.argv[4].split(",")]
    for k in sel:
        name, dm, lk, q, t, pm = cells[k]; rng = np.random.default_rng([20261030, 3, k, ci]); rej = 0; cov = 0; ests = []; true = []
        for _ in range(R):
            c, d, th, qm, f = gen(rng, dm, lk, q, t); ft = center_two_way(f); V = L/(2*(L-1))*(ft**2).sum()
            cls = classes(c, pm, rng); Uo = u_mon(ft, th, cls, 4)
            Um = np.empty(M, complex)
            for m in range(M):
                a = rng.uniform(0, 2*np.pi, (K, 1)); thm = np.take_along_axis(a+base, np.argsort(rng.random((K, L)), axis=1), axis=1)
                Um[m] = u_mon(ft, thm, cls, 4)
            cnt = (np.abs(Um) >= abs(Uo)).sum(); rej += (1+cnt)/(M+1) <= 0.05
            s = np.sqrt((np.abs(Um)**2).mean()/2); ub = abs(Uo)+np.sqrt(CHI2_2_95)*s
            n = np.bincount(cls.ravel(), minlength=4).astype(float); F = np.bincount(cls.ravel(), weights=ft.ravel(), minlength=4)
            qc = np.bincount(cls.ravel(), weights=qm.ravel(), minlength=4)/n; Bmon = abs((F*qc).sum())
            cov += Bmon <= ub; ests.append(abs(Uo)/np.sqrt(V)); true.append(Bmon/np.sqrt(V))
        lo, hi = wilson(rej, R); clo, chi = wilson(cov, R)
        print("MON %-34s | reject %.3f [%.3f,%.3f] | mean est delta_mon %.2f vs true delta_mon %.2f | UB cover %.3f [%.3f,%.3f]" % (name, rej/R, lo, hi, np.mean(ests), np.mean(true), cov/R, clo, chi), flush=True)
