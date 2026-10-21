# r21 (V2.0.08) Part 1: Lemma 3.  Observed type-I vs alpha(delta_r) with delta_r = |sum f~ q_true|/sqrt(V), per replication.
# Also Cauchy-Schwarz bound delta <= beta_max(q_rms) and alignment omega.
import sys, numpy as np
sys.path.insert(0, "code")
from snd.leak import type1_of_beta, beta_of_B, CHI2_2_95
from snd.stats import center_two_way
from snd.accept import wilson
L, K, R = 100, 400, int(sys.argv[1]) if len(sys.argv) > 1 else 2000
j = np.arange(L); base = 2*np.pi*j/L; g = (j/L-0.5)*0.4
def run(cellname, dm, lk, q, center, ci, R=R):
    rng = np.random.default_rng([20261030, 1, ci]); rej = 0; pred = 0.0; om = []; dl = []; viol = 0
    for _ in range(R):
        c = rng.choice([-1.0, 1.0], (K, 1))
        d = np.broadcast_to(g, (K, L)) if dm == "common" else c*g[None, :]
        p = np.clip(0.2+d, .01, .99); f = 1.25*(rng.random((K, L)) < p)-.25
        a = rng.uniform(0, 2*np.pi, (K, 1)); th = np.take_along_axis(a+base, np.argsort(rng.random((K, L)), axis=1), axis=1)
        qm = np.broadcast_to(np.where(g > 0, q, 0.0), (K, L)) if lk == "static" else np.where(d > 0, q, 0.0)
        th = np.where(rng.random((K, L)) < qm, 0.0, th)
        ft = center_two_way(f) if center == "two" else f-f.mean(axis=1, keepdims=True)
        U = (ft*np.exp(1j*th)).sum(); V = L/(2*(L-1))*(ft**2).sum()
        rej += abs(U)**2/V > CHI2_2_95
        dd, w, qr = beta_of_B(ft, qm.astype(complex)); pred += float(type1_of_beta(dd)); om.append(w); dl.append(dd)
        viol += dd > np.sqrt(2*(L-1)*K*L/L)*qr*(1+1e-9)
    lo, hi = wilson(rej, R)
    print("%-26s %-3s q=%.3f  delta=%.2f omega=%.3f | obs %.3f [%.3f,%.3f]  pred alpha(delta) %.3f  CSviol=%d" %
          (cellname, center, q, np.mean(dl), np.mean(om), rej/R, lo, hi, pred/R, viol), flush=True)
ci = 0
for q in (0.0, 0.01, 0.02, 0.03, 0.05):
    for name, dm, lk, cen in (("B static x common", "common", "static", "two"), ("B static x common", "common", "static", "row"),
                              ("C blockvar x blockdrift", "block", "blockvar", "two")):
        run(name, dm, lk, q, cen, ci); ci += 1
