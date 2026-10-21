# r21c (V2.0.08): gated discovery  P(main rejects AND U_mon does not reject) under monitored-channel C-class leakage.
import sys, numpy as np
sys.path.insert(0, "code")
from snd.leak import u_mon, CHI2_2_95
from snd.stats import center_two_way
from snd.accept import wilson
L, K, M = 100, 400, 199
j = np.arange(L); base = 2*np.pi*j/L; g = (j/L-0.5)*0.4
R = int(sys.argv[1]); mode = sys.argv[2]; qs = [float(x) for x in sys.argv[3].split(",")]
for ci, q in enumerate(qs):
    rng = np.random.default_rng([20261030, 4, ord(mode[0]), ci]); n_main = n_mon = n_both = 0
    for _ in range(R):
        c = rng.choice([-1.0, 1.0], (K, 1)); d = c*g[None, :]
        a = rng.uniform(0, 2*np.pi, (K, 1)); th = np.take_along_axis(a+base, np.argsort(rng.random((K, L)), axis=1), axis=1)
        qm = np.where(d > 0, q, 0.0); th = np.where(rng.random((K, L)) < qm, 0.0, th)
        f = 1.25*(rng.random((K, L)) < np.clip(0.2+d, .01, .99))-.25; ft = center_two_way(f); V = L/(2*(L-1))*(ft**2).sum()
        s = (c[:, 0] > 0).astype(int)
        sg = s if mode == "perfect" else np.where(rng.random(K) < 0.25, 1-s, s)
        cls = (2*sg[:, None]+np.broadcast_to((g > 0).astype(int), (K, L))).astype(int)
        Uo = (ft*np.exp(1j*th)).sum(); Mo = u_mon(ft, th, cls, 4); Um = np.empty(M, complex); Mm = np.empty(M, complex)
        for m in range(M):
            a = rng.uniform(0, 2*np.pi, (K, 1)); thm = np.take_along_axis(a+base, np.argsort(rng.random((K, L)), axis=1), axis=1)
            Um[m] = (ft*np.exp(1j*thm)).sum(); Mm[m] = u_mon(ft, thm, cls, 4)
        pm = (1+(np.abs(Um) >= abs(Uo)).sum())/(M+1); pn = (1+(np.abs(Mm) >= abs(Mo)).sum())/(M+1)
        sd = np.sqrt((np.abs(Mm)**2).mean()/2); ub = (abs(Mo)+np.sqrt(CHI2_2_95)*sd)/np.sqrt(V)
        main = pm <= .05; mon = pn <= .05; n_main += main; n_mon += mon; n_both += main and not mon
    f_ = lambda k: "%.3f [%.3f,%.3f]" % (k/R, *wilson(k, R))
    print("GATE proxy=%-7s q=%.3f | main rej %s | mon rej %s | main & ~mon %s" % (mode, q, f_(n_main), f_(n_mon), f_(n_both)), flush=True)
