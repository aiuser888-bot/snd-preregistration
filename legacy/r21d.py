# r21d (V2.0.08): resolution of U_mon:  s_mon/sqrt(V) ~ eta,  eta^2 = sum_c F_c^2/n_c / sum f~^2  (class share of f~ energy).
import sys, numpy as np
sys.path.insert(0, "code")
from snd.leak import u_mon
from snd.stats import center_two_way
L, K, M = 100, 400, 1000
j = np.arange(L); base = 2*np.pi*j/L; g = (j/L-0.5)*0.4
rng = np.random.default_rng([20261030, 5])
for trial in range(4):
    c = rng.choice([-1.0, 1.0], (K, 1)); d = c*g[None, :]
    f = 1.25*(rng.random((K, L)) < np.clip(0.2+d, .01, .99))-.25; ft = center_two_way(f); V = L/(2*(L-1))*(ft**2).sum()
    for nm, cls in (("4 classes(sign x posHalf)", (2*(c > 0).astype(int)+np.broadcast_to((g > 0).astype(int), (K, L))).astype(int)),
                    ("16 classes(sign x pos quarter)", (4*(c > 0).astype(int)+np.broadcast_to((4*j//L), (K, L))).astype(int))):
        nc = cls.max()+1; n = np.bincount(cls.ravel(), minlength=nc); F = np.bincount(cls.ravel(), weights=ft.ravel(), minlength=nc)
        eta2 = (F**2/n).sum()/(ft**2).sum()
        Um = np.empty(M, complex)
        for m in range(M):
            a = rng.uniform(0, 2*np.pi, (K, 1)); th = np.take_along_axis(a+base, np.argsort(rng.random((K, L)), axis=1), axis=1)
            Um[m] = u_mon(ft, th, cls, nc)
        s2 = (np.abs(Um)**2).mean()/2
        print("%-32s eta_pred=%.4f   s_mon/sqrt(V)=%.4f  (M=%d)" % (nm, np.sqrt(eta2), np.sqrt(s2/V), M))
