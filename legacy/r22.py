# r22 (V2.0.09): effects of raising V_target 1.375e9 -> 1.5e9 (T = 1e-4*sqrt(V): 3.708 -> 3.873).
# Part A: gross/net exclusion probability and power (gate false alarm 0.05, Gaussian model).  Part B: reproduction success probabilities.
import sys, numpy as np
sys.path.insert(0, "code")
from scipy import stats
from scipy.optimize import brentq
from snd.fcrice import load_table
xs, ucb = load_table("results/fc_rice_table.npz")
c95 = stats.chi2.ppf(.95, 2)
def p_excl(T):                      # H0: P(UCB<r_min) = P(Rayleigh xhat < x_min(T))
    xm = xs[np.argmax(ucb >= T)]; return 1 - np.exp(-xm**2/2)
def power(T): return float(stats.ncx2.sf(c95, 2, T*T))
print("T      Vtarget   P(excl|H0) gross   power(r_min) gross")
for V in (1.375e9, 1.5e9):
    T = 1e-4*np.sqrt(V); print("%.4f  %.3e  %.4f             %.4f" % (T, V, p_excl(T), power(T)))
Vg = (brentq(lambda t: p_excl(t) - .80, 3.0, 4.5, xtol=1e-5)*1e4)**2
Vn = (brentq(lambda t: 0.95*p_excl(t) - .80, 3.0, 4.5, xtol=1e-5)*1e4)**2
print("critical V gross(0.80) = %.4e ; net(0.80, gate x0.95 at eta->0) = %.4e" % (Vg, Vn))
for V in (1.375e9, 1.5e9): print("V_target %.3e: margin vs gross %.1f%%, vs net %.1f%%" % (V, 100*(V/Vg-1), 100*(V/Vn-1)))
rng = np.random.default_rng([20261032, 1]); n = 3_000_000
za = rng.standard_normal((n, 2)); zb = rng.standard_normal((n, 2))
print("\nPart A: net values with gate (alpha_mon=0.05), per eta")
for V in (1.375e9, 1.5e9):
    T = 1e-4*np.sqrt(V); xm = xs[np.argmax(ucb >= T)]
    for eta in (0.0, 0.25, 0.5, 0.75):
        e = max(eta, 1e-9); um = eta*za; alarm = ((um**2).sum(1)/e**2 > c95) if eta > 0 else (rng.random(n) < 0.05)
        ex = np.sqrt((((um + np.sqrt(1-eta**2)*zb)**2).sum(1))) < xm
        tot = um + np.sqrt(1-eta**2)*zb + np.array([T, 0]); rej = (tot**2).sum(1) > c95
        print("  V=%.3e eta=%.2f  net P(excl|H0)=%.4f   net power(r_min)=%.4f" % (V, eta, (ex & ~alarm).mean(), (rej & ~alarm).mean()))
# Part B: reproduction (Gaussian large-K model; as r12 but nu0 = T)
print("\nPart B: reproduction success (i)&(ii), conditional on original significant; chi2_{2,.95}")
xc = np.sqrt(c95)
def run(nu0, mult, truth, rng):
    srep = 1/np.sqrt(mult); zo = rng.normal(size=(4_000_000, 2)) + [nu0, 0]
    zr = rng.normal(size=(4_000_000, 2))*srep + np.array(truth)
    so = np.hypot(*zo.T) >= xc; sr = np.hypot(*zr.T)/srep >= xc; d2 = ((zr-zo)**2).sum(axis=1)/(1+srep**2)
    return ((d2 <= c95) & sr & so).sum()/so.sum()
for V in (1.375e9, 1.5e9):
    nu0 = 1e-4*np.sqrt(V); rr = np.random.default_rng([20261032, 2, int(V/1e8)])
    cases = [("same effect", [nu0, 0]), ("amp 0.5x", [0.5*nu0, 0]), ("phase +60", [nu0*np.cos(np.pi/3), nu0*np.sin(np.pi/3)]), ("phase +90", [0, nu0])]
    print(" V_target=%.3e (nu0=%.3f)" % (V, nu0))
    for mult in (1.0, 2.0):
        print("   V_rep=%dx: " % mult + "  ".join("%s %.3f" % (nm, run(nu0, mult, t, rr)) for nm, t in cases))
    print("   c grid (power(i)*0.95): " + "  ".join("c=%g:%.3f" % (c, power(nu0*np.sqrt(c))*0.95) for c in (0.25, 0.5, 1, 2, 4)))
    print("   c grid (i) power: " + "  ".join("c=%g:%.3f" % (c, power(nu0*np.sqrt(c))) for c in (0.25, 0.5, 1, 2, 4)))
# Part C: duration
for V in (1.375e9, 1.5e9):
    for s2 in (0.25, 0.15, 0.10, 0.05):
        N = 2*V/s2; days = N/(4*100)/86400; print("V=%.3e sigma_f^2=%.2f  N=%.3e  months=%.1f" % (V, s2, N, days/30.4375))
    print("gate sigma_f^2 for N_max=3e10: %.4f  (months at N_max %.1f)" % (2*V/3e10, 3e10/(4*100)/86400/30.4375))
