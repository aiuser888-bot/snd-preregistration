"""Recompute the closed-form numbers quoted in protocol V2.0.10 (sections 6, 10, 12) and assert them. Exit 0 = all consistent."""
import os, sys
import numpy as np
from scipy import stats
from snd.fcrice import load_table
from snd.leak import (type1_of_beta, coverage_min_over_nu, coverage_worst_phase, q_rms_for_beta, hash_D_bound,
                      devices_D, beta_ub_from_D, size_ub_from_D, best_of_k_size, seed_entropy_ok, retry_size_ub)
root = os.path.join(os.path.dirname(__file__), "..", "..", "results")
xs, ucb = load_table(os.path.join(root, "fc_rice_table.npz"))
c95 = stats.chi2.ppf(.95, 2); bad = 0
def chk(name, val, ref, tol):
    global bad
    ok = abs(val - ref) <= tol; bad += (not ok)
    print("%-4s %-52s %.4f (protocol %.4f, tol %.4f)" % ("OK" if ok else "FAIL", name, val, ref, tol))
V = 1.5e9; T = 1e-4*np.sqrt(V)
chk("T = r_min sqrt(V_target)", T, 3.873, 5e-4)
xm = xs[np.argmax(ucb >= T)]
chk("P(UCB<r_min|H0) gross", 1-np.exp(-xm**2/2), 0.891, 5e-4)
chk("power at r_min gross", stats.ncx2.sf(c95, 2, T*T), 0.944, 5e-4)
chk("net P(excl) lower bound 0.95*gross", 0.95*(1-np.exp(-xm**2/2)), 0.846, 1e-3)
chk("net power 0.95*gross", 0.95*stats.ncx2.sf(c95, 2, T*T), 0.897, 1e-3)
chk("N = 2V/sigma_f^2 (1e10)", 2*V/0.25/1e10, 1.2, 1e-9)
months = lambda N: N/(4*100)/86400/30.4375
chk("months at sigma_f^2=0.25", months(2*V/0.25), 11.4, 0.05)
chk("months at sigma_f^2=0.15", months(2*V/0.15), 19.0, 0.05)
chk("months at sigma_f^2=0.10", months(2*V/0.10), 28.5, 0.05)
chk("months at sigma_f^2=0.05", months(2*V/0.05), 57.0, 0.1)
chk("gate sigma_f^2 for N_max=3e10", 2*V/3e10, 0.100, 1e-9)
chk("CPU core-days M=1999, N=1.2e10", 1.2e10*1999/4.8e7/86400, 5.8, 0.05)
for b, e, c, cT in ((0, .0500, .9707, .9741), (.25, .0547, .9515, .9557), (.5, .0693, .9237, .9279), (1.0, .1327, .8309, .8360), (2.0, .4154, .4836, .5133)):
    chk("e_I(beta=%g)" % b, float(type1_of_beta(b)), e, 6e-5)
    chk("c_min(beta=%g)" % b, coverage_min_over_nu(b, xs, ucb)[0], c, 6e-4)
    chk("c(nu=T,beta=%g)" % b, coverage_worst_phase(T, b, xs, ucb), cT, 6e-4)
chk("gamma for beta=0.25 (omega=1, N=1.2e10) x1e6", q_rms_for_beta(0.25, 1.2e10, 10**4)*1e6, 1.6, 0.03)
# --- 12.5 / Lemma 5 (V2.0.10): seed channel under the registered requirement H2(R|G) >= lambda + 32, lambda = 128, 4 devices
D1 = hash_D_bound(128 + 32); D4 = devices_D(D1, 4)
chk("D per device x2^32", D1 * 2**32, 1.0, 1e-12)
chk("D (4 devices) x1e10", D4 * 1e10, 9.3132, 5e-4)
chk("beta <= sqrt(D), 1 device x1e5", beta_ub_from_D(D1) * 1e5, 1.5259, 5e-5)
chk("beta <= sqrt(D), 4 devices x1e5", beta_ub_from_D(D4) * 1e5, 3.0518, 5e-5)
chk("size <= 0.05 + sqrt(D)/2, 4 devices: excess x1e5", (size_ub_from_D(D4) - 0.05) * 1e5, 1.5259, 5e-5)
chk("raw bits x h_min >= 160 at 512 bits, h_min=0.3125", float(seed_entropy_ok(512, 0.3125)), 1.0, 0)
chk("raw bits x h_min >= 160 at 512 bits, h_min=0.30", float(seed_entropy_ok(512, 0.30)), 0.0, 0)
for k, v in ((1, .0500), (2, .0975), (3, .1426), (20, .6415)):
    chk("best-of-%d size 1-0.95^k" % k, best_of_k_size(k), v, 5e-5)
chk("retry on theta-only failure, f=1e-3: size bound", retry_size_ub(1e-3), 0.05005, 5e-6)
sys.exit(1 if bad else 0)
