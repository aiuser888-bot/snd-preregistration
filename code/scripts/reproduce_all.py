"""Reproduce the headline reduced-scale numbers with the reference implementation. Writes results/reproduce.json."""
import json, os, sys, time, numpy as np
from scipy import stats
from snd.fcrice import load_table
from snd.kappa import kappa_exact, coverage_prediction
from snd.accept import simulate
from snd.design import seed_to_key
from snd.stats import center_two_way
from snd.randtest import randomization_test

root = os.path.join(os.path.dirname(__file__), "..", "..", "results")
xs, ucb = load_table(os.path.join(root, "fc_rice_table.npz"))
out = {"kappa": [], "coverage": [], "size": None}
t0 = time.time()
for L in (100, 400):
    for rho in (-0.9, -0.6, 0.0, 0.6, 0.9):
        k, S = kappa_exact(L, rho); out["kappa"].append(dict(L=L, rho=rho, kappa=k, S_L=S))
R = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
zs = []
for i, rho in enumerate((-0.9, -0.6, 0.0, 0.6, 0.9)):
    kap, _ = kappa_exact(100, rho)
    for nu0 in (2.0, 3.0, 4.0):
        r = simulate(100, 400, rho, nu0, R, 100 + 10 * i + int(nu0), xs, ucb)
        nu = nu0 * np.sqrt(kap)                         # nu = r / sigma_hat (actual V)
        pred = coverage_prediction(nu, kap, xs, ucb)
        z = (r["coverage"] - pred) / np.sqrt(pred * (1 - pred) / R); zs.append(z)
        r.update(pred=pred, z=float(z), kappa_exact=kap, wilson_cov=list(r["wilson_cov"])); out["coverage"].append(r)
        print("rho=%+.1f nu0=%g cov=%.4f [%.4f,%.4f] pred=%.4f z=%+.2f" % (rho, nu0, r["coverage"], *r["wilson_cov"], pred, z), flush=True)
out["coverage_summary"] = dict(R=R, cells=len(zs), min_wilson_lower=min(c["wilson_cov"][0] for c in out["coverage"]),
                               min_coverage=min(c["coverage"] for c in out["coverage"]), max_abs_z=float(np.max(np.abs(zs))), mean_z=float(np.mean(zs)))
# size of the full randomization test (reduced scale, M=199 -> exact level 0.05 at #<=9)
rng = np.random.default_rng(99); Rsz, K, L = 400, 20, 50; rej = 0
for r in range(Rsz):
    ft = center_two_way(1.25 * (rng.random((K, L)) < 0.2) - 0.25)
    rej += randomization_test(ft, seed_to_key(b"size-%08d-" % r + bytes(range(56))), M=199)["reject"]
out["size"] = dict(R=Rsz, K=K, L=L, M=199, reject_rate=rej / Rsz, se=float(np.sqrt(.05 * .95 / Rsz)))
out["seconds"] = time.time() - t0
json.dump(out, open(os.path.join(root, "reproduce.json"), "w"), indent=1, default=float)
print(json.dumps(out["coverage_summary"]), out["size"])
