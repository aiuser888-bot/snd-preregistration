"""Section 12 (V2.0.08): analytic leakage curves.  Writes results/leak_curves.json and prints the protocol's table."""
import json, os
import numpy as np
from snd.fcrice import load_table
from snd.leak import type1_of_beta, coverage_min_over_nu, coverage_worst_phase, q_rms_for_beta

root = os.path.join(os.path.dirname(__file__), "..", "..", "results")
xs, ucb = load_table(os.path.join(root, "fc_rice_table.npz"))
T, N, L = 3.873, 1.2e10, 10 ** 4   # V_target = 1.5e9 (V2.0.09); sigma_f^2 = 0.25
rows = []
print("beta   e_I      c_min[0,15] (nu*)     c(nu=T)    gamma_max (omega=1, N=1.2e10, L=1e4)")
for b in (0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    c, n = coverage_min_over_nu(b, xs, ucb)
    row = dict(beta=b, e_I=float(type1_of_beta(b)), c_min=c, nu_at_min=n, c_at_T=coverage_worst_phase(T, b, xs, ucb),
               gamma_max_omega1=float(q_rms_for_beta(b, N, L, 1.0)))
    rows.append(row)
    print("%-5g %.4f   %.4f (%.2f)        %.4f     %.2e" % (b, row["e_I"], c, n, row["c_at_T"], row["gamma_max_omega1"]))
json.dump(dict(T=T, N=N, L=L, rows=rows), open(os.path.join(root, "leak_curves.json"), "w"), indent=1)
