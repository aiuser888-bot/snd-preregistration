"""Build the FC-Rice table, check the protocol's checkpoints, save npz and SHA-256."""
import hashlib, sys, numpy as np
from snd.fcrice import build_table, TABLE_CHECKPOINTS
out = sys.argv[1] if len(sys.argv) > 1 else "fc_rice_table.npz"
xs, ucb, mass = build_table()
np.savez(out, xs=xs, ucb=ucb)
print("acceptance mass range: %.5f - %.5f" % (mass.min(), mass.max()))
ok = True
for x, v in TABLE_CHECKPOINTS.items():
    got = float(np.interp(x, xs, ucb)); good = abs(got - v) <= 0.01; ok &= good
    print("u(%g) = %.3f (protocol %.3f) %s" % (x, got, v, "OK" if good else "MISMATCH"))
print("sha256", hashlib.sha256(open(out, "rb").read()).hexdigest())
sys.exit(0 if ok and mass.min() >= 0.9500 and mass.max() <= 0.9505 else 1)
