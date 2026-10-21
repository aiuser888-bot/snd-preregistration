# Reproduction criterion: success probability under different "truths" for the replication, thresholds chi2_{2,.95} vs .99.
# Gaussian large-K model, seed [20261021,1], M=4e6. Units: sigma of ORIGINAL experiment (V_target). Orig true |zeta|=nu0*s, real axis.
import numpy as np
from scipy import stats
rng=np.random.default_rng([20261021,1]); M=4_000_000; xc=np.sqrt(stats.chi2.ppf(.95,2))
th={"0.95":stats.chi2.ppf(.95,2),"0.99":stats.chi2.ppf(.99,2)}
def run(nu0,mult,rep_truth):
    srep=1/np.sqrt(mult); zo=rng.normal(size=(M,2))+[nu0,0]
    zr=rng.normal(size=(M,2))*srep+np.array(rep_truth)*1.0
    so=np.hypot(*zo.T)>=xc; sr=np.hypot(*zr.T)/srep>=xc; d2=((zr-zo)**2).sum(axis=1)/(1+srep**2)
    return {k:((d2<=v)&sr&so).sum()/so.sum() for k,v in th.items()}
nu0=3.708
cases=[("replication truth = same effect",[nu0,0]),("rep amplitude 0.75x, same phase",[0.75*nu0,0]),("rep amplitude 0.5x, same phase",[0.5*nu0,0]),("rep same amplitude, phase +60deg",[nu0*np.cos(np.pi/3),nu0*np.sin(np.pi/3)]),("rep same amplitude, phase +90deg",[0,nu0])]
for mult in (1.0,2.0):
    print("\nV_rep = %dx V_orig, original true effect = r_min (nu=%.3f), conditioned on original significant"%(mult,nu0))
    for name,t in cases:
        # rep truth in the SAME units as orig sigma; rep estimator noise scaled by srep
        r=run(nu0,mult,t); print("  %-36s success@0.95=%.3f  @0.99=%.3f"%(name,r["0.95"],r["0.99"]))
print("\nchi2 thresholds: 0.95 -> %.3f (|diff| <= %.2f sigma_diff), 0.99 -> %.3f (%.2f sigma_diff)"%(th["0.95"],np.sqrt(th["0.95"]),th["0.99"],np.sqrt(th["0.99"])))
