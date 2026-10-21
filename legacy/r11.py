# Reproduction criterion (v10 sec.10): (i) p_rep<=.05 and (ii) |mu_rep-mu_orig|^2/(s_rep^2+s_orig^2) <= chi2_{2,.95}.
# Gaussian large-K model: mu_hat ~ N(zeta, s^2 I_2), s=1/sqrt(V). Original is only replicated if it was significant (selection!).
import numpy as np
from scipy import stats
rng=np.random.default_rng([20261019,1]); M=4_000_000
c95=stats.chi2.ppf(.95,2); c99=stats.chi2.ppf(.99,2); xcrit=np.sqrt(c95)
def run(nu_true,Vrep_mult,label):
    s=1.0; srep=s/np.sqrt(Vrep_mult)             # sigma units of the original; true |zeta| = nu_true*s
    zo=rng.normal(size=(M,2))*s+np.array([nu_true*s,0]); zr=rng.normal(size=(M,2))*srep+np.array([nu_true*s,0])
    sig_o=np.hypot(*zo.T)/s>=xcrit; sig_r=np.hypot(*zr.T)/srep>=xcrit
    d2=((zr-zo)**2).sum(axis=1)/(s*s+srep*srep)
    p_sel=sig_o.mean()
    res={"P(orig sig)":p_sel,"P(rep sig | orig sig)":(sig_r&sig_o).sum()/sig_o.sum(),
         "P(ii@.95 | orig sig)":((d2<=c95)&sig_o).sum()/sig_o.sum(),
         "P(success .95 | orig sig)":((d2<=c95)&sig_r&sig_o).sum()/sig_o.sum(),
         "P(success .99 | orig sig)":((d2<=c99)&sig_r&sig_o).sum()/sig_o.sum()}
    print("%-34s "%label+"  ".join("%s=%.3f"%(k,v) for k,v in res.items()))
print("nu = true amplitude in units of the ORIGINAL experiment's sigma (r_min => nu=3.708 at V_target); rep sigma smaller by sqrt(mult)")
for nu in (3.708,5.0,8.0):
    for mult in (1.0,2.0): run(nu,mult,"nu=%.2f V_rep=%dxV_orig"%(nu,mult))
print("(null check) criterion (ii) unconditional, no selection, same true zeta: P=%.3f (expected .95)"%(lambda zo,zr:(((zr-zo)**2).sum(axis=1)/2<=c95).mean())(rng.normal(size=(M,2)),rng.normal(size=(M,2))))
