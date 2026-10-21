# Exact (Rice) coverage and power of the UCB table, vs the pooled simulation cells in r7.py output.
import numpy as np
from scipy import stats
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy"); ok=xs<=18
crit=stats.chi2.ppf(.95,2)
# u is nondecreasing on ok-range -> {u(x)>=nu} = {x>=xmin(nu)}
def xmin(nu): return np.interp(nu,ucb[ok],xs[ok])
def cov_exact(nu): return stats.ncx2.sf(xmin(nu)**2,2,nu**2)
def pow_exact(nu): return stats.ncx2.sf(crit,2,nu**2)
print("monotone ucb on x<=18:",bool(np.all(np.diff(ucb[ok])>=0)))
nus=np.arange(0.05,16,0.05); cv=np.array([cov_exact(n) for n in nus])
print("exact coverage curve: min %.5f at nu=%.2f ; at nu=16 -> %.5f"%(cv.min(),nus[cv.argmin()],cv[-1]))
for n in (0.5,1,2,3,4,5,8,12): print("  nu=%.1f exact coverage %.4f  exact power %.4f"%(n,cov_exact(n),pow_exact(n)))
# pooled simulation (r7.py fixed-alpha, 8 phis x 3000 = 24000 per nu): coverage / power from the printed cells
sim={2.0:dict(cov=[.9787,.9817,.9823,.9783,.9757,.9773,.9777,.9817],pw=[.408,.425,.421,.424,.423,.415,.401,.421]),
     3.0:dict(cov=[.9727,.9747,.9713,.9760,.9710,.9743,.9673,.9737],pw=[.778,.783,.765,.784,.785,.778,.760,.759]),
     4.0:dict(cov=[.9740,.9700,.9727,.9720,.9743,.9763,.9753,.9750],pw=[.956,.952,.954,.953,.957,.956,.954,.954])}
for nu,d in sim.items():
    n=24000
    for name,vals,ex in (("coverage",d["cov"],cov_exact(nu)),("power",d["pw"],pow_exact(nu))):
        p=np.mean(vals); z=(p-ex)/np.sqrt(ex*(1-ex)/n)
        print("nu=%.0f %s: pooled sim %.4f  exact %.4f  z=%+.2f  (SE %.4f)"%(nu,name,p,ex,z,np.sqrt(ex*(1-ex)/n)))
