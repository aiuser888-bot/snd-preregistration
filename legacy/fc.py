import numpy as np
from scipy.special import ive
from scipy import stats
dx=0.002; xs=np.arange(dx,15,dx); dn=0.005; nus=np.arange(0,13+dn/2,dn)
# log p(x|nu) for Rice(sigma=1)
LP=np.empty((len(nus),len(xs)))
for i,nu in enumerate(nus):
    LP[i]=np.log(xs)-(xs-nu)**2/2+np.log(ive(0,xs*nu))
LPmax=LP.max(axis=0)          # profile over nu>=0
Rc=np.empty(len(nus)); acc=np.zeros(LP.shape,bool)
for i in range(len(nus)):
    R=LP[i]-LPmax; p=np.exp(LP[i])*dx
    o=np.argsort(-R); cm=np.cumsum(p[o]); k=np.searchsorted(cm,0.95)
    Rc[i]=R[o[k]]; acc[i]=R>=Rc[i]
cov=np.array([ (np.exp(LP[i])*dx)[acc[i]].sum() for i in range(len(nus))])
print("grid coverage min/max:",cov.min(),cov.max())
ucb=np.array([nus[acc[:,j]].max() if acc[:,j].any() else np.nan for j in range(len(xs))])
lcb=np.array([nus[acc[:,j]].min() if acc[:,j].any() else np.nan for j in range(len(xs))])
np.save("fc_x.npy",xs);np.save("fc_ucb.npy",ucb)
print("UCB at x=0.1,1,1.794,2.448,3,5:",[round(float(np.interp(v,xs,ucb)),3) for v in (0.1,1,1.794,2.448,3,5)])
print("monotone UCB:",bool(np.all(np.diff(ucb)>=-0.011)))
# exclusion prob under H0
for tag,T in (("find",None),):
    pass
def pex(T):
    xstar=np.interp(T,ucb,xs)   # ucb nondecreasing
    return 1-np.exp(-xstar**2/2)
from scipy import optimize
T80=optimize.brentq(lambda T:pex(T)-0.8,2.7,4.5); print("T80 (FC)=",T80)
for E in (0.25,):
    k=(T80/1e-4)**2; print("FC: k=%.3e  N=%.3e (E f2=%.2f)"%(k,2*k/E,E))
k=1.15e10*0.125; T=1e-4*np.sqrt(k); print("N=1.15e10: T=%.3f Pexcl(FC)=%.3f"%(T,pex(T)))
crit=stats.chi2.ppf(.95,2)
for kk in (1.0e9,1.1e9,1.25e9):
    print("k=%.2e power=%.3f Pexcl=%.3f"%(kk,1-stats.ncx2.cdf(crit,2,1e-8*kk),pex(1e-4*np.sqrt(kk))))
# permutation-test vs FC consistency: p<=.05 <-> x>=2.448 ; does FC interval exclude 0 exactly there?
j=np.where(lcb>0)[0][0]; print("first x with LCB>0: %.3f (sqrt(5.99)=%.3f)"%(xs[j],np.sqrt(crit)))
