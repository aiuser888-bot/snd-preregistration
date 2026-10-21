import numpy as np
from scipy.special import ive
from scipy import stats
dx=0.004; xs=np.arange(dx,22,dx); dn=0.005; nus=np.arange(0,20+dn/2,dn)
LP=np.empty((len(nus),len(xs)),dtype=np.float64)
for i,nu in enumerate(nus): LP[i]=np.log(xs)-(xs-nu)**2/2+np.log(ive(0,xs*nu))
LPm=LP.max(axis=0); acc=np.zeros(LP.shape,bool); cov=np.empty(len(nus))
for i in range(len(nus)):
    R=LP[i]-LPm; p=np.exp(LP[i])*dx; o=np.argsort(-R); cm=np.cumsum(p[o]); k=np.searchsorted(cm,0.95)
    acc[i]=R>=R[o[k]]; cov[i]=p[acc[i]].sum()
ucb=np.array([nus[acc[:,j]].max() for j in range(len(xs))])
print("cov range",cov.min(),cov.max())
for v in (0.1,1,1.794,2.448,3,5,10,15,19,21): print(v,round(float(np.interp(v,xs,ucb)),3),round(float(np.interp(v,xs,ucb))-v,3))
np.save("fc_x2.npy",xs);np.save("fc_ucb2.npy",ucb)
# --- robustness of centering: position-dependent bias delta_j correlated with drift
rng=np.random.default_rng(5); L,K=100,400; j=np.arange(L); base=2*np.pi*j/L
drift=(j/L-0.5)*0.4          # common within-block position drift in P(x)
dj=0.03*np.sign(drift)       # position-dependent angle leakage, correlated with drift
def run(mode,R=3000):
    rej=0;crit=stats.chi2.ppf(.95,2)
    for _ in range(R):
        s=rng.uniform(-0.03,0.03,(K,1)); p=np.clip(0.2+s+drift,.01,.99); f=1.25*(rng.random((K,L))<p)-.25
        a=rng.uniform(0,2*np.pi,(K,1)); th=a+base; th=np.take_along_axis(th,np.argsort(rng.random((K,L)),axis=1),axis=1)
        th=np.where(rng.random((K,L))<np.abs(dj)[None,:]*np.ones((K,1)),0.0,th) if False else np.where(rng.random((K,L))<np.abs(np.take(dj,np.argsort(np.argsort(np.zeros(L))))) ,0.0,th)
        if mode=="raw": ft=f
        elif mode=="row": ft=f-f.mean(axis=1,keepdims=True)
        else:
            ft=f-f.mean(axis=1,keepdims=True); ft=ft-ft.mean(axis=0,keepdims=True)
        ft=ft-ft.mean(axis=1,keepdims=True)
        U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum(); rej+=abs(U)**2/V>crit
    return rej/R
print("pos-dep bias x drift type I (chi2 ref): raw=%.3f rowcentred=%.3f two-way=%.3f"%(run("raw"),run("row"),run("two")))
