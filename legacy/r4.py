# Valid test: position-dependent angle leakage q_j correlated with a common within-block drift d_j.
# seed fixed per call; L,K reduced (scaling note in report).
import numpy as np
from scipy import stats
L,K=100,400; j=np.arange(L); base=2*np.pi*j/L
d=(j/L-0.5)*0.4                      # drift in P(x), zero-mean over positions
def sim(mode,q,R,seed):
    rng=np.random.default_rng(seed); crit=stats.chi2.ppf(.95,2)
    qj=np.where(d>0,q,0.0); Lam=np.empty(R)
    for r in range(R):
        s=rng.uniform(-0.03,0.03,(K,1))
        f=1.25*(rng.random((K,L))<np.clip(0.2+s+d,.01,.99))-.25
        a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
        th=np.where(rng.random((K,L))<qj[None,:],0.0,th)    # leakage: E exp(i theta_j)=q_j
        if mode=="raw": ft=f
        elif mode=="row": ft=f-f.mean(axis=1,keepdims=True)
        else:
            ft=f-f.mean(axis=1,keepdims=True); ft=ft-ft.mean(axis=0,keepdims=True)
        ft=ft-ft.mean(axis=1,keepdims=True)
        U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum(); Lam[r]=abs(U)**2/V
    return (Lam>crit).mean(), Lam.mean()
R=2000
for q in (0.0,0.01,0.03):
    out={m:sim(m,q,R,seed=100+i) for i,m in enumerate(("raw","row","two"))}
    print("q=%.2f  type-I raw=%.3f row=%.3f two-way=%.3f   (SE %.3f) mean Lam raw=%.2f row=%.2f two=%.2f"%(q,out["raw"][0],out["row"][0],out["two"][0],np.sqrt(.05*.95/R),out["raw"][1],out["row"][1],out["two"][1]))
