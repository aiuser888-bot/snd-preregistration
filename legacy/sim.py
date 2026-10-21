import numpy as np
from scipy import stats
rng=np.random.default_rng(7)
xs=np.load("fc_x.npy"); ucb=np.load("fc_ucb.npy")
L,K=200,200; N=L*K
j=np.arange(L); base=2*np.pi*j/L
crit=stats.chi2.ppf(.95,2)
def blocks_U(x,theta):          # x: (K,L) in {0,1}; theta (K,L) angles aligned to units
    f=1.25*x-0.25
    U=(f*np.exp(1j*theta)).sum(axis=(0,1))
    S2=((f-f.mean(axis=1,keepdims=True))**2).sum(axis=1)
    V=(L/(2*(L-1))*S2).sum()
    return U,V
# ---- (a) H0 with within-block drift + block shifts, random permutation of equispaced angles
R=4000; Lam=np.empty(R); Ua=np.empty(R); Ub=np.empty(R); Vs=np.empty(R)
for r in range(R):
    d=rng.uniform(-0.2,0.2,(K,1)); s=rng.uniform(-0.05,0.05,(K,1))
    p=np.clip(0.2+s+d*(j/L-0.5),0.01,0.99)
    x=(rng.random((K,L))<p).astype(float)
    alpha=rng.uniform(0,2*np.pi,(K,1)); th=alpha+base
    perm=np.argsort(rng.random((K,L)),axis=1); th=np.take_along_axis(th,perm,axis=1)
    U,V=blocks_U(x,th); Ua[r],Ub[r],Vs[r]=U.real,U.imag,V; Lam[r]=abs(U)**2/V
print("(a) H0 drift: mean Lam=%.3f (2), P(Lam>5.99)=%.4f  Var(Ua)/E[V]=%.3f Var(Ub)/E[V]=%.3f corr=%.3f"%(Lam.mean(),(Lam>crit).mean(),Ua.var()/Vs.mean(),Ub.var()/Vs.mean(),np.corrcoef(Ua,Ub)[0,1]))
print("    naive iid-model variance ratio (N*.25/2)/E[V]=%.3f"%(N*0.25/2/Vs.mean()))
KS=stats.kstest(Lam,"chi2",args=(2,)); print("    KS vs chi2_2 p=%.3f"%KS.pvalue)
# ---- (b) H1 iid x, coverage of FC-Rice UCB and exclusion prob
def one(rtrue_sigma,Rn,sig):
    out=[];exc=[];rej=[]
    for _ in range(Rn):
        alpha=rng.uniform(0,2*np.pi,(K,1)); th=alpha+base
        phi=rng.uniform(0,2*np.pi); c=(rtrue_sigma*sig)*np.cos(th-phi)
        p=0.2*(1+c); x=(rng.random((K,L))<p).astype(float)
        U,V=blocks_U(x,th); mu=U/V; sh=1/np.sqrt(V); xh=abs(mu)/sh
        u=np.interp(xh,xs,ucb)*sh
        out.append(u>=rtrue_sigma*sig); exc.append(u<RMIN); rej.append(abs(U)**2/V>crit)
    return np.mean(out),np.mean(exc),np.mean(rej)
EV=Vs.mean(); sig=1/np.sqrt(EV); RMIN=3.791*sig
print("(b) sigma=%.4e  RMIN=%.4e"%(sig,RMIN))
for nu in (0,0.5,1,2,3,3.791,5):
    cv,ex,rj=one(nu,3000,sig); print("  r=%.3f sigma: coverage=%.4f (se %.4f) P(UCB<rmin)=%.3f P(reject)=%.3f"%(nu,cv,np.sqrt(.95*.05/3000),ex,rj))
