# (1) CLT at real scale: excess kurtosis of Re U_b/sqrt(V_b) for one block L=1e4 under within-block permutation (H0, conditional on x). seeds [20261029,1]
import numpy as np
from scipy import stats
rng=np.random.default_rng([20261029,1]); L=10000
for name,f in (("Bernoulli.2",1.25*(rng.random(L)<.2)-.25),("heavy(t3)",np.clip(rng.standard_t(3,L)/6,-1,1))):
    f=f-f.mean(); V=L/(2*(L-1))*(f**2).sum(); base=2*np.pi*np.arange(L)/L; a=rng.uniform(); m=20000; out=np.empty(m)
    for i in range(m): out[i]=(f*np.cos(a+base[rng.permutation(L)])).sum()
    k=stats.kurtosis(out); se=np.sqrt(24/m)
    print("%s: Var(ReU)/V=%.4f  excess kurt of block score=%+.4f (SE %.3f) -> total over K=1e6 blocks: %+.2e"%(name,out.var()/V,k,se,k/1e6))
# (2) cross-block dependence: shared latent factor W_j (same for all blocks at slot j) and shared per-block-group factor. seeds [20261029,2..]
from scipy.signal import lfilter
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy")
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def tw(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
def xmin(nu): return xs[np.argmax(ucb>=nu)]
L,K,R=100,400,1500; base=2*np.pi*np.arange(L)/L; EV0=K*L*.25/2*(1-1/K); sig0=1/np.sqrt(EV0); rho=0.6; nu=3.0
for cell,(cs,cb) in enumerate(((0.0,0.0),(0.5,0.0),(0.0,0.5),(0.5,0.5)),1):   # cs: shared-by-slot weight, cb: shared-across-all-slots-per-time-cluster (K rows share a drifting W over blocks)
    rng=np.random.default_rng([20261029,10+cell]); cov=0; rej=0; Vs=[]
    for _ in range(R):
        e=rng.normal(size=(K,L)); e[:,0]/=np.sqrt(1-rho**2); Zi=lfilter([np.sqrt(1-rho**2)],[1,-rho],e,axis=1)
        W=lfilter([np.sqrt(1-rho**2)],[1,-rho],rng.normal(size=L)*np.r_[1/np.sqrt(1-rho**2),np.ones(L-1)])   # slot-shared AR(1)
        B=lfilter([np.sqrt(1-rho**2)],[1,-rho],rng.normal(size=K)*np.r_[1/np.sqrt(1-rho**2),np.ones(K-1)])    # block-level AR(1) across blocks (cross-device/time correlation)
        Z=np.sqrt(1-cs-cb)*Zi+np.sqrt(cs)*W[None,:]+np.sqrt(cb)*B[:,None]
        a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1); phi=rng.uniform(0,2*np.pi)
        pp=0.2*(1+nu*sig0*np.cos(th-phi)); f=1.25*(stats.norm.cdf(Z)>1-pp)-.25
        ft=tw(f); U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum(); xh=abs(U)/np.sqrt(V)
        cov+=np.interp(xh,xs,ucb)/np.sqrt(V)>=nu*sig0; Vs.append(V)
    lo,hi=wilson(cov,R); kap=np.mean(Vs)/EV0; pred=stats.ncx2.sf(xmin(nu)**2,2,nu**2/kap)
    print("shared(slot)=%.1f shared(block-AR)=%.1f latent rho=.6 nu=3 R=%d: kappa_emp=%.3f cov=%.4f [%.4f,%.4f] pred=%.4f z=%+.2f"%(cs,cb,R,kap,cov/R,lo,hi,pred,(cov/R-pred)/np.sqrt(pred*(1-pred)/R)),flush=True)
