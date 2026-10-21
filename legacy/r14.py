# Closed-form check: kappa_theory=(L-1-2*S_L)/(L-1), S_L=sum_{k=1}^{L-1}(1-k/L) rho_k  (rho_k = within-block autocorrelation of f, H0)
# vs kappa_emp = E[V]/E[V_iid]. f from Gaussian-copula Bernoulli(0.2), AR(1) latent. rho_k estimated from independent large sample (different seed).
import numpy as np
from scipy import stats
from scipy.signal import lfilter
def gen(rng,K,L,rho):
    e=rng.normal(size=(K,L)); e[:,0]/=np.sqrt(1-rho**2) if rho>0 else 1
    Z=lfilter([np.sqrt(1-rho**2)],[1,-rho],e,axis=1) if rho>0 else e
    return 1.25*(stats.norm.cdf(Z)>0.8)-.25
def twoway(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
for L,K in ((100,400),(400,100)):
    for rho in (0.0,0.6,0.9):
        rng=np.random.default_rng([20261023,L+int(10*rho)])
        # rho_k from a big independent sample (long single sequences, mean 0 since E f=0): use many rows
        big=gen(rng,4000,L,rho); big=big-0.0
        sig2=(big**2).mean(); rk=np.array([ (big[:,:L-k]*big[:,k:]).mean()/sig2 for k in range(1,L)])
        S=np.sum((1-np.arange(1,L)/L)*rk); kap_th=(L-1-2*S)/(L-1)
        Vs=[];
        for _ in range(1500):
            ft=twoway(gen(rng,K,L,rho)); Vs.append((L/(2*(L-1))*(ft**2).sum(axis=1)).sum())
        EV_iid=K*L*sig2/2*(1-1/K)
        print("L=%d K=%d rho_latent=%.1f: S_L=%.3f kappa_theory=%.4f | kappa_emp=E[V]/EV_iid=%.4f (SE %.4f)  [sigma2=%.4f]"%(L,K,rho,S,kap_th,np.mean(Vs)/EV_iid,np.std(Vs)/np.sqrt(len(Vs))/EV_iid,sig2))
