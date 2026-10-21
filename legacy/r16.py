# (a) quantify the old AR(1) start-up bug; (b) EXACT kappa from copula orthant probabilities vs simulated E[V] (true sigma2=0.25), rho_latent in -0.9..0.9
import numpy as np
from scipy import stats
from scipy.signal import lfilter
from scipy.stats import multivariate_normal as mvn
z=stats.norm.ppf(0.8); p=0.2
def rho_x(r):                      # exact within-block correlation of Bernoulli(.2) copula at latent corr r
    if abs(r)<1e-12: return 0.0
    c=np.array([[1,r],[r,1]]); both=1-2*stats.norm.cdf(z)+mvn(mean=[0,0],cov=c).cdf([z,z])   # P(Z1>z,Z2>z)
    return float(np.squeeze(both)-p*p)/(p*(1-p))
def kappa_exact(L,rho):
    ks=np.arange(1,L); rk=np.array([rho_x(rho**k) if rho>=0 else rho_x((-abs(rho))**k) for k in ks])
    S=np.sum((1-ks/L)*rk); return (L-1-2*S)/(L-1),S
def gen(rng,K,L,rho,buggy=False):
    e=rng.normal(size=(K,L))
    if not buggy: e[:,0]/=np.sqrt(1-rho**2)
    Z=lfilter([np.sqrt(1-rho**2)],[1,-rho],e,axis=1)
    if buggy: Z[:,0]=e[:,0]                         # old code: overwrite Z_0 with e_0 (unit variance, inconsistent with Z_1)
    return Z
def f_of(Z): return 1.25*(stats.norm.cdf(Z)>0.8)-.25
def twoway(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
# (a) bug effect on lag-1 latent correlation at positions (0,1) and on E[V]
rng=np.random.default_rng([20261025,1]); L,K=100,400
for rho in (0.9,-0.9):
    for buggy in (True,False):
        Z=gen(rng,200000,2,rho,buggy)
        print("AR(1) rho=%+.1f %s: corr(Z0,Z1)=%+.4f var(Z0)=%.3f var(Z1)=%.3f (target corr %+.1f)"%(rho,"OLD-BUGGY" if buggy else "fixed    ",np.corrcoef(Z[:,0],Z[:,1])[0,1],Z[:,0].var(),Z[:,1].var(),rho))
# (b) exact kappa vs simulated, true sigma2=0.25 (p=.2 exactly), K=400
print()
for L,K in ((100,400),(400,100)):
    for rho in (-0.9,-0.6,0.0,0.6,0.9):
        kt,S=kappa_exact(L,rho); rng=np.random.default_rng([20261025,L+int(100*rho)+1000]); Vs=[]
        for _ in range(1500): ft=twoway(f_of(gen(rng,K,L,rho))); Vs.append((L/(2*(L-1))*(ft**2).sum(axis=1)).sum())
        emp=np.mean(Vs)/(K*L*0.25/2*(1-1/K)); se=np.std(Vs)/np.sqrt(len(Vs))/(K*L*0.25/2*(1-1/K))
        print("L=%d K=%d rho_latent=%+.1f: S_L(exact)=%+.3f kappa_exact=%.4f | kappa_sim=%.4f (MC SE %.4f) diff=%+.4f (%+.1f SE) | bound L/(L-1)=%.4f"%(L,K,rho,S,kt,emp,se,emp-kt,(emp-kt)/se,L/(L-1)))
