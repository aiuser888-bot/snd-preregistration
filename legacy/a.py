import numpy as np
from scipy import stats, optimize
rng=np.random.default_rng(1)
# A. C_smooth sensitivity: f in {1 w.p. .2, -.25 w.p. .8}, L=1e4
L=10**4
f=np.where(rng.random(L)<0.2,1.0,-0.25)
Z=f-f.mean()
Y=0.99*Z/np.sqrt(0.99**2+np.sum(Z**2))
print("A: sum Z^2 =",np.sum(Z**2)," sum Y^2 (smooth)=",np.sum(Y**2)," E[Y^2]=",np.mean(Y**2))
K=int(1.15e6)
k_smooth=K*(L/(2*(L-1)))*np.sum(Y**2); k_plain=K*(L/(2*(L-1)))*np.sum(Z**2)
print("   k_smooth=%.3e sigma=%.2e | k_plain=%.3e sigma=%.2e"%(k_smooth,1/np.sqrt(k_smooth),k_plain,1/np.sqrt(k_plain)))
# C. Rice pivot upper limit
def nup_table(xs):
    out=np.zeros_like(xs)
    for i,x in enumerate(xs):
        if stats.ncx2.cdf(x*x,2,0.0)<=0.05: out[i]=0;continue
        g=lambda lam: stats.ncx2.cdf(x*x,2,lam)-0.05
        hi=1.0
        while g(hi)>0: hi*=2
        out[i]=np.sqrt(optimize.brentq(g,0,hi))
    return out
xs=np.linspace(0,12,2401); nu=nup_table(xs)
np.save("xs.npy",xs);np.save("nu.npy",nu)
print("C: nu_up at xhat=0.32,1,1.794,2,2.5:",[float(np.interp(v,xs,nu)) for v in (0.32,1,1.794,2,2.5)])
def excl_prob(T,method):
    if method=="rice":
        xstar=np.interp(T,nu,xs)  # nu monotone increasing in x beyond floor
    else:
        xstar=T-np.sqrt(stats.chi2.ppf(.95,2))
    return 1-np.exp(-max(xstar,0)**2/2)
for m in("rice","disk"):
    T=optimize.brentq(lambda T:excl_prob(T,m)-0.8,3,6)
    sig=1e-4/T;k=1/sig**2
    print("   %s: T=%.3f sigma=%.3e k=%.3e N(E[Y2]=.25)=%.3e"%(m,T,sig,k,8*k))
# at N=1.15e10
k=1.15e10*0.125; sig=1/np.sqrt(k); T=1e-4/sig
print("   N=1.15e10 sigma=%.3e T=%.3f P_excl rice=%.3f disk=%.3f"%(sig,T,excl_prob(T,"rice"),excl_prob(T,"disk")))
# power at eps=1e-4: noncentral chi2 df2 lam=r^2 k
crit=stats.chi2.ppf(.95,2)
for N in (1e10,1.15e10):
    lam=1e-8*N*0.125; print("   N=%g lambda=%.2f power=%.3f"%(N,lam,1-stats.ncx2.cdf(crit,2,lam)))
# D. Holm resolution / MCSE
print("D: Holm smallest threshold",1e-3/66," min bootstrap p with 1e4 reps",1/(1e4+1))
print("   LB induced by block de-meaning rho=-1/(L-1)=",-1/(L-1)," SE at N=1.15e10:",1/np.sqrt(1.15e10)," Q(20)~",1.15e10*20*(1/(L-1))**2)
print("   MCSE of I_ij with 1e6 sims: sd~sqrt(2)*L/2*... =",np.sqrt(2)*2500*1/np.sqrt(1e6))
