# H1 under serial correlation: marginal-tilt Gaussian-copula model. x_{b,j}=1[Phi(Z_bj) > 1-p_bj], p_bj=0.2(1+nu*sig0*cos(theta_bj-phi)),
# Z AR(1) along position j within block (independent across blocks), theta random permutation of equispaced angles (alpha_b uniform).
# Statistic = v11 two-way centred U, V.  Prediction: E[V]=kappa*EV_iid, nu_eff=nu/sqrt(kappa); coverage/power(Rice at nu_eff).
import numpy as np, time
from scipy import stats
from scipy.signal import lfilter
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy"); ok=xs<=18; crit=stats.chi2.ppf(.95,2)
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def twoway(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
def cov_exact(nu): return stats.ncx2.sf(np.interp(nu,ucb[ok],xs[ok])**2,2,nu**2)
def pow_exact(nu): return stats.ncx2.sf(crit,2,nu**2)
def sim(L,K,rho,nu,R,seed):
    rng=np.random.default_rng([20261022,seed]); base=2*np.pi*np.arange(L)/L; EV0=K*L*0.25/2*(1-1/K); sig0=1/np.sqrt(EV0)
    rej=cov=0; Vs=[]
    for _ in range(R):
        e=rng.normal(size=(K,L)); e[:,0]/=np.sqrt(1-rho**2)          # so that Z_0 ~ N(0,1): stationary AR(1)
        Z=lfilter([np.sqrt(1-rho**2)],[1,-rho],e,axis=1)
        a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
        phi=rng.uniform(0,2*np.pi); p=0.2*(1+nu*sig0*np.cos(th-phi))
        f=1.25*(stats.norm.cdf(Z)>1-p)-.25
        ft=twoway(f); U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum()
        xh=abs(U)/np.sqrt(V); rej+=xh*xh>crit; cov+=np.interp(xh,xs,ucb)*(1/np.sqrt(V))>=nu*sig0; Vs.append(V)
    return rej,cov,np.mean(Vs),EV0
t0=time.time(); R=2500
for (L,K) in ((100,400),(400,100)):
    for rho in (0.0,0.6,0.9):
        out=[]
        _,_,Vh,EV0=sim(L,K,rho,0.0,300,1)          # kappa from H0 mean V
        kap=Vh/EV0
        line="L=%d K=%d rho=%.1f kappa_emp=%.3f (nu_eff/nu=%.3f)"%(L,K,rho,kap,1/np.sqrt(kap))
        for nu in (0.0,3.0):
            r,c,_,_=sim(L,K,rho,nu,R,int(100*rho)+int(nu*10)+L)
            ne=nu/np.sqrt(kap)
            line+=" | nu=%.0f: reject %.3f[%.3f,%.3f] cov %.4f[%.4f,%.4f]"%(nu,r/R,*wilson(r,R),c/R,*wilson(c,R))
            if nu>0: line+="  pred@nu_eff: power %.3f cov %.4f ; pred@nu: power %.3f cov %.4f"%(pow_exact(ne),cov_exact(ne),pow_exact(nu),cov_exact(nu))
        print(line,flush=True)
print("t=%.0fs"%(time.time()-t0))
