# Phase dependence of H1 power AND UCB coverage, v10 statistic (two-way centering), L=100,K=400.
# Independent seed per cell: default_rng([20261015, cell]). Wilson 95% CIs. Fixed-alpha grid + random-alpha control.
import numpy as np, sys, time
from scipy import stats
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy"); crit=stats.chi2.ppf(.95,2)
L,K=100,400; base=2*np.pi*np.arange(L)/L; EV=K*L*0.25/2; sig=1/np.sqrt(EV)
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def twoway(f):
    ft=f-f.mean(axis=1,keepdims=True); ft=ft-ft.mean(axis=0,keepdims=True); return ft-ft.mean(axis=1,keepdims=True)
def cell(nu,phi,fixed,R,cid):
    rng=np.random.default_rng([20261015,cid]); rej=cov=0
    for _ in range(R):
        a=np.zeros((K,1)) if fixed else rng.uniform(0,2*np.pi,(K,1))
        th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
        f=1.25*(rng.random((K,L))<0.2*(1+nu*sig*np.cos(th-phi)))-.25
        ft=twoway(f); U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum()
        xh=abs(U)/np.sqrt(V); rej+=xh*xh>crit; cov+=np.interp(xh,xs,ucb)>=nu
    return rej,cov
R=3000; phis=[0.0,0.3,np.pi/4,np.pi/2,2.0,np.pi,4.0,5.2]; t0=time.time(); cid=0; res={}
for fixed in (True,False):
    for nu in ((2.0,3.0,4.0) if fixed else (3.0,)):
        for phi in phis:
            res[(fixed,nu,phi)]=cell(nu,phi,fixed,R,cid); cid+=1
    print("done fixed=%s  t=%.0fs"%(fixed,time.time()-t0),flush=True)
for fixed in (True,False):
    for nu in ((2.0,3.0,4.0) if fixed else (3.0,)):
        print("\n== alpha %s, nu=%.1f, R=%d per phi =="%("FIXED(0)" if fixed else "RANDOM",nu,R))
        P=[];C=[]
        for phi in phis:
            r,c=res[(fixed,nu,phi)]; pl,ph=wilson(r,R); cl,ch=wilson(c,R); P.append(r);C.append(c)
            print("  phi=%.3f  power %.3f [%.3f,%.3f]  coverage %.4f [%.4f,%.4f]"%(phi,r/R,pl,ph,c/R,cl,ch))
        P=np.array(P);C=np.array(C)
        for name,X in (("power",P),("coverage",C)):
            tab=np.vstack([X,R-X]); chi2,pv,_,_=stats.chi2_contingency(tab,correction=False)
            print("  %s: range %.4f-%.4f  spread %.4f | homogeneity chi2(7)=%.2f p=%.3f"%(name,X.min()/R,X.max()/R,(X.max()-X.min())/R,chi2,pv))
        print("  worst-phi coverage %.4f, Wilson lower %.4f (unadjusted; 8 phis)"%(C.min()/R,wilson(C.min(),R)[0]))
# H0 sanity (nu=0) with same code
for fixed in (True,False):
    r,c=cell(0.0,0.0,fixed,R,900+int(fixed)); print("H0 alpha %s: P(reject)=%.4f [%.4f,%.4f] coverage=%.4f"%("fixed" if fixed else "random",r/R,*wilson(r,R),c/R))
print("total %.0fs"%(time.time()-t0))
