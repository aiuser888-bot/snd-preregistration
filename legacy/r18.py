# Coverage grid incl. NEGATIVE serial correlation. seeds [20261027,L,cell]. L=100,K=400,R=2000. pred uses exact kappa (r16).
import numpy as np, time, sys
exec(open('r16.py').read().split('# (a) bug effect')[0])
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy")
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def xmin(nu): return xs[np.argmax(ucb>=nu)]
import sys
L,K,R=int(sys.argv[1]),int(sys.argv[2]),1500; base=2*np.pi*np.arange(L)/L; EV0=K*L*0.25/2*(1-1/K); sig0=1/np.sqrt(EV0)
t0=time.time(); worst=(1,None); cell=0
for rho in (-0.9,0.0,0.9):
    kap=kappa_exact(L,rho)[0]
    for nu in (3.0,):
        cell+=1; rng=np.random.default_rng([20261027,L,cell]); cov=0
        for _ in range(R):
            Z=gen(rng,K,L,rho); a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
            phi=rng.uniform(0,2*np.pi); pp=0.2*(1+nu*sig0*np.cos(th-phi)); f=1.25*(stats.norm.cdf(Z)>1-pp)-.25
            ft=twoway(f); U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum()
            xh=abs(U)/np.sqrt(V); cov+=np.interp(xh,xs,ucb)/np.sqrt(V)>=nu*sig0
        lo,hi=wilson(cov,R); pred=stats.ncx2.sf(xmin(nu)**2,2,nu**2/kap); zz=(cov/R-pred)/np.sqrt(pred*(1-pred)/R)
        print("rho=%+.1f kappa=%.4f nu=%.0f L=%d K=%d R=%d nominal .95 | cov=%.4f [%.4f,%.4f] pred=%.4f z=%+.2f"%(rho,kap,nu,L,K,R,cov/R,lo,hi,pred,zz),flush=True)
        if cov/R<worst[0]: worst=(cov/R,(rho,nu,lo))
print("worst",worst,"t=%.0fs"%(time.time()-t0))
