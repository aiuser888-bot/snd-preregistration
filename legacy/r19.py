# Conditional-tilt H1: x_j ~ Bernoulli(s(u_j)*(1+delta*cos(th-phi))), u=Phi(Z) AR(1) latent, s=.8 if u>.8 else .05 (marginal mean .2). seeds [20261028,cell]
import numpy as np, time
from scipy import stats
from scipy.signal import lfilter
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy")
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def tw(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
def xmin(nu): return xs[np.argmax(ucb>=nu)]
L,K,R=100,400,1500; base=2*np.pi*np.arange(L)/L; EV0=K*L*.25/2*(1-1/K); sig0=1/np.sqrt(EV0)
def run(rho,nu,R,rng):
    cov=0;Vs=[]
    for _ in range(R):
        e=rng.normal(size=(K,L)); e[:,0]/=np.sqrt(1-rho**2); Z=lfilter([np.sqrt(1-rho**2)],[1,-rho],e,axis=1); u=stats.norm.cdf(Z)
        a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1); phi=rng.uniform(0,2*np.pi)
        q=np.where(u>.8,.8,.05)*(1+nu*sig0*np.cos(th-phi)); f=1.25*(rng.random((K,L))<q)-.25
        ft=tw(f); U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum()
        xh=abs(U)/np.sqrt(V); cov+=np.interp(xh,xs,ucb)/np.sqrt(V)>=nu*sig0; Vs.append(V)
    return cov,np.mean(Vs)/EV0
cell=0
for rho in (-0.9,0.0,0.9):
    rng=np.random.default_rng([20261028,100+int(10*rho)]); _,kap=run(rho,0.0,600,rng)
    for nu in (3.0,):
        cell+=1; rng=np.random.default_rng([20261028,cell]); c,_=run(rho,nu,R,rng); lo,hi=wilson(c,R)
        pred=stats.ncx2.sf(xmin(nu)**2,2,nu**2/kap)
        print("cond-tilt rho=%+.1f kappa_emp(H0,600)=%.4f nu=%.0f R=%d cov=%.4f [%.4f,%.4f] pred=%.4f z=%+.2f"%(rho,kap,nu,R,c/R,lo,hi,pred,(c/R-pred)/np.sqrt(pred*(1-pred)/R)),flush=True)
