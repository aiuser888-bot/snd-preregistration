import numpy as np
from scipy import stats
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy"); crit=stats.chi2.ppf(.95,2)
L,K=100,400; j=np.arange(L); base=2*np.pi*j/L; nu=3.0
for phi in (0.0,np.pi/4,np.pi/2,2.0):
    rng=np.random.default_rng(400); EV=K*L*0.25/2; sig=1/np.sqrt(EV); rej=cov=0;R=2500
    for _ in range(R):
        th=np.take_along_axis(np.broadcast_to(base,(K,L)),np.argsort(rng.random((K,L)),axis=1),axis=1)   # alpha_b = 0 (fixed)
        p=0.2*(1+nu*sig*np.cos(th-phi)); f=1.25*(rng.random((K,L))<p)-.25; ft=f-f.mean(axis=1,keepdims=True)
        U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum(); xh=abs(U)/np.sqrt(V)
        rej+=xh**2>crit; cov+=np.interp(xh,xs,ucb)>=nu
    print("fixed alpha, phi=%.3f: power=%.3f coverage=%.3f (SE %.3f)"%(phi,rej/R,cov/R,0.008))
