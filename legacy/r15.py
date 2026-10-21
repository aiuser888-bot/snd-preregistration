# Direct check of Lemma 2(ii): H0, serial correlation, one block at a time (row-centred). Var(Re U_b), Var(Im U_b), Cov vs E[V_b]; also cross-block cov with correlated x.
import numpy as np
from scipy import stats
from scipy.signal import lfilter
L=100; M=200000; base=2*np.pi*np.arange(L)/L
def gen(rng,n,rho,cross=0.0):
    e=rng.normal(size=(n,L)); e[:,0]/=np.sqrt(1-rho**2) if rho>0 else 1
    Z=lfilter([np.sqrt(1-rho**2)],[1,-rho],e,axis=1) if rho>0 else e
    return 1.25*(stats.norm.cdf(Z)>0.8)-.25
for rho in (0.0,0.6,0.9):
    rng=np.random.default_rng([20261024,int(10*rho)]); ur=[];ui=[];Vb=[]
    for _ in range(M//20000):
        f=gen(rng,20000,rho); f=f-f.mean(axis=1,keepdims=True)
        a=rng.uniform(0,2*np.pi,(20000,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((20000,L)),axis=1),axis=1)
        U=(f*np.exp(1j*th)).sum(axis=1); ur.append(U.real); ui.append(U.imag); Vb.append(L/(2*(L-1))*(f**2).sum(axis=1))
    ur=np.concatenate(ur);ui=np.concatenate(ui);Vb=np.concatenate(Vb); EV=Vb.mean(); se=np.sqrt(2/len(ur))*1.0
    print("rho=%.1f: E[V_b]=%.3f  Var(Re)/E[V]=%.4f Var(Im)/E[V]=%.4f  (SE~%.4f) corr(Re,Im)=%+.4f"%(rho,EV,ur.var()/EV,ui.var()/EV,se*1.2,np.corrcoef(ur,ui)[0,1]))
# cross-block covariance with cross-block-correlated x: two blocks, x of block 2 strongly correlated with block 1 (same Z), independent alpha_b
rng=np.random.default_rng([20261024,99]); n=200000
e=rng.normal(size=(n,L)); Z=lfilter([np.sqrt(1-0.9**2)],[1,-0.9],e,axis=1)
f1=1.25*(stats.norm.cdf(Z)>0.8)-.25; f2=1.25*(stats.norm.cdf(0.95*Z+np.sqrt(1-0.95**2)*rng.normal(size=(n,L)))>0.8)-.25
def U(f):
    a=rng.uniform(0,2*np.pi,(n,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((n,L)),axis=1),axis=1); return (f*np.exp(1j*th)).sum(axis=1)
u1,u2=U(f1),U(f2)
print("cross-block (x correlated 0.95 between blocks): corr(Re U1,Re U2)=%+.4f corr(Re U1,Im U2)=%+.4f (SE %.4f); corr of f's block means=%+.3f"%(np.corrcoef(u1.real,u2.real)[0,1],np.corrcoef(u1.real,u2.imag)[0,1],1/np.sqrt(n),np.corrcoef(f1.mean(axis=1),f2.mean(axis=1))[0,1]))
