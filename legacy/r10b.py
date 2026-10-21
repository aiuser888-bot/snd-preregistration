# Recheck: is mean(Re U)/sqrtV for the heavy-tail ftilde really 0? New seeds, M=4e5 each, 3 independent randomization seeds.
import numpy as np
L,K=100,20; base=2*np.pi*np.arange(L)/L
def twoway(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
rng=np.random.default_rng([20261018,1]); _=twoway(1.25*(rng.random((K,L))<0.2)-.25)      # same x_a draw as r10, to reproduce x_b
ft=twoway(np.clip(rng.standard_cauchy((K,L))*0.05,-1,1)); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum()
print("same dataset as r10: V=%.4f  (r10 printed 58.1588)"%V)
M=400000
for seed in (21,22,23):
    r2=np.random.default_rng([20261020,seed]); ur=np.empty(M); ui=np.empty(M)
    for s in range(0,M,2000):
        m=min(2000,M-s); a=r2.uniform(0,2*np.pi,(m,K,1)); perm=np.argsort(r2.random((m,K,L)),axis=2)
        th=np.take_along_axis(a+base[None,None,:],perm,axis=2); U=(ft[None]*np.exp(1j*th)).sum(axis=(1,2)); ur[s:s+m]=U.real; ui[s:s+m]=U.imag
    print("seed %d: mean Re/sqrtV=%+.4f  mean Im/sqrtV=%+.4f  (SE %.4f)  Var(Re)/V=%.4f Var(Im)/V=%.4f (SE %.4f)"%(seed,ur.mean()/np.sqrt(V),ui.mean()/np.sqrt(V),1/np.sqrt(M),ur.var()/V,ui.var()/V,np.sqrt(2/M)))
