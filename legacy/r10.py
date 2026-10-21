# (i) Conditional-on-x identity: over the randomization (alpha_b uniform + uniform permutation of equispaced angles),
#     Var(Re U)=Var(Im U)=V, Cov=0, E U=0, for ANY fixed x-derived ftilde (row sums 0). Test on two very different ftilde.
# (ii) E[V] = N*sigma_f^2*(1-1/K)/2 for iid x.   Seeds: default_rng([20261018, tag]).
import numpy as np
L,K,M=100,20,200000; base=2*np.pi*np.arange(L)/L
def twoway(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
rng=np.random.default_rng([20261018,1])
x_a=twoway(1.25*(rng.random((K,L))<0.2)-.25)                       # two-point f
x_b=twoway(np.clip(rng.standard_cauchy((K,L))*0.05,-1,1))          # heavy-tailed, clipped to |f|<=1
for name,ft in (("two-point f ",x_a),("heavy-tail f",x_b)):
    V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum(); ur=np.empty(M); ui=np.empty(M); r2=np.random.default_rng([20261018,2])
    for s in range(0,M,2000):
        m=min(2000,M-s); a=r2.uniform(0,2*np.pi,(m,K,1)); perm=np.argsort(r2.random((m,K,L)),axis=2)
        th=np.take_along_axis(a+base[None,None,:],perm,axis=2); U=(ft[None]*np.exp(1j*th)).sum(axis=(1,2)); ur[s:s+m]=U.real; ui[s:s+m]=U.imag
    se=np.sqrt(2/M)
    print("%s: V=%.4f | Var(Re)/V=%.4f Var(Im)/V=%.4f (SE %.4f) corr=%+.4f mean(Re,Im)/sqrtV=(%+.4f,%+.4f) (SE %.4f)"%(name,V,ur.var()/V,ui.var()/V,se,np.corrcoef(ur,ui)[0,1],ur.mean()/np.sqrt(V),ui.mean()/np.sqrt(V),1/np.sqrt(M)))
# (ii) E[V]
for (L2,K2) in ((100,400),(50,20)):
    r3=np.random.default_rng([20261018,3]); Vs=[]
    for _ in range(3000):
        f=1.25*(r3.random((K2,L2))<0.2)-.25; ft=twoway(f); Vs.append((L2/(2*(L2-1))*(ft**2).sum(axis=1)).sum())
    N=K2*L2; pred=N*0.25/2*(1-1/K2); print("L=%d K=%d: mean V=%.2f (SE %.2f) predicted N*sigma^2/2*(1-1/K)=%.2f ; uncorrected N*sigma^2/2=%.2f"%(L2,K2,np.mean(Vs),np.std(Vs)/np.sqrt(len(Vs)),pred,N*0.25/2))
