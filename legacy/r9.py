# Follow-up on cell B two-way type-I excess (0.058 at R=4000): larger R, plus variance diagnostics.
import numpy as np
from scipy import stats
L,K=100,400; j=np.arange(L); base=2*np.pi*j/L; crit=stats.chi2.ppf(.95,2); g=(j/L-0.5)*0.4
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def two(f):
    r=f-f.mean(axis=1,keepdims=True); t=r-r.mean(axis=0,keepdims=True); return t-t.mean(axis=1,keepdims=True)
def run(leak,R,seed):
    rng=np.random.default_rng([20261017,seed]); rej=0; re=[];im=[];Lam=[];Vs=[]
    q=np.where(g>0,0.03,0.0) if leak=="static_half" else (np.full(L,0.03) if leak=="uniform" else np.zeros(L))
    for _ in range(R):
        f=1.25*(rng.random((K,L))<np.clip(0.2+g,.01,.99)[None,:])-.25
        a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
        th=np.where(rng.random((K,L))<q[None,:],0.0,th)
        ft=two(f); U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum()
        Lam.append(abs(U)**2/V); re.append(U.real); im.append(U.imag); Vs.append(V)
    Lam=np.array(Lam); V=np.mean(Vs)
    print("leak=%-11s R=%d: type-I %.4f %s  meanLam %.3f  Var(ReU)/V=%.3f Var(ImU)/V=%.3f  mean ReU/sqrtV=%+.3f"%(leak,R,(Lam>crit).mean(),np.round(wilson((Lam>crit).sum(),R),4),Lam.mean(),np.var(re)/V,np.var(im)/V,np.mean(re)/np.sqrt(V)))
for leak,seed in (("none",1),("uniform",2),("static_half",3)): run(leak,12000,seed)
