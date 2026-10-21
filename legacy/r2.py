import numpy as np
from scipy import stats
rng=np.random.default_rng(11)
xs=np.load("fc_x.npy"); ucb=np.load("fc_ucb.npy")
print("large-x extrapolation: u(x)-x at x=10,14,14.9:",[round(float(np.interp(v,xs,ucb)-v),3) for v in (10,14,14.9)]," (spec uses +1.96)")
L,K=100,400; j=np.arange(L); base=2*np.pi*j/L
def angles(delta=0.0,fixed_alpha=False):
    a=np.zeros((K,1)) if fixed_alpha else rng.uniform(0,2*np.pi,(K,1))
    th=a+base; pm=np.argsort(rng.random((K,L)),axis=1); th=np.take_along_axis(th,pm,axis=1)
    if delta>0:   # bias: position-dependent first-harmonic leakage (positions j<L/2 favour angle 0)
        mask=rng.random((K,L))<delta; th=np.where(mask,0.0,th)
    return th
def Tpos(th):
    z=np.exp(1j*th).mean(axis=0); return K*np.sum(abs(z)**2)
# T_pos null: Gamma(L,1)?
T0=np.array([Tpos(angles()) for _ in range(3000)])
print("Tpos null: mean %.2f var %.2f (Gamma(L-1?) L=%d)"%(T0.mean(),T0.var(),L))
print("   exact mean: expected K*L*(1/K)*(1-1/(L)) ..."," ; KS vs Gamma(L-1): p=%.3f"%stats.kstest(T0,"gamma",args=(L-1,)).pvalue)
# detection power of Tpos vs bias delta
thr=np.quantile(T0,1-1e-3)
for d in (0.0,0.001,0.003,0.01):
    T=np.array([Tpos(angles(d)) for _ in range(300)]); print("   delta=%.3f: P(Tpos>thr(1e-3))=%.2f"%(d,(T>thr).mean()))
# type-I under delta x drift (position-dependent drift m_j)
def typeI(delta,m,R=3000):
    rej=0;crit=stats.chi2.ppf(.95,2)
    for _ in range(R):
        p=np.clip(0.2+m*(j/L-0.5)*np.ones((K,1)),.01,.99); x=(rng.random((K,L))<p); f=1.25*x-.25
        th=angles(delta); U=(f*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*((f-f.mean(axis=1,keepdims=True))**2).sum(axis=1)).sum()
        rej+= abs(U)**2/V>crit
    return rej/R
for d,m in ((0,0.3),(0.003,0.3),(0.01,0.3),(0.03,0.3)):
    print("delta=%.3f m=%.1f  type-I (chi2 ref)=%.3f"%(d,m,typeI(d,m)))
# non-uniform alpha (fixed 0): does Lambda still ~chi2_2?
Lam=[]
for _ in range(4000):
    x=(rng.random((K,L))<0.2); f=1.25*x-.25; th=angles(0,True); U=(f*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*((f-f.mean(axis=1,keepdims=True))**2).sum(axis=1)).sum(); Lam.append(abs(U)**2/V)
print("fixed alpha: P(Lam>5.99)=%.4f KS p=%.3f"%((np.array(Lam)>5.991).mean(),stats.kstest(Lam,'chi2',args=(2,)).pvalue))
