import numpy as np
from scipy import stats
xs=np.load("fc_x.npy"); ucb=np.load("fc_ucb.npy")
def pex(T): return 1-np.exp(-np.interp(T,ucb,xs)**2/2)
crit=stats.chi2.ppf(.95,2)
for k in (1.25e9,1.375e9,1.5e9):
    T=1e-4*np.sqrt(k); N=2*k/0.25; mo=N*0.01/4/(30.44*86400)
    print("k=%.3e T=%.3f Pexcl=%.3f power=%.3f  N(Ef2=.25)=%.3e months(4 dev)=%.1f  N(Ef2=.15)=%.2e months=%.1f"%(k,T,pex(T),1-stats.ncx2.cdf(crit,2,T*T),N,mo,2*k/0.15,2*k/0.15*0.01/4/(30.44*86400)))
print("ops perm: N=1.1e10 * M=1e5 =",1.1e10*1e5)
print("binomial B=2e4 SE=",np.sqrt(.05*.95/2e4)," 3.5SE interval:",.05-3.5*np.sqrt(.05*.95/2e4),.05+3.5*np.sqrt(.05*.95/2e4))
print("Xhat range check: p<=.05 <-> xhat>=",np.sqrt(crit))
# exact randomization test small-scale type I under drift+shift
rng=np.random.default_rng(3)
K,L,M,R=12,60,399,1500; j=np.arange(L); base=2*np.pi*j/L
def lam(f,th):
    U=(f*np.exp(1j*th)).sum(); S2=((f-f.mean(axis=1,keepdims=True))**2).sum(); V=L/(2*(L-1))*S2
    return abs(U)**2/V
rej=0
for r in range(R):
    d=rng.uniform(-0.3,0.3,(K,1)); s=rng.uniform(-0.08,0.08,(K,1))
    p=np.clip(0.2+s+d*(j/L-.5),.01,.99); x=(rng.random((K,L))<p); f=1.25*x-.25
    def draw():
        a=rng.uniform(0,2*np.pi,(K,1)); th=a+base; pm=np.argsort(rng.random((K,L)),axis=1); return np.take_along_axis(th,pm,axis=1)
    lo=lam(f,draw()); cnt=sum(lam(f,draw())>=lo for _ in range(M)); rej+=((1+cnt)/(M+1)<=0.05)
print("randomization test type I (drift+shift, M=%d): %.4f  (SE %.4f)"%(M,rej/R,np.sqrt(.05*.95/R)))
# modulo bias illustration: Fisher-Yates with 32-bit modulo for L=10^4
n=10**4; print("modulo bias 32-bit, L=1e4: relative bias of smallest residues ~",n/2**32)
