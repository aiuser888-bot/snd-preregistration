import numpy as np
from scipy import stats, optimize
xs=np.load("fc_x2.npy"); ucb=np.load("fc_ucb2.npy"); ok=xs<=18
def pex(T): return 1-np.exp(-np.interp(T,ucb[ok],xs[ok])**2/2)
T80=optimize.brentq(lambda T:pex(T)-0.8,3,4); print("T80=%.4f  V80=%.4e"%(T80,(T80/1e-4)**2))
Vt=1.375e9;T=1e-4*np.sqrt(Vt);crit=stats.chi2.ppf(.95,2)
print("V_target T=%.4f Pexcl=%.4f power(r=rmin)=%.4f"%(T,pex(T),1-stats.ncx2.cdf(crit,2,T*T)))
# power / coverage of two-way-centred statistic under H1 (iid x, no drift); seed 200
L,K=100,400; j=np.arange(L); base=2*np.pi*j/L
def run(center,nu,R=3000,seed=200):
    rng=np.random.default_rng(seed); EV=K*L*0.25/2; sig=1/np.sqrt(EV); rej=cov=0
    for _ in range(R):
        a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
        phi=rng.uniform(0,2*np.pi); p=0.2*(1+nu*sig*np.cos(th-phi)); f=1.25*(rng.random((K,L))<p)-.25
        ft=f-f.mean(axis=1,keepdims=True)
        if center: ft=ft-ft.mean(axis=0,keepdims=True); ft=ft-ft.mean(axis=1,keepdims=True)
        U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*(ft**2).sum(axis=1)).sum()
        xh=abs(U)/np.sqrt(V); rej+= xh**2>crit; cov+= (np.interp(xh,xs,ucb)>=nu)
    return rej/R, cov/R
for nu in (0.0,2.0,3.708):
    print("nu=%.3f  rowonly: power=%.3f cov=%.3f | two-way: power=%.3f cov=%.3f (SE~%.3f)"%((nu,)+run(False,nu)+run(True,nu)+(0.008,)))
# T_pos sensitivity: lambda=K*sum q_j^2 = K*(L/2)*q^2 ; seed 300
def Tpos(rng,q):
    a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
    th=np.where(rng.random((K,L))<np.where(np.arange(L)<L//2,q,0.0)[None,:],0.0,th)
    return K*np.sum(abs(np.exp(1j*th).mean(axis=0))**2)
rng=np.random.default_rng(300); T0=np.array([Tpos(rng,0) for _ in range(2000)])
print("Tpos null mean %.2f sd %.2f (indep. expectation mean L=%d sd sqrt(L)=%.1f)"%(T0.mean(),T0.std(),L,np.sqrt(L)))
for q in (0.03,0.06,0.1):
    Tq=np.array([Tpos(rng,q) for _ in range(300)]); print("q=%.2f mean shift %.1f predicted K*(L/2)*q^2=%.1f"%(q,Tq.mean()-T0.mean(),K*(L/2)*q*q))
# analytic extrapolation to real scale (not simulated)
K0,L0=1e6,1e4; N0=K0*L0
print("real-scale: detect q (3sd of Tpos) ~ %.2e ; harmful q (0.1 sigma, d=0.1 avg, half positions) ~ %.2e"%(np.sqrt(3*np.sqrt(L0)/(K0*L0/2)), 0.1/(1.25*np.sqrt(8*K0/L0)*(L0/2)*0.1)))
