"""r23 (V2.0.10): exact checks of Lemma 5 - design moments, beta <= sqrt(D) (tight for aligned tilt), TV bound, GF(2) hash identity. Superseded by code/tests/test_snd.py::test_lemma5_*."""
import numpy as np, itertools, math
rng=np.random.default_rng(20261021)
# --- (1) exact design moments: E|U_b|^2 = L/(L-1) sum f^2 (=2V_b), E U_b^2 = 0, E Re(w^* U)^2 = V
for L in (4,5,6):
    f=rng.normal(size=L); f-=f.mean()
    A=7; alphas=2*np.pi*np.arange(A)/A
    perms=list(itertools.permutations(range(L)))
    U=np.array([[ (f*np.exp(1j*(a+2*np.pi*np.array(p)/L))).sum() for p in perms] for a in alphas]).ravel()
    Vb=L/(2*(L-1))*(f**2).sum()
    print("L",L,"E|U|^2/(2V)=%.12f  |E U|=%.2e  |E U^2|=%.2e  E Re(wU)^2/V=%.12f"%(
        (abs(U)**2).mean()/(2*Vb),abs(U.mean()),abs((U**2).mean()), (( (np.exp(-0.7j)*U).real**2).mean())/Vb))
# --- (2) discrete seed space K=2,L=4, A=6 : beta <= sqrt(D), tightness for aligned tilt
L=4; K=2; A=6
alphas=2*np.pi*np.arange(A)/A; perms=list(itertools.permutations(range(L)))
def Ub_table(f):
    return np.array([[ (f*np.exp(1j*(a+2*np.pi*np.array(p)/L))).sum() for p in perms] for a in alphas]).ravel()
fs=[]; tabs=[]
for b in range(K):
    f=rng.normal(size=L); f-=f.mean(); fs.append(f); tabs.append(Ub_table(f))
V=sum(L/(2*(L-1))*(f**2).sum() for f in fs)
U=(tabs[0][:,None]+tabs[1][None,:]).ravel()
print("states",U.size,"E U=%.1e E|U|^2/(2V)=%.12f"%(abs(U.mean()),(abs(U)**2).mean()/(2*V)))
def beta_D(l):
    B=(l*U).mean(); D=((l-1)**2).mean(); return abs(B)/np.sqrt(V), np.sqrt(D), ((l-1)**2).mean(), 0.5*np.abs(l-1).mean()
w=np.exp(1j*0.9)
for eps in (0.05,0.2,0.4):
    l=1+eps*(np.conj(w)*U).real/np.sqrt(V)
    assert l.min()>=0
    b,sD,D,tv=beta_D(l); print("aligned eps=%.2f beta=%.9f sqrtD=%.9f ratio=%.9f  TV=%.5f <= .5sqrtD=%.5f"%(eps,b,sD,b/sD,tv,0.5*sD))
worst=0
for t in range(300):
    g=rng.normal(size=U.size)*rng.uniform(0.01,1.5); l=np.exp(g); l/=l.mean()
    b,sD,D,tv=beta_D(l); worst=max(worst,b/sD); assert tv<=0.5*sD+1e-12
    # random event A: P_G(A)-P_0(A) <= TV
    Aset=rng.random(U.size)<rng.uniform(0.05,0.95)
    assert (l*Aset).mean()-Aset.mean()<=tv+1e-12
print("random tilts max beta/sqrtD = %.4f (<=1)"%worst)
# --- (3) GF(2) hash: E_A D_out = (2^m-1) CP(R)   (exhaustive, n=6, m=3)
n,m=6,3
p1=0.5+0.2; bits=((np.arange(2**n)[:,None]>>np.arange(n))&1)
pr=np.prod(np.where(bits==1,p1,1-p1),axis=1); CP=(pr**2).sum()
par=np.array([[bin(x&r).count('1')&1 for r in range(2**n)] for x in range(2**n)])  # x,r
tot=0.0; cnt=0
R=np.arange(2**n)
for r1 in range(2**n):
    for r2 in range(2**n):
        z12=par[:,r1]+2*par[:,r2]            # (64,)
        z=z12[None,:]+4*par[:,R].T           # (r3,x)
        pm=np.stack([((z==k)*pr[None,:]).sum(1) for k in range(2**m)],axis=1)
        tot+=(2**m*(pm**2).sum(1)-1).sum(); cnt+=2**n
print("E_A D_out=%.12f  (2^m-1)CP=%.12f  2^(m-H2)=%.6f"%(tot/cnt,(2**m-1)*CP,2**m*CP))
