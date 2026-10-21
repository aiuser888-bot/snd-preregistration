# Factorial design: leakage structure x drift structure x centering (raw / row / two-way), paired on identical data.
# theta leakage: with prob q_bj the angle is replaced by 0  => E exp(i theta_bj) = q_bj  (angle bias).  L=100, K=400, R=4000.
import numpy as np
from scipy import stats
L,K,R=100,400,4000; j=np.arange(L); base=2*np.pi*j/L; crit=stats.chi2.ppf(.95,2)
g=(j/L-0.5)*0.4                       # position profile, zero mean over j, range [-0.2,0.2)
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def methods(f):
    row=f-f.mean(axis=1,keepdims=True)
    two=row-row.mean(axis=0,keepdims=True); two=two-two.mean(axis=1,keepdims=True)
    return {"raw":f,"row":row,"two":two}
# cell: (name, baseline shift sbar, position drift mode, leakage mode)
# drift: 'none' | 'common' (d_j=g_j for all blocks) | 'block' (d_bj=c_b g_j, c_b=+-1 random sign per block)
# leak : 'none' | 'uniform' (q everywhere) | 'static' (q_j=q where g_j>0, same every block) | 'blockvar' (q_bj=q where d_bj>0)
cells=[("E null control                        ",0.00,"none","none",0.0),
       ("A uniform leak x systematic baseline  ",0.03,"none","uniform",0.03),
       ("B static pos-leak x common drift      ",0.00,"common","static",0.03),
       ("D static pos-leak x block-specific drift",0.00,"block","static",0.03),
       ("C block-varying leak x block-specific drift",0.00,"block","blockvar",0.03),
       ("F leak only, no drift/shift (uniform)   ",0.00,"none","uniform",0.03)]
print("L=%d K=%d R=%d  sigma_U=sqrt(V)~%.1f ; type-I with chi2_2 reference (nominal .05, SE %.4f)"%(L,K,R,np.sqrt(K*L/8),np.sqrt(.05*.95/R)))
for ci,(name,sbar,dm,lk,q) in enumerate(cells):
    rng=np.random.default_rng([20261016,ci]); rej={m:0 for m in("raw","row","two")}; reU={m:[] for m in rej}; Vs=[]
    for _ in range(R):
        c=rng.choice([-1.0,1.0],(K,1)); d=(0*g[None,:]+0) if dm=="none" else (np.broadcast_to(g,(K,L)) if dm=="common" else c*g[None,:])
        p=np.clip(0.2+sbar+d,.01,.99); f=1.25*(rng.random((K,L))<p)-.25
        a=rng.uniform(0,2*np.pi,(K,1)); th=np.take_along_axis(a+base,np.argsort(rng.random((K,L)),axis=1),axis=1)
        if lk=="uniform": qm=np.full((K,L),q)
        elif lk=="static": qm=np.broadcast_to(np.where(g>0,q,0.0),(K,L))
        elif lk=="blockvar": qm=np.where(d>0,q,0.0)
        else: qm=np.zeros((K,L))
        th=np.where(rng.random((K,L))<qm,0.0,th)
        for m,ft in methods(f).items():
            U=(ft*np.exp(1j*th)).sum(); V=(L/(2*(L-1))*((ft-ft.mean(axis=1,keepdims=True))**2).sum(axis=1)).sum()
            rej[m]+=abs(U)**2/V>crit; reU[m].append(U.real/np.sqrt(V))
    out="  ".join("%s %.3f[%.3f,%.3f]"%(m,rej[m]/R,*wilson(rej[m],R)) for m in("raw","row","two"))
    bias="  ".join("%s %+.2f"%(m,np.mean(reU[m])) for m in("raw","row","two"))
    print("%s\n   type-I: %s\n   mean Re(U)/sqrt(V) (bias in sigma units): %s"%(name,out,bias))
# analytic bias predictions (raw): B: K*sum_j 1.25*g_j*q_j / sigma ; A: q*K*L*1.25*sbar/sigma
sig=np.sqrt(K*L*0.25/2)
print("\nanalytic raw bias/sigma:  A=%.2f  B=%.2f"%(0.03*K*L*1.25*0.03/sig, K*np.sum(1.25*g*np.where(g>0,0.03,0.0))/sig))
