import numpy as np, math
from scipy.linalg import expm
rng=np.random.default_rng(5)
L=10; l=5; J=1.0; Delta=1.0
dim=1<<L; occ=((np.arange(dim)[:,None]>>np.arange(L))&1)
def run(W,nmax=3):
    h=rng.uniform(-W,W,L)
    E=occ@h + J*Delta*(occ[:,:-1]*occ[:,1:]).sum(1)
    T=np.zeros((dim,dim))
    for i in range(L-1):
        for c in range(dim):
            if occ[c,i]!=occ[c,i+1]:
                T[c^((1<<i)|(1<<(i+1))),c]=1
    NA=np.diag(occ[:,:l].sum(1).astype(float))
    H=np.diag(E)+J*T
    Y=[NA]; dE=E[:,None]-E[None,:]
    gam_used=np.inf; out=[]
    for k in range(1,nmax+1):
        R=J*(T@Y[-1]-Y[-1]@T)
        assert np.abs(np.diag(R)).max()<1e-12
        mask=np.abs(R)>1e-14
        gam_used=min(gam_used,np.abs(dE[mask]).min())
        Yk=np.zeros_like(R); Yk[mask]=-R[mask]/dE[mask]
        Y.append(Yk)
        Nt=sum(Y); comm=H@Nt-Nt@H
        lhs=np.linalg.norm(comm,2); rhs=np.linalg.norm(J*(T@Yk-Yk@T),2)
        dn=np.linalg.norm(Nt-NA,2)
        rowsum=np.abs(Yk).sum(1).max()
        # support check: sites where Yk acts nontrivially
        sup=[]
        for i in range(L):
            # Yk commutes with sigma^x_i ?  check flip invariance of matrix elements
            P=np.arange(dim)^(1<<i)
            if np.abs(Yk[np.ix_(P,P)]-Yk).max()>1e-12: sup.append(i)
        bound_Y=(J/gam_used)*np.prod([2*(4*j+1)*J/gam_used for j in range(1,k)])
        out.append((k,dn,lhs,abs(lhs-rhs),rowsum,bound_Y,sup,math.factorial(k)*(10*J/gam_used)**k*J,gam_used))
    return out,H,NA,Y
for W in [50,200]:
    out,H,NA,Y=run(W)
    for o in out:
        k,dn,lhs,err,rs,bY,sup,epsb,g=o
        print(f"W={W} n={k}: ||Ntilde-N_A||={dn:.3e}  ||[H,Nt]||={lhs:.3e} (=J[T,Y_n] err {err:.1e})  |||Y_n|||={rs:.3e} <= bound {bY:.3e}; eps bound {epsb:.3e}; gamma_min={g:.3f}; supp(Y_n) sites={sup}")
    # dynamics check for n=2
    Nt=sum(Y[:3]); d=np.linalg.norm(Nt-NA,2); e=np.linalg.norm(H@Nt-Nt@H,2)
    viol=0
    for t in [1,10,100,1000]:
        Ut=expm(-1j*H*t); NAt=Ut.conj().T@NA@Ut
        lhs=np.linalg.norm(NAt-NA,2); 
        if lhs>2*d+e*t+1e-9: viol+=1
        print(f"   t={t}: ||N_A(t)-N_A||={lhs:.3e} <= 2d+eps t = {2*d+e*t:.3e}")
