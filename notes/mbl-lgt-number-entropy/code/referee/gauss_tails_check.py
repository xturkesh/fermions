import numpy as np
from scipy.linalg import expm
g=0.3
for alpha in [1.0,1.5,2.0]:
    for L in [64,256]:
        x=np.arange(L); R=np.abs(x[:,None]-x[None,:]).astype(float)
        th=np.where(R>0,g*np.where(R>0,R,1.0)**(-alpha),0.0)
        u=expm(1j*th); n=(x%2==0).astype(float); A=x<L//2
        M=u.conj().T@np.diag(A.astype(float))@u          # U^dag N_A U single-particle
        C0=np.diag(n)
        var=np.real(np.trace(C0@M@(np.eye(L)-C0)@M))
        K=M-np.diag(A.astype(float)); ev=np.linalg.eigvalsh((K+K.conj().T)/2)
        # Hankel block of theta: trace norm vs trace (positivity argument)
        B=th[np.ix_(A,~A)]; sv=np.linalg.svd(B,compute_uv=False)
        print(f"alpha={alpha} L={L} Var={var:.4f} ||Delta||={ev[ev>0].sum():.4f}  ||theta_AB||_1={sv.sum():.4f}  first moment={B.sum():.3f}")
