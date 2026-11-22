import numpy as np
from scipy.special import gammaln
from scipy.linalg import eigh_tridiagonal
# (a) T4 power-series coefficients d_k = 1/(2k(2k-1)) - ln2 * c_k, c_k from sqrt(1-s^2)=1-sum c_k s^{2k}
k=np.arange(1,200001,dtype=float)
logc=gammaln(2*k+1)-np.log(2*k-1)-2*gammaln(k+1)-k*np.log(4)
d=1/(2*k*(2*k-1))-np.log(2)*np.exp(logc)
print("d_1=",d[0]," max d_k (k>=2)=",d[1:].max()," at k=",np.argmax(d[1:])+2, " partial sum to 2e5:",d.sum())
# (b) T3 with negative h: bound h^2/(N(J+h)^2) vs exact
J=1.0; rng=np.random.default_rng(0)
for h in [-0.9,-0.5,0.5]:
    N=40; acc=0; S=400
    for s in range(S):
        q=rng.choice([-1,1],N); E,V=eigh_tridiagonal(2*h*q.astype(float),-J*np.ones(N-1))
        st=(-1.0)**np.arange(1,N+1); acc+=((st@V**2)**2).sum()/N/S
    print(f"h={h}: exact E Ibar={acc:.4f}, note's global bound h^2/(N(J+h)^2)={h*h/(N*(J+h)**2):.4f}, fixed bound with |h|={h*h/(N*(J+abs(h))**2):.4f}")
# (c) dimer model: critical energy E_c=2h; rotation of (-,-) dimer at E_c: trace = 16 h^2 - 2 (J=1)
for h in [0.2,0.25,0.3,1/(2*np.sqrt(2)),np.sqrt(3)/4,0.45]:
    tr=16*h*h-2; print(f"h={h:.4f}: trace T_(--)^2(E_c)={tr:+.4f}, rotation angle/pi={np.arccos(tr/2)/np.pi:.4f}")
# Lyapunov exponent of the dimer model near E_c (transfer-matrix products)
def lyap(h,E,L=400000,seed=1):
    r=np.random.default_rng(seed); v=2*h*r.choice([-1,1],L//2).repeat(2)
    x=np.array([1.0,0.3]); g=0.0
    for vn in v:
        x=np.array([(vn-E)*x[0]-x[1],x[0]]); nr=np.hypot(*x); g+=np.log(nr); x/=nr
    return g/L
for h in [0.2,0.25,1/(2*np.sqrt(2))]:
    Ec=2*h
    out=[(eps,lyap(h,Ec+eps)) for eps in [0.0,0.02,0.04,0.08]]
    print(f"h={h:.4f}: "+"  ".join(f"eps={e}: gamma={g:.2e} ratio g/eps^2={g/e**2 if e>0 else float('nan'):.3f}" for e,g in out))
