"""Independent referee check of T1 (duality incl. dictionary), C1, P5(1), P5(3), T4(2)
1D Z2 LGT with boundary link: sites 1..N, links b_1..b_N (b_N boundary), sigma^x_{b_0}:=1.
Qubit order: fermion JW qubits 0..N-1, then link qubits N..2N-1."""
import numpy as np, itertools
from scipy.linalg import expm, eigh
I2=np.eye(2); X=np.array([[0,1],[1,0]],float); Z=np.diag([1.,-1.]); a=np.array([[0,1],[0,0]],float) # a|1>=|0>, occupied=|1>
Hd=np.array([[1,1],[1,-1]])/np.sqrt(2)
def op(lst,n):
    M=np.array([[1.]])
    for k in range(n): M=np.kron(M,lst.get(k,I2))
    return M
def fermion(j,n):  # JW annihilation, occupied = bit 1 ; Z on earlier qubits
    d={k:Z for k in range(j)}; d[j]=np.array([[0,1],[0,0]],float); return op(d,n)
def run(N,J,h,ell,times):
    n=2*N; D=2**n
    f=[fermion(j,n) for j in range(N)]
    nop=[fk.T@fk for fk in f]
    sz=[op({N+j:Z},n) for j in range(N)]; sx=[op({N+j:X},n) for j in range(N)]
    one=np.eye(D)
    sxb=lambda j: one if j<0 else sx[j]
    H=np.zeros((D,D))
    for j in range(N-1): H+=-J*sz[j]@(f[j].T@f[j+1]+f[j+1].T@f[j])
    for j in range(N): H+=-h*sxb(j-1)@sxb(j)
    par=[one-2*nop[j] for j in range(N)]
    G=[sxb(j-1)@par[j]@sx[j] for j in range(N)]
    print(f"N={N}: max|[H,G]|=",max(np.abs(H@g-g@H).max() for g in G))
    # explicit U: dual basis (n bits, q bits with bit0<->q=+1) -> (n bits, x bits) then Hadamard on links
    P=np.zeros((D,D))
    for nb in itertools.product([0,1],repeat=N):
        for qb in itertools.product([0,1],repeat=N):
            q=[1-2*b for b in qb]; x=[]; s=1
            for j in range(N):
                s*=q[j]*(-1)**nb[j]; x.append(s)
            xb=[(1-v)//2 for v in x]
            i_in=int(''.join(map(str,list(nb)+list(qb))),2); i_out=int(''.join(map(str,list(nb)+xb)),2)
            P[i_out,i_in]=1
    U=op({N+j:Hd for j in range(N)},n)@P
    print(" unitary err", np.abs(U.T@U-one).max())
    # dual operators: c_j = JW on fermion qubits (same matrices), tau^x_j = Z on qubit N+j, tau^z_j = X on qubit N+j
    tx=[op({N+j:Z},n) for j in range(N)]; tz=[op({N+j:X},n) for j in range(N)]
    Hdual=np.zeros((D,D))
    for j in range(N-1): Hdual+=-J*(f[j].T@f[j+1]+f[j+1].T@f[j])
    for j in range(N): Hdual+=2*h*tx[j]@nop[j]-h*tx[j]
    print(" |U^T H U - Hdual| =",np.abs(U.T@H@U-Hdual).max())
    print(" |U^T G_j U - tau^x_j| =",max(np.abs(U.T@G[j]@U-tx[j]).max() for j in range(N)))
    print(" |U^T f_j U - c_j tau^z_j| =",max(np.abs(U.T@f[j]@U-f[j]@tz[j]).max() for j in range(N)))
    tzN=lambda j: one if j>=N else tz[j]
    print(" |U^T sz_j U - tz_j tz_j+1| =",max(np.abs(U.T@sz[j]@U-tz[j]@tzN(j+1)).max() for j in range(N)))
    # C1: Psi0 = |up..up>_z (link bit 0) (x) CDW (occupied on even sites 2,4,.. -> index 1,3,..)
    nb=[1 if (j+1)%2==0 else 0 for j in range(N)]
    idx=int(''.join(map(str,nb+[0]*N)),2); Psi0=np.zeros(D); Psi0[idx]=1
    dual=U.T@Psi0
    exp_dual=np.zeros(D)
    for qb in itertools.product([0,1],repeat=N):
        exp_dual[int(''.join(map(str,nb+list(qb))),2)]=2**(-N/2)
    print(" C1 err",np.abs(dual-exp_dual).max())
    # entropies
    E,V=eigh(H)
    def S_of(psi,keep):
        psi=psi.reshape([2]*n); rest=[k for k in range(n) if k not in keep]
        M=np.transpose(psi,keep+rest).reshape(2**len(keep),-1)
        s=np.linalg.svd(M,compute_uv=False)**2; s=s[s>1e-15]; return -(s*np.log(s)).sum()
    A=list(range(ell))+[N+j for j in range(ell)]
    NA=sum(nop[j] for j in range(ell)).diagonal()
    for t in times:
        psi=V@(np.exp(-1j*E*t)*(V.T@Psi0))
        Sf=S_of(psi,A); Sd=S_of(U.T@psi,A)
        pN=np.bincount(NA.astype(int),weights=np.abs(psi)**2); pN=pN[pN>1e-15]; SN=-(pN*np.log(pN)).sum()
        # S_G : sector-wise entanglement of dual state with q fixed (fermions 1..ell)
        dpsi=(U.T@psi).reshape(2**N,2**N)  # rows: fermion config, cols: q config
        SG=0.0
        for col in range(2**N):
            v=dpsi[:,col]*2**(N/2)
            M=v.reshape(2**ell,-1); s=np.linalg.svd(M,compute_uv=False)**2; s=s[s>1e-15]; SG+=-(s*np.log(s)).sum()/2**N
        print(f"  t={t:7.2f} S_full={Sf:.4f} S_dual={Sd:.4f} diff={Sf-Sd:+.4f} (<=ln2={np.log(2):.3f})  S_G={SG:.4f}  S_N={SN:.4f}  S_G<S_N? {SG<SN}")
    # T4(2): definite-charge state |x>_sigma (x) CDW, x random
    rng=np.random.default_rng(1)
    xb=list(rng.integers(0,2,N))
    idx=int(''.join(map(str,nb+xb)),2); e=np.zeros(D); e[idx]=1
    Psi1=op({N+j:Hd for j in range(N)},n)@e
    for t in [0.7,3.1]:
        psi=V@(np.exp(-1j*E*t)*(V.T@Psi1))
        print(f"  T4(2) t={t}: S_full={S_of(psi,A):.6f}  S_dual={S_of(U.T@psi,A):.6f}  S_fermions-only dual={S_of(U.T@psi,list(range(ell))):.6f}")
for N,ell in [(4,2),(6,3)]:
    run(N,1.0,1.0,ell,[0.3,1.0,2.0,5.0,20.0])
