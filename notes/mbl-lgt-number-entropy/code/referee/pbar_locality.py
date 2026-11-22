# Independent check of Prop 4.1: dephased number distribution depends only on gates in light cone R
import numpy as np
from itertools import combinations
def basis(L,N):
    st=[]
    for occ in combinations(range(L),N):
        s=0
        for i in occ: s|=1<<i
        st.append(s)
    return np.array(sorted(st))
def rand_gate(rng,f):
    a,e=rng.normal(size=2); M=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); M=(M+M.conj().T)/2
    w,v=np.linalg.eigh(M); G=(v*np.exp(-1j*f*w))@v.conj().T
    return np.exp(-1j*f*a),np.exp(-1j*f*e),G
def build_U(L,gates,du):
    N=L//2; st=basis(L,N); idx={s:k for k,s in enumerate(st)}; D=len(st); l=L//2
    U=np.eye(D,dtype=complex)
    for layer in range(2*du):
        for i in range(L-1):
            r=i-(l-1)
            if (r-layer)%2: continue
            g00,g11,G=gates[(layer,r)]
            Gm=np.zeros((D,D),complex)
            for k,s in enumerate(st):
                bi=(s>>i)&1; bj=(s>>(i+1))&1
                if bi==0 and bj==0: Gm[k,k]=g00
                elif bi==1 and bj==1: Gm[k,k]=g11
                else:
                    s01 = s if (bi==0) else s^((1<<i)|(1<<(i+1)))   # n_i=0,n_{i+1}=1
                    s10 = s01^((1<<i)|(1<<(i+1)))
                    a=0 if bi==0 else 1
                    Gm[idx[s01],k]+=G[0,a]; Gm[idx[s10],k]+=G[1,a]
            U=Gm@U
    NA=np.array([bin(s&((1<<l)-1)).count('1') for s in st])
    return st,idx,U,NA
def pbar(L,gates,du,par):
    st,idx,U,NA=build_U(L,gates,du); l=L//2
    # Neel pattern fixed relative to the cut: site l-1+r occupied iff r%2==par
    s0=sum(1<<i for i in range(L) if (i-(l-1))%2==par)
    psi=np.zeros(len(st),complex); psi[idx[s0]]=1; n0=NA[idx[s0]]
    assert np.allclose(U.conj().T@U,np.eye(len(st)))
    phi=U.conj().T@psi
    P=np.zeros((l+1,len(st)))
    for n in range(l+1): P[n]=(np.abs(U[NA==n,:])**2).sum(0)   # <c|U^dag P_n U|c>
    pb=P@(np.abs(phi)**2)
    return np.array([pb[n0+k] if 0<=n0+k<=l else 0 for k in range(-4,5)]), U, phi, NA, psi
rng=np.random.default_rng(7)
for du,Ls in [(1,[4,6,8,10]),(2,[8,10,12,14])]:
    gates={(lay,r):rand_gate(rng,0.9) for lay in range(2*du) for r in range(-12,13)}
    for par in [0,1]:
        ref=None
        for L in Ls:
            v,_,_,_,_=pbar(L,gates,du,par)
            d=None if ref is None else np.abs(v-ref).max()
            print(f"du={du} par={par} L={L} pbar(n0-4..n0+4)={np.round(v,8)}  maxdiff vs prev={d}")
            ref=v
# brute-force time average vs formula for a random diagonal D (L=8,du=1)
gates={(lay,r):rand_gate(rng,0.9) for lay in range(2) for r in range(-12,13)}
v,U,phi,NA,psi=pbar(8,gates,1,0)
st=basis(8,4); Dc=rng.normal(size=len(st))
ts=rng.uniform(0,1e6,4000); acc=np.zeros(5)
l=4; n0=NA[np.argmax(np.abs(psi))]
for t in ts:
    ps=U@(np.exp(-1j*Dc*t)*phi)
    p=np.bincount(NA,weights=np.abs(ps)**2,minlength=l+1); acc+=p
acc/=len(ts)
print("time-avg p (n=0..4):",np.round(acc,4)); print("formula      :",np.round([v[4+k-n0] for k in range(5)],4))
