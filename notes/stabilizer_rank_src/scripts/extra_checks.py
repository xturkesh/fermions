import numpy as np, itertools as it
from functools import reduce
kron = lambda vs: reduce(np.kron, vs)
c, s = np.cos(np.pi/8), np.sin(np.pi/8); H1 = np.array([c, s])
phi = np.array([1,0,0,1])/np.sqrt(2); g = np.array([1,1,1,-1])/2
print("H^2 == (Phi+ + CZ|++>)/sqrt2 :", np.allclose(kron([H1,H1]), (phi+g)/np.sqrt(2)), " <Phi+|CZ++> =", phi@g)
# stabilizer states for n=4 via BFS, then F(W_4) (claim 9/16)
n=4; N=2**n; I2=np.eye(2); Hd=np.array([[1,1],[1,-1]])/np.sqrt(2); Sg=np.diag([1,1j])
on=lambda G,q: kron([G if i==q else I2 for i in range(n)])
def cnot(a,b):
    U=np.zeros((N,N))
    for x in range(N):
        bits=[(x>>(n-1-i))&1 for i in range(n)]
        if bits[a]: bits[b]^=1
        U[int("".join(map(str,bits)),2),x]=1
    return U
gens=[on(Hd,q) for q in range(n)]+[on(Sg,q) for q in range(n)]+[cnot(a,b) for a in range(n) for b in range(n) if a!=b]
def key(v):
    k=np.argmax(np.abs(v)>1e-9); v=v*np.conj(v[k])/abs(v[k]); return tuple(np.round(v,6).view(float))
st=np.zeros(N,complex); st[0]=1; seen={key(st):st}; fr=[st]
while fr:
    new=[]
    for v in fr:
        for G in gens:
            w=G@v; kk=key(w)
            if kk not in seen: seen[kk]=w; new.append(w)
    fr=new
S=np.array(list(seen.values())); print("n=4 stabilizer states:", len(S))
W=np.zeros(N); 
for i in range(n): W[1<<i]=1
W/=np.linalg.norm(W); print("F(W_4) = %.6f  (9/16 = %.6f)" % (np.max(np.abs(S.conj()@W)**2), 9/16))
