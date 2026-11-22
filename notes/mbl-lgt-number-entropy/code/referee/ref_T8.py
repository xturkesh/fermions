import numpy as np, itertools
I2=np.eye(2); X=np.array([[0,1],[1,0]],float); Z=np.diag([1.,-1.])
def op(d,n):
    M=np.array([[1.]])
    for k in range(n): M=np.kron(M,d.get(k,I2))
    return M
ns=4; links=[(0,1),(2,3),(0,2),(1,3)]; nl=len(links); n=ns+nl; J=1.0; h=0.73
f=[op({**{k:Z for k in range(j)}, j:np.array([[0,1],[0,0]],float)},n) for j in range(ns)]
sz=[op({ns+l:Z},n) for l in range(nl)]; sx=[op({ns+l:X},n) for l in range(nl)]
D=2**n; H=np.zeros((D,D))
for l,(i,j) in enumerate(links): H+=-J*sz[l]@(f[i].T@f[j]+f[j].T@f[i])
for i in range(ns):
    P=np.eye(D)
    for l,(a,b) in enumerate(links):
        if i in (a,b): P=P@sx[l]
    H+=-h*P
Eed=np.sort(np.linalg.eigvalsh(H))
spec=[]
for q in itertools.product([1,-1],repeat=ns):
    for W in [1,-1]:
        s={(0,1):1.,(2,3):W,(0,2):1.,(1,3):1.}
        hq=np.diag(2*h*np.array(q,float))
        for (i,j),v in s.items(): hq[i,j]=hq[j,i]=-J*v
        e=np.linalg.eigvalsh(hq)
        for occ in itertools.product([0,1],repeat=ns):
            if (-1)**sum(occ)!=np.prod(q): continue
            spec.append(np.dot(occ,e)-h*sum(q))
spec=np.sort(spec); print(len(spec),D,"max|dE|=",np.abs(spec-Eed).max())
