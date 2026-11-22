import numpy as np, itertools
from lgt1d_ed import *
from scipy.linalg import expm
N=4; J=1.; h=0.7; occ=[0,2]; t=2.0
m = build(N,J,h,'open_bl')
vf = np.zeros(2**N); idx=int(''.join('1' if j in occ else '0' for j in range(N)),2); vf[idx]=1
vs = np.zeros(2**N); vs[0]=1
psi0=np.kron(vf,vs).astype(complex)
psi = expm(-1j*t*m['H'].toarray())@psi0
ED_n=np.array([np.vdot(psi,m['n'][x]@psi).real for x in range(N)])
Pn=np.zeros(N)
for q in itertools.product([1,-1],repeat=N):
    U=expm(-1j*t*hq(N,J,h,q)); Phi=U[:,occ]; Pn+=np.real(np.diag(Phi@Phi.conj().T))/2**N
print(ED_n, Pn)
# check Q expectation distribution
for q in list(itertools.product([1,-1],repeat=N))[:4]:
    P=np.eye(2**(2*N))
    for j in range(N): P=P@(np.eye(2**(2*N))+q[j]*m['Q'][j].toarray())/2
    print(q, np.vdot(psi0,P@psi0).real)
