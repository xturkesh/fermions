# validate: entropies of region A computed from the full ED state (sigma^z spin basis) vs ent_small.run
import numpy as np
from scipy.linalg import expm
from lgt1d_ed import build
from ent_small import run, entropy_svd, entropy_from_rho
N=6; J=1.; h=1.0; l=3; t=2.7
m=build(N,J,h,'open_bl')
occ=[0,2,4]
vf=np.zeros(2**N); vf[int(''.join('1' if j in occ else '0' for j in range(N)),2)]=1
vs=np.zeros(2**N); vs[0]=1
psi=expm(-1j*t*m['H'].toarray())@np.kron(vf,vs).astype(complex)
P=psi.reshape(2**l,2**(N-l),2**l,2**(N-l)).transpose(0,2,1,3).reshape(4**l,-1)
Sfull_ED=entropy_svd(P)
Mf=psi.reshape(2**l,-1); Sf_ED=entropy_from_rho(Mf@Mf.conj().T)
r=run(N,J,h,[t])[0]
print('S_full ED =',Sfull_ED,' dual-reconstructed =',r[1])
print('S_f    ED =',Sf_ED,' dual-reconstructed =',r[5])
