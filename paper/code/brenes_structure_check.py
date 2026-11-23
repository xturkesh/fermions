"""Brenes et al. states: (i) full gauge-field simulation vs sector average; (ii) H_ebar = H_clean + g^2 sum_k A_k rho_k + const."""
import itertools, numpy as np
from scipy.linalg import expm
L=4; w=1.0; m=0.3; g2=1.7; theta=0.2; e0=0
odd=lambda j: j%2==1
vac=tuple(1 if odd(j) else 0 for j in range(1,L+1)); Nf=sum(vac)
def rho(n): return np.array([n[j-1]-(1 if odd(j) else 0) for j in range(1,L+1)])
# ---------- (ii) algebraic identity in every sector
confs=[n for n in itertools.product([0,1],repeat=L) if sum(n)==Nf]
F=np.array([np.cumsum(rho(n)) for n in confs]); R=np.array([rho(n) for n in confs])
mass=np.array([sum(m*(-1)**j*n[j-1] for j in range(1,L+1)) for n in confs])
err=0
for eb in itertools.product([-1,0,1],repeat=L):
    eb=np.array(eb); E=eb+F                       # fields in sector Q(eb): E_j = eb_j + F_j
    d=mass+0.5*g2*((E+theta)**2).sum(1)
    A=np.cumsum((eb+theta)[::-1])[::-1]           # A_k = sum_{j>=k}(eb_j+theta)
    dclean=mass+0.5*g2*(F**2).sum(1)
    d2=dclean+g2*(R*A).sum(1)+0.5*g2*((eb+theta)**2).sum()
    err=max(err,np.abs(d-d2).max())
print('identity H_ebar = H_clean + g^2 sum_k A_k rho_k + const: max error',err)
# ---------- (i) full LGT with explicit link fields (truncated |E|<=K), Brenes initial state
K=L+1
fvals=list(range(-K,K+1))
states=[(n,E) for n in itertools.product([0,1],repeat=L) if sum(n)==Nf for E in itertools.product(fvals,repeat=L)]
# keep only gauge-invariant states reachable from the initial superposition: G_j = E_j-E_{j-1}-rho_j = eb_j-eb_{j-1}
idx={s:i for i,s in enumerate(states)}
D=len(states); print('full space dim',D)
import scipy.sparse as sp
rows=[];cols=[];vals=[];diag=np.zeros(D)
for i,(n,E) in enumerate(states):
    diag[i]=sum(m*(-1)**j*n[j-1] for j in range(1,L+1))+0.5*g2*sum((x+theta)**2 for x in E)
    for j in range(1,L):
        if n[j]==1 and n[j-1]==0:                  # f_j^dag U_j f_{j+1}: E_j -> E_j+1
            nn=list(n); nn[j]=0; nn[j-1]=1; EE=list(E); EE[j-1]+=1
            s2=(tuple(nn),tuple(EE))
            if s2 in idx:
                k=idx[s2]; rows+= [k,i]; cols+=[i,k]; vals+=[-w,-w]
H=sp.coo_matrix((vals,(rows,cols)),shape=(D,D)).tocsr()+sp.diags(diag)
psi0=np.zeros(D,complex)
for eb in itertools.product([-1,0,1],repeat=L):
    psi0[idx[(vac,eb)]]=3**(-L/2)
from scipy.sparse.linalg import expm_multiply
npa_full=[];ts=[0.5,2.0,7.0]
for t in ts:
    psi=expm_multiply(-1j*t*H,psi0)
    P=np.abs(psi)**2
    npa=sum(P[i]*np.abs(rho(n)).sum()/L for i,(n,E) in enumerate(states))
    npa_full.append(npa)
# sector average
Kh=np.zeros((len(confs),len(confs))); ci={n:i for i,n in enumerate(confs)}
for i,n in enumerate(confs):
    for j in range(1,L):
        if n[j]==1 and n[j-1]==0:
            nn=list(n); nn[j]=0; nn[j-1]=1; k=ci[tuple(nn)]; Kh[k,i]+=-w; Kh[i,k]+=-w
npa_sec=np.zeros(len(ts)); v=ci[vac]
for eb in itertools.product([-1,0,1],repeat=L):
    E=np.array(eb)+F; d=mass+0.5*g2*((E+theta)**2).sum(1)
    Hs=Kh+np.diag(d)
    for a,t in enumerate(ts):
        ps=expm(-1j*t*Hs)[:,v]; npa_sec[a]+=3**(-L)*(np.abs(ps)**2@(np.abs(R).sum(1)/L))
print('particle density, full LGT vs sector average:',np.array(npa_full),npa_sec,'max diff',np.abs(np.array(npa_full)-npa_sec).max())
