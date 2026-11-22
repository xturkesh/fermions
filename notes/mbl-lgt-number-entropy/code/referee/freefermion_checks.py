import numpy as np
rng=np.random.default_rng(0)
def pb_dist(nu):
    p=np.array([1.0])
    for v in nu: p=np.convolve(p,[1-v,v])
    return p
def H(p):
    p=p[p>1e-300]; return -(p*np.log(p)).sum()
# 1. Counterexample to "S_N ~ 1/2 ln S" two-sided for free fermions: M modes with nu=1/M
print("M, S, Var, S_N, 0.5 ln S, S2")
for M in [10,100,1000,10000,100000]:
    nu=np.full(M,1.0/M)
    S=np.sum(-nu*np.log(nu)-(1-nu)*np.log(1-nu)); Var=np.sum(nu*(1-nu))
    S2=-np.sum(np.log(nu**2+(1-nu)**2))
    p=pb_dist(nu) if M<=10000 else None
    SN=H(p) if p is not None else float('nan')
    print(M, round(S,3), round(Var,4), round(SN,4), round(0.5*np.log(S),3), round(S2,4))
# 2. random checks of Thm 5.1 inequalities
worst=[-1e9]*5
for _ in range(3000):
    k=rng.integers(1,60); nu=rng.uniform(0,1,k)**rng.uniform(0.2,5)
    p=pb_dist(nu); Var=np.sum(nu*(1-nu)); S=np.sum([-v*np.log(v)-(1-v)*np.log(1-v) for v in nu if 0<v<1])
    S2=-np.sum(np.log(nu**2+(1-nu)**2)); SN=H(p); Hinf=-np.log(p.max())
    worst[0]=max(worst[0], 4*np.log(2)*Var-S)
    worst[1]=max(worst[1], 2*Var-S2); worst[2]=max(worst[2], S2-4*np.log(2)*Var)
    if Var>0: worst[3]=max(worst[3], 0.5*np.log(8*Var/np.pi)-Hinf)
    worst[4]=max(worst[4], SN-0.5*np.log(2*np.pi*np.e*(Var+1/12)))
print("max violations (should be <=0):",worst)
# 3. Prop 5.2 numbers
from scipy.stats import binom
for n,q in [(10000,0.5),(10000,0.05)]:
    k=np.arange(n+1); p=binom.pmf(k,n,q); Var=n*q*(1-q)
    print("binom",n,q,"Var",Var,"H-1/2ln(2pieVar)=",H(p)-0.5*np.log(2*np.pi*np.e*Var))
