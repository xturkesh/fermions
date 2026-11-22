"""T3 check on the user's model: free fermions, hoppings J_i = u_i^delta (u~U(0,1)),
zero potential, periodic ring (as in /home/user/fermions with bc=.true.), initial
state = the code's density wave (sites 4,8,...,L occupied, 1-indexed: quarter filling),
region A = first L/2 sites.  For each t: S (von Neumann), Var N_A, exact S_N from the
Poisson-binomial law of N_A, and the T3 bounds
   1/2 ln(8 Var/pi) <= S_N <= 1/2 ln(2 pi e (Var+1/12)) <= 1/2 ln(2 pi e (S/(4 ln 2)+1/12)).
usage: python3 random_hopping_T3.py L delta nreal
"""
import sys
import numpy as np

L, delta, nreal = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
ell = L//2
tgrid = 10.0**np.arange(-1, 12.01, 0.5)
rng = np.random.default_rng(int(1000*delta) + L)
occ = np.arange(3, L, 4)            # 0-indexed sites 3,7,11,... == 1-indexed 4,8,12,...

def hbin(nu):
    nu = np.clip(nu, 1e-15, 1-1e-15)
    return -nu*np.log(nu) - (1-nu)*np.log(1-nu)

def pb(nus):
    p = np.array([1.0])
    for v in nus:
        p = np.convolve(p, [1-v, v])
    return p

S = np.zeros((nreal, len(tgrid))); V = np.zeros_like(S); SN = np.zeros_like(S)
viol_up = viol_lo = viol_S = 0.0
for r in range(nreal):
    J = rng.random(L)**delta
    h = np.zeros((L, L))
    for i in range(L):
        h[i, (i+1) % L] = h[(i+1) % L, i] = J[i]
    E, W = np.linalg.eigh(h)
    Wocc = W[occ, :].conj()                     # rows: occupied sites
    for k, t in enumerate(tgrid):
        Phi = (W*np.exp(-1j*E*t)) @ Wocc.T      # U(t)[:, occ]
        CA = Phi[:ell] @ Phi[:ell].conj().T
        nu = np.clip(np.linalg.eigvalsh(CA), 0, 1)
        s = hbin(nu).sum(); var = (nu*(1-nu)).sum()
        p = pb(nu); p = p[p > 1e-300]
        sn = -(p*np.log(p)).sum()
        S[r, k], V[r, k], SN[r, k] = s, var, sn
        viol_up = max(viol_up, sn - 0.5*np.log(2*np.pi*np.e*(var + 1/12)))
        viol_S = max(viol_S, 0.5*np.log(2*np.pi*np.e*(var + 1/12)) - 0.5*np.log(2*np.pi*np.e*(s/(4*np.log(2)) + 1/12)))
        if var > 0:
            viol_lo = max(viol_lo, 0.5*np.log(8*var/np.pi) - sn)
print(f"L={L} delta={delta} nreal={nreal}; max violations (all must be <=0): "
      f"upper {viol_up:.2e}, S-form {viol_S:.2e}, lower {viol_lo:.2e}")
print("   t        <S>      <Var>    <S_N>   1/2ln(2pie(<Var>+1/12))  1/2ln(2pie(<S>/4ln2+1/12))  1/2ln(8<Var>/pi)")
for k, t in enumerate(tgrid):
    if k % 2 == 0:
        s, v, sn = S[:, k].mean(), V[:, k].mean(), SN[:, k].mean()
        print(f"{t:9.1e} {s:8.4f} {v:8.4f} {sn:8.4f}   {0.5*np.log(2*np.pi*np.e*(v+1/12)):8.4f}"
              f"              {0.5*np.log(2*np.pi*np.e*(s/(4*np.log(2))+1/12)):8.4f}              {0.5*np.log(8*v/np.pi) if v>0 else -np.inf:8.4f}")
np.savez(f"rh_L{L}_d{delta}.npz", tgrid=tgrid, S=S, V=V, SN=SN)
