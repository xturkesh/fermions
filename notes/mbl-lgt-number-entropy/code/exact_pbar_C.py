import numpy as np
from lbit import *
rng = np.random.default_rng(2024)
L, d_u, f = 8, 1, 0.8
sec = Sector(L)
Hs = []; Vs = []
for r in range(4000):
    gates = make_circuit_gates(rng, d_u, f, rel_range=L+2)
    clist = circuit_list(sec, gates, d_u)
    U = np.eye(sec.D, dtype=complex)
    for i, g in clist: apply_gate(sec, U, i, g)
    par = rng.integers(2); s0 = sec.neel(par)
    psi = np.zeros(sec.D, complex); psi[sec.index(s0)] = 1
    phi = U.conj().T @ psi
    pb = number_distribution(sec, U) @ np.abs(phi)**2
    Hs.append(shannon(pb)); Vs.append(variance(pb))
Hs = np.array(Hs); Vs = np.array(Vs)
print(f"exact E H(pbar) (d_u=1,f=0.8, L=8; L-independent for L>=4) = {Hs.mean():.4f} +- {Hs.std()/np.sqrt(len(Hs)):.4f};  E Var_pbar = {Vs.mean():.4f}")
