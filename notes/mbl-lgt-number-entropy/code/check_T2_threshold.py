import numpy as np
from lbit import *
rng = np.random.default_rng(11)
def dense_U(sec, clist):
    M = np.eye(sec.D, dtype=complex)
    for i, g in clist: apply_gate(sec, M, i, g)
    return M
def pbar_exact(sec, clist, psi):
    U = dense_U(sec, clist); phi = U.conj().T @ psi
    return number_distribution(sec, U) @ (np.abs(phi)**2)
for d_u, Ls in [(1, [4, 6, 8, 10]), (2, [6, 8, 10, 12, 14]), (3, [12, 14, 16])]:
    gates = make_circuit_gates(rng, d_u, 0.9, rel_range=20)
    ref = None
    for L in Ls:
        sec = Sector(L)
        if sec.D > 13000: continue
        clist = circuit_list(sec, gates, d_u)
        par = sec.l % 2
        s0 = sec.neel(par); psi = np.zeros(sec.D, complex); psi[sec.index(s0)] = 1
        n0 = sec.NA[sec.index(s0)]
        pb = pbar_exact(sec, clist, psi)
        vec = np.array([pb[n0 + k] if 0 <= n0 + k < len(pb) else 0 for k in range(-4, 5)])
        print(f"d_u={d_u} L={L:2d} (8d_u-4={8*d_u-4:2d}) p_bar(n0-4..n0+4) =", np.array2string(vec, precision=10, suppress_small=True), "H=", round(shannon(pb), 10))
