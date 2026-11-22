"""T2: l-bit model H = U D U^dag (U = U(1) brickwork, d_u gates per bond, mixing f;
D = all-k-body diagonal l-bit Hamiltonian with exp(-diam/kappa) couplings).
Quench from a Neel state (random one of the two), half-chain number entropy S_N(t)
on a log grid; infinite-time quantities from random phases (Kronecker-Weyl):
   Shat = E_theta S_N(p_theta)  (= infinite-time average of S_N, D_c rationally indep.)
   Hbar = H(E_theta p_theta)    (MC estimate of H(p_bar), p_bar = dephased distribution)
usage: python3 lbit_run.py tag d_u f kappa L:nreal L:nreal ...
"""
import sys, time
import numpy as np
from multiprocessing import Pool
from lbit import *

TGRID = 10.0**np.arange(-1, 13.01, 0.25)
MPH = 96
TWO_PI = np.longdouble(2)*np.pi*np.longdouble(1)
_SEC = {}

def get_sec(L):
    if L not in _SEC:
        _SEC[L] = Sector(L)
    return _SEC[L]

def one_realization(args):
    L, d_u, f, kappa, seed = args
    rng = np.random.default_rng(seed)
    sec = get_sec(L)
    gates = make_circuit_gates(rng, d_u, f, rel_range=L + 2)
    clist = circuit_list(sec, gates, d_u)
    Dc = lbit_energies_full(rng, L, W=1.0, J0=1.0, kappa=kappa)[sec.states]
    par = rng.integers(2)
    s0 = sec.neel(par)
    psi = np.zeros(sec.D, complex); psi[sec.index(s0)] = 1.0
    n0 = sec.NA[sec.index(s0)]
    phi = apply_Udag(sec, psi, clist)
    Dl = Dc.astype(np.longdouble)
    ph = np.empty((sec.D, len(TGRID)))
    for k, t in enumerate(TGRID):
        ph[:, k] = np.asarray(np.mod(Dl*np.longdouble(t), TWO_PI), dtype=float)
    th = rng.uniform(0, 2*np.pi, size=(sec.D, MPH))
    Psi = np.exp(-1j*np.concatenate([ph, th], axis=1))*phi[:, None]
    Psi = apply_U(sec, Psi, clist)
    P = number_distribution(sec, Psi)
    Pt, Pth = P[:, :len(TGRID)], P[:, len(TGRID):]
    SN_t = shannon(Pt); Var_t = variance(Pt)
    SN_th = shannon(Pth)
    pbar = Pth.mean(1)
    cut = [gates[(lay, 0)] for lay in range(0, 2*d_u, 2)]
    delta_circ = sum(min(1.0, gate_dist_from_identity(g)) for g in cut)
    return dict(SN_t=SN_t, Var_t=Var_t, Shat=SN_th.mean(), Hbar=shannon(pbar),
                Varbar=variance(pbar), delta=delta_circ, n0=n0, pbar=pbar)

if __name__ == "__main__":
    tag, d_u, f, kappa = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    jobs = [(int(a.split(':')[0]), int(a.split(':')[1])) for a in sys.argv[5:]]
    out = {}
    for L, nreal in jobs:
        t0 = time.time()
        args = [(L, d_u, f, kappa, 1000003*L + 7919*d_u + int(1000*f) + s) for s in range(nreal)]
        with Pool(4) as pool:
            res = pool.map(one_realization, args, chunksize=max(1, nreal//32))
        for key in ["SN_t", "Var_t", "Shat", "Hbar", "Varbar", "delta"]:
            out[f"L{L}_{key}"] = np.array([r[key] for r in res])
        print(f"[{tag}] L={L} nreal={nreal} done in {time.time()-t0:.1f}s", flush=True)
    np.savez(f"res_{tag}.npz", tgrid=TGRID, **out)
