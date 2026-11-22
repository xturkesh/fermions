"""
Entropies for the 1D Z2 LGT (open chain + boundary link) after a quench from
Psi0 = |up...up>_z (gauge) (x) |CDW> (matter), computed EXACTLY via the dual representation
   Phi(t) = sum_q 2^{-N/2} |q>_tau (x) psi_q(t),  psi_q(t) = exp(-it[dGamma(h_q) - h sum q]) psi
and the explicit permutation  U: |n>|q> -> |n>|x(n,q)>  back to original variables.

Quantities (half chain A = sites 0..l-1 + links 0..l-1, l = N/2):
  S_full  : vN entropy of region A (matter + links) in ORIGINAL variables (pure global state)
  S_dual  : vN entropy of A' = (c_0..c_{l-1}, tau_0..tau_{l-1}) in dual variables; |S_full-S_dual|<=ln2 (theorem)
  S_G     : E_q S(psi_q restricted to c_A)  (sector-averaged Gaussian entanglement)
  S_mix   : S( E_q rho^q_{c_A} )  (dressed-matter state with charges traced out)
  S_f     : S(rho^f_A), bare matter fermions f_A with all gauge spins traced out (original variables)
  S_N     : Shannon entropy of the distribution of N_A (same in all pictures)
  quenched: E_q S_full for single-sector initial states |x(q)> (x) CDW  (electric-basis gauge state)
"""
import numpy as np, itertools, sys
from math import comb

def entropy_from_rho(r):
    w = np.linalg.eigvalsh(r); w = w[w > 1e-14]
    return float(-(w * np.log(w)).sum())

def entropy_svd(M):
    s = np.linalg.svd(M, compute_uv=False); p = s ** 2; p = p[p > 1e-14]
    return float(-(p * np.log(p)).sum())

def run(N, J, h, times, seed=0):
    l = N // 2
    occ = [j for j in range(0, N, 2)]; Nf = len(occ)
    # number-sector basis
    basis = [c for c in itertools.combinations(range(N), Nf)]
    bidx = {c: i for i, c in enumerate(basis)}
    D = len(basis)
    def fockint(c):  # qubit 0 most significant
        return sum(1 << (N - 1 - j) for j in c)
    full_index = np.array([fockint(c) for c in basis])
    # hopping matrix in number sector (JW signs; nearest neighbour -> sign +1 for open chain hop j<->j+1)
    T = np.zeros((D, D))
    for i, c in enumerate(basis):
        s = set(c)
        for j in range(N - 1):
            if j + 1 in s and j not in s:
                c2 = tuple(sorted((s - {j + 1}) | {j})); T[bidx[c2], i] += -J
                T[i, bidx[c2]] += -J
    occmat = np.array([[1.0 if j in c else 0.0 for j in range(N)] for c in basis])  # D x N
    psi0 = np.zeros(D, complex); psi0[bidx[tuple(occ)]] = 1
    qs = list(itertools.product([1, -1], repeat=N))
    eig = []
    for q in qs:
        Hq = T + np.diag(occmat @ (2 * h * np.array(q, float)))
        E, V = np.linalg.eigh(Hq)
        eig.append((E - h * sum(q), V))
    NA = occmat[:, :l].sum(1).astype(int)
    bits = ((np.arange(2 ** N)[:, None] >> (N - 1 - np.arange(N))[None, :]) & 1)   # [int, j], bit=1 <-> q=-1 / x=-1
    sgn_n = 1 - 2 * bits                                                           # (-1)^{n_j}
    sgn_q = 1 - 2 * bits                                                           # q_j
    xx = np.cumprod(sgn_n[:, None, :] * sgn_q[None, :, :], axis=2)                 # [n, q, b]
    PERM = (((1 - xx) // 2) * (1 << (N - 1 - np.arange(N)))[None, None, :]).sum(2)  # x index
    out = []
    for t in times:
        # dual state Phi[nfull, qindex]
        Phi = np.zeros((2 ** N, 2 ** N), complex)
        SG = 0.0; rho_mix = np.zeros((2 ** l, 2 ** l), complex); pN = np.zeros(Nf + 1)
        for k, q in enumerate(qs):
            E, V = eig[k]
            v = V @ (np.exp(-1j * E * t) * (V.conj().T @ psi0))
            vf = np.zeros(2 ** N, complex); vf[full_index] = v
            qi = int(''.join('0' if a == 1 else '1' for a in q), 2)
            Phi[:, qi] = vf / 2 ** (N / 2)
            M = vf.reshape(2 ** l, 2 ** (N - l))
            SG += entropy_svd(M) / 2 ** N
            rho_mix += M @ M.conj().T / 2 ** N
            np.add.at(pN, NA, np.abs(v) ** 2 / 2 ** N)
        Smix = entropy_from_rho(rho_mix)
        pNn = pN[pN > 1e-15]; SN = float(-(pNn * np.log(pNn)).sum())
        # dual entanglement: index (n bits)(q bits) -> reshape (nA, nB, qA, qB)
        Tn = Phi.reshape(2 ** l, 2 ** (N - l), 2 ** l, 2 ** (N - l)).transpose(0, 2, 1, 3).reshape(4 ** l, -1)
        Sdual = entropy_svd(Tn)
        # original variables: Psi[n, x] with x_b = prod_{k<=b} q_k (-1)^{n_k}  (vectorized permutation)
        Psi = np.zeros((2 ** N, 2 ** N), complex)
        rows = np.repeat(np.arange(2 ** N)[:, None], 2 ** N, axis=1)
        Psi[rows, PERM] = Phi
        Tf = Psi.reshape(2 ** l, 2 ** (N - l), 2 ** l, 2 ** (N - l)).transpose(0, 2, 1, 3).reshape(4 ** l, -1)
        Sfull = entropy_svd(Tf)
        # bare matter f_A: trace spins and f_B
        Mf = Psi.reshape(2 ** l, -1)
        Sf = entropy_from_rho(Mf @ Mf.conj().T)
        out.append((t, Sfull, Sdual, SG, Smix, Sf, SN))
    return out

def quenched(N, J, h, times, nsamp=64, seed=1):
    """single-sector initial states: gauge spins in a sigma^x product state -> definite q.
    By the duality theorem S_full = S(psi_q, c_A) +- ln2; we compute both exactly."""
    rng = np.random.default_rng(seed)
    l = N // 2; occ = list(range(0, N, 2)); Nf = len(occ)
    basis = [c for c in itertools.combinations(range(N), Nf)]; bidx = {c: i for i, c in enumerate(basis)}; D = len(basis)
    full_index = np.array([sum(1 << (N - 1 - j) for j in c) for c in basis])
    T = np.zeros((D, D))
    for i, c in enumerate(basis):
        s = set(c)
        for j in range(N - 1):
            if j + 1 in s and j not in s:
                c2 = tuple(sorted((s - {j + 1}) | {j})); T[bidx[c2], i] += -J; T[i, bidx[c2]] += -J
    occmat = np.array([[1.0 if j in c else 0.0 for j in range(N)] for c in basis])
    psi0 = np.zeros(D, complex); psi0[bidx[tuple(occ)]] = 1
    nbits = ((np.arange(2 ** N)[:, None] >> (N - 1 - np.arange(N))[None, :]) & 1)
    res = np.zeros((len(times), 3))
    for s_ in range(nsamp):
        q = rng.choice([1, -1], size=N)
        E, V = np.linalg.eigh(T + np.diag(occmat @ (2 * h * q)))
        for it, t in enumerate(times):
            v = V @ (np.exp(-1j * E * t) * (V.conj().T @ psi0))
            vf = np.zeros(2 ** N, complex); vf[full_index] = v
            Sc = entropy_svd(vf.reshape(2 ** l, -1))
            Psi = np.zeros((2 ** N, 2 ** N), complex)
            for ni in np.nonzero(np.abs(vf) > 0)[0]:
                x = np.cumprod(q * (1 - 2 * nbits[ni]))
                xi = int(''.join('0' if a == 1 else '1' for a in x), 2)
                Psi[ni, xi] += vf[ni]
            Sfull = entropy_svd(Psi.reshape(2 ** l, 2 ** (N - l), 2 ** l, 2 ** (N - l)).transpose(0, 2, 1, 3).reshape(4 ** l, -1))
            res[it] += np.array([Sc, Sfull, abs(Sfull - Sc)]) / nsamp
    return res

if __name__ == '__main__':
    J = 1.0
    times = [0, 1, 3, 10, 30, 100, 300, 1e3, 1e4, 1e5]
    for N, h in [(8, 1.0), (10, 1.0), (10, 2.0)]:
        print(f'=== N={N}, h/J={h} (binary potential +-{2*h}), half chain l={N//2}; ln2={np.log(2):.3f}, l*ln2={N//2*np.log(2):.3f}')
        print('      t   S_full  S_dual |diff|   S_G    S_mix   S_f(bare) S_N')
        for (t, Sfull, Sdual, SG, Smix, Sf, SN) in run(N, J, h, times):
            print(f'{t:8.0f} {Sfull:7.3f} {Sdual:7.3f} {abs(Sfull-Sdual):5.3f} {SG:7.3f} {Smix:7.3f} {Sf:8.3f} {SN:7.3f}')
        sys.stdout.flush()
    print('=== quenched (single-sector, sigma^x-product gauge state) N=8 h=1: E_q S(psi_q,c_A), E_q S_full, E_q|diff|')
    r = quenched(8, J, 1.0, times, nsamp=32)
    for t, row in zip(times, r):
        print(f'{t:8.0f} {row[0]:7.3f} {row[1]:7.3f} {row[2]:7.3f}')
