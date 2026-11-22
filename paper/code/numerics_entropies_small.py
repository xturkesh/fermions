"""
Numerics for Section sec:numerics, Figure fig:entropies(a): exact entropies of the Z2 LGT for N = 10.

Conventions (paper/BRIEF.md, Sec. 2): sites x = 1..N, links b_j = (j, j+1), j = 1..N (b_N = boundary link),
    H = -J sum_{j<N} sigma^z_{b_j} (f_j^+ f_{j+1} + h.c.) - h sum_j sigma^x_{b_{j-1}} sigma^x_{b_j},  sigma^x_{b_0} = 1,
Psi_0 = |up...up>_z (x) CDW (n_x = 1 iff x even), region A = sites 1..l and links b_1..b_l, l = N/2.

Method: exact reconstruction through the duality (thm:duality, cor:sector).  For every one of the 2^N
sectors q, psi_q(t) = Gamma(U_q(t)) psi is a Slater determinant (amplitudes = determinants of U_q(t)), and
    Phi(t) = 2^{-N/2} sum_q exp(i h t sum_j q_j) psi_q(t) (x) |q>,       Psi(t) = U Phi(t),
    U(|n> (x) |q>) = |n> (x) |x(n,q)>_sigma,   x_{b_j}(n,q) = prod_{k<=j} q_k (-1)^{n_k}   (eq:U, eq:xnq).
In bits (1 <-> value -1, site/link 1 = most significant bit): q = n XOR x XOR (x >> 1).
Quantities (prop:entropies):
    S_full = S_A(Psi(t)) (matter + links, original variables), S_dual = S_{A'}(Phi(t)),
    S_G = E_q S(rho^q_A), S_mix = S(E_q rho^q_A), S_f (bare matter, all links traced out),
    S_N = H(p_t), p_t = law of N_A = sum_{x<=l} n_x;  also E_q H(p^q_t).
The pipeline is checked against direct exact diagonalization of the full LGT (N = 4 and 6): the state
Psi(t) itself is compared entrywise (links rotated to the sigma^z basis).  For N = 6 we also check
S_full = S(rho^q_A) for all 64 definite-charge initial states (thm:arealaw(c)).

Usage:  python3 numerics_entropies_small.py  > numerics_entropies_small.out
        -> data/numerics_entropies_small.npz.   Runtime ~2 min.  Nothing here is used in a proof.
"""
import itertools
import os
import time

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.sparse.linalg import expm_multiply

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
LN2 = np.log(2.0)


def P(*a):
    print(*a, flush=True)


def hbin(nu):
    nu = np.clip(nu, 1e-300, 1.0 - 1e-16)
    return -nu * np.log(nu) - (1.0 - nu) * np.log1p(-nu)


def shannon(p):
    p = np.asarray(p, float)
    p = p[p > 1e-300]
    return float(-(p * np.log(p)).sum())


def vn_from_svals(s):
    p = s ** 2
    p = p[p > 1e-16]
    return float(-(p * np.log(p)).sum())


def vn_from_rho(r):
    w = np.linalg.eigvalsh(r)
    w = w[w > 1e-16]
    return float(-(w * np.log(w)).sum())


def poibin(nu):
    p = np.zeros(len(nu) + 1)
    p[0] = 1.0
    for k, v in enumerate(nu):
        p[1:k + 2] = p[1:k + 2] * (1 - v) + p[0:k + 1] * v
        p[0] *= (1 - v)
    return p


def popcount(a):
    a = np.asarray(a, np.int64)
    c = np.zeros_like(a)
    while np.any(a):
        c += a & 1
        a = a >> 1
    return c


class Pipeline:
    """Exact dual reconstruction for a chain of N sites, all 2^N sectors."""

    def __init__(self, N, J, h, ell, occ):
        self.N, self.J, self.h, self.ell = N, J, h, ell
        self.occ = list(occ)                                    # occupied sites, 0-based (x-1)
        self.M = len(self.occ)
        confs = list(itertools.combinations(range(N), self.M))
        self.rows = np.array(confs)
        self.cidx = np.array([sum(1 << (N - 1 - j) for j in c) for c in confs])
        self.qs = np.array(list(itertools.product([1.0, -1.0], repeat=N)))   # index k: bit 1 <-> q=-1, MSB = site 1
        self.eig = [eigh_tridiagonal(2.0 * h * q, -J * np.ones(N - 1)) for q in self.qs]
        n = np.arange(2 ** N)
        self.nidx = n
        self.Qi = n[:, None] ^ n[None, :] ^ (n[None, :] >> 1)    # q-index for (n, x)
        self.NAof = popcount(n >> (N - ell))                     # N_A of a matter configuration
        self.NBof = popcount(n & ((1 << (N - ell)) - 1))

    def U(self, k, t):
        E, V = self.eig[k]
        return (V * np.exp(-1j * E * t)) @ V.T

    def sector_vector(self, k, t):
        U = self.U(k, t)
        v = np.zeros(2 ** self.N, complex)
        v[self.cidx] = np.linalg.det(U[self.rows][:, :, self.occ])
        return v, U

    def blocks_entropy(self, T4):
        """entanglement of a pure state T4[nA, nB, xA, xB] across (nA,xA)|(nB,xB), block-diagonal in N_A."""
        N, ell = self.N, self.ell
        nA = np.arange(2 ** ell)
        nB = np.arange(2 ** (N - ell))
        cA = popcount(nA)
        cB = popcount(nB)
        sv = []
        for kA in range(self.M + 1):
            ia = np.nonzero(cA == kA)[0]
            ib = np.nonzero(cB == self.M - kA)[0]
            if len(ia) == 0 or len(ib) == 0:
                continue
            Mb = T4[np.ix_(ia, ib)]                               # [ia, ib, xA, xB]
            Mb = Mb.transpose(0, 2, 1, 3).reshape(len(ia) * 2 ** ell, len(ib) * 2 ** (N - ell))
            sv.append(np.linalg.svd(Mb, compute_uv=False))
        return vn_from_svals(np.concatenate(sv))

    def Phi(self, t, sectors=None):
        N = self.N
        Phi = np.zeros((2 ** N, 2 ** N), complex)
        ks = range(2 ** N) if sectors is None else sectors
        for k in ks:
            v, _ = self.sector_vector(k, t)
            Phi[:, k] = np.exp(1j * self.h * t * self.qs[k].sum()) * v / 2 ** (N / 2)
        return Phi

    def Psi_from_Phi(self, Phi):
        return Phi[self.nidx[:, None], self.Qi]

    def all_entropies(self, t):
        N, ell = self.N, self.ell
        Phi = np.zeros((2 ** N, 2 ** N), complex)
        SG = 0.0
        EHq = 0.0
        rho_mix = np.zeros((2 ** ell, 2 ** ell), complex)
        dev_slater = 0.0
        for k in range(2 ** N):
            v, U = self.sector_vector(k, t)
            Phi[:, k] = np.exp(1j * self.h * t * self.qs[k].sum()) * v / 2 ** (N / 2)
            K = U[:, self.occ] @ U[:, self.occ].conj().T
            nu = np.clip(np.linalg.eigvalsh(K[:ell, :ell]), 0.0, 1.0)
            sq = float(hbin(nu).sum())
            SG += sq / 2 ** N
            Mq = v.reshape(2 ** ell, 2 ** (N - ell))
            rq = Mq @ Mq.conj().T
            rho_mix += rq / 2 ** N
            pq = np.bincount(self.NAof, weights=np.abs(v) ** 2, minlength=ell + 1)[:ell + 1]
            EHq += shannon(pq) / 2 ** N
            if k % 97 == 0:          # lem:slater(a): many-body reduced state vs correlation matrix
                dev_slater = max(dev_slater, abs(vn_from_rho(rq) - sq),
                                 np.max(np.abs(poibin(nu) - pq)))
        Smix = vn_from_rho(rho_mix)
        T4 = Phi.reshape(2 ** ell, 2 ** (N - ell), 2 ** ell, 2 ** (N - ell))
        Sdual = self.blocks_entropy(T4)
        Psi = self.Psi_from_Phi(Phi)
        Sfull = self.blocks_entropy(Psi.reshape(2 ** ell, 2 ** (N - ell), 2 ** ell, 2 ** (N - ell)))
        Mf = Psi.reshape(2 ** ell, -1)
        Sf = vn_from_rho(Mf @ Mf.conj().T)
        pN = np.bincount(self.NAof, weights=(np.abs(Psi) ** 2).sum(1), minlength=ell + 1)[:ell + 1]
        SN = shannon(pN)
        norm = float((np.abs(Psi) ** 2).sum())
        return dict(Sfull=Sfull, Sdual=Sdual, SG=SG, Smix=Smix, Sf=Sf, SN=SN, EHq=EHq,
                    norm=norm, dev_slater=dev_slater)

    def SG_SN(self, t):
        """S_G and S_N only, from correlation matrices and Poisson-binomial laws (lem:slater)."""
        ell = self.ell
        SG = 0.0
        p = np.zeros(ell + 1)
        for k in range(2 ** self.N):
            U = self.U(k, t)
            Us = U[:, self.occ]
            nu = np.clip(np.linalg.eigvalsh((Us @ Us.conj().T)[:ell, :ell]), 0.0, 1.0)
            SG += hbin(nu).sum() / 2 ** self.N
            p += poibin(nu) / 2 ** self.N
        return float(SG), shannon(p)


# --------------------------------------------------------------------------- direct ED of the full LGT
def ed_lgt(N, J, h):
    """sparse H on (fermions f_1..f_N in JW order) (x) (links b_1..b_N in the sigma^z basis, 0 = up)."""
    I2 = sp.identity(2, format="csr")
    Z = sp.csr_matrix(np.diag([1.0, -1.0]))
    A = sp.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))         # annihilation on (|0>, |1>)
    X = sp.csr_matrix(np.array([[0.0, 1.0], [1.0, 0.0]]))

    def kron_list(ops):
        out = ops[0]
        for o in ops[1:]:
            out = sp.kron(out, o, format="csr")
        return out

    def f(j):          # j = 1..N
        ops = [Z] * (j - 1) + [A] + [I2] * (N - j) + [I2] * N
        return kron_list(ops)

    def link(op, j):   # j = 1..N
        ops = [I2] * N + [I2] * (j - 1) + [op] + [I2] * (N - j)
        return kron_list(ops)

    D = 4 ** N
    H = sp.csr_matrix((D, D))
    for j in range(1, N):
        hop = f(j).T @ f(j + 1)
        H = H - J * link(Z, j) @ (hop + hop.T)
    for j in range(1, N + 1):
        term = link(X, j) if j == 1 else link(X, j - 1) @ link(X, j)
        H = H - h * term
    return H.tocsr()


def ed_check(N, J, h, ell, times):
    occ = [x - 1 for x in range(1, N + 1) if x % 2 == 0]
    pl = Pipeline(N, J, h, ell, occ)
    H = ed_lgt(N, J, h)
    ncdw = sum(1 << (N - 1 - j) for j in occ)
    psi0 = np.zeros(4 ** N, complex)
    psi0[ncdw * 2 ** N + 0] = 1.0                                      # links all up (s = 0)
    Had = np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2.0)          # T[x, s] = <s|x>
    T = Had
    for _ in range(N - 1):
        T = np.kron(T, Had)
    errs = []
    for t in times:
        ped = expm_multiply(-1j * t * H, psi0).reshape(2 ** N, 2 ** N)
        Psi_x = pl.Psi_from_Phi(pl.Phi(t))
        Psi_z = Psi_x @ T
        e_state = float(np.max(np.abs(Psi_z - ped)))
        Ted = ped.reshape(2 ** ell, 2 ** (N - ell), 2 ** ell, 2 ** (N - ell))
        S_ed = pl.blocks_entropy(Ted)
        r = pl.all_entropies(t)
        Mf = ped.reshape(2 ** ell, -1)
        Sf_ed = vn_from_rho(Mf @ Mf.conj().T)
        pN = np.bincount(pl.NAof, weights=(np.abs(ped) ** 2).sum(1), minlength=ell + 1)[:ell + 1]
        errs.append((t, e_state, abs(S_ed - r["Sfull"]), abs(Sf_ed - r["Sf"]), abs(shannon(pN) - r["SN"])))
    return errs


def quenched_check(N, J, h, ell, times):
    """S_full of U(psi_q(t) (x) |q>) versus S(rho^q_A), all q (definite-charge initial states)."""
    occ = [x - 1 for x in range(1, N + 1) if x % 2 == 0]
    pl = Pipeline(N, J, h, ell, occ)
    worst = 0.0
    for t in times:
        for k in range(2 ** N):
            v, U = pl.sector_vector(k, t)
            Phi = np.zeros((2 ** N, 2 ** N), complex)
            Phi[:, k] = np.exp(1j * h * t * pl.qs[k].sum()) * v
            Psi = pl.Psi_from_Phi(Phi)
            Sf = pl.blocks_entropy(Psi.reshape(2 ** ell, 2 ** (N - ell), 2 ** ell, 2 ** (N - ell)))
            K = U[:, occ] @ U[:, occ].conj().T
            sq = float(hbin(np.clip(np.linalg.eigvalsh(K[:ell, :ell]), 0, 1)).sum())
            worst = max(worst, abs(Sf - sq))
    return worst


def main():
    os.makedirs(DATA, exist_ok=True)
    J = 1.0
    P("=== cross-check of the dual reconstruction against direct ED of the full LGT ===")
    for (N, h, ell) in ((4, 1.3, 2), (6, 1.0, 3), (6, -0.7, 2)):
        for (t, es, e1, e2, e3) in ed_check(N, J, h, ell, (0.3, 1.7, 6.0)):
            P(f"  N={N} h={h:+.1f} l={ell} t={t:4.1f}: max|Psi_pipeline - Psi_ED| = {es:.1e}; "
              f"|dS_full| = {e1:.1e}, |dS_f| = {e2:.1e}, |dS_N| = {e3:.1e}")
    wq = quenched_check(6, J, 1.0, 3, (0.7, 5.0))
    P(f"  N=6, h=J, l=3, all 64 definite-charge states, t=0.7,5: max|S_full - S(rho^q_A)| = {wq:.1e}")

    N, h, ell = 10, 1.0, 5
    occ = [x - 1 for x in range(1, N + 1) if x % 2 == 0]
    P(f"=== exact entropies: N={N}, J={J}, h={h} (lambda={2*h/J}), l={ell}, CDW on even sites, "
      f"all {2**N} sectors ===")
    t0 = time.time()
    pl = Pipeline(N, J, h, ell, occ)
    times = np.concatenate([[0.0], np.logspace(-2, 6, 97)])
    keys = ("Sfull", "Sdual", "SG", "Smix", "Sf", "SN", "EHq")
    res = {k: np.zeros(len(times)) for k in keys}
    viol = dict(cut=-np.inf, Gmix=-np.inf, mixmax=-np.inf, SN=-np.inf, SGq=-np.inf, EHqSN=-np.inf,
                norm=0.0, slater=0.0)
    for it, t in enumerate(times):
        r = pl.all_entropies(t)
        for k in keys:
            res[k][it] = r[k]
        if t == 0.0:
            P(f"  t=0: all entropies vanish: max = {max(abs(r[k]) for k in keys):.1e}")
            continue
        viol["cut"] = max(viol["cut"], abs(r["Sfull"] - r["Sdual"]) - LN2)
        viol["Gmix"] = max(viol["Gmix"], r["SG"] - r["Smix"])
        viol["mixmax"] = max(viol["mixmax"], r["Smix"] - ell * LN2)
        viol["SN"] = max(viol["SN"], r["SN"] - min(r["Sfull"], r["Smix"], r["Sf"]))
        viol["SGq"] = max(viol["SGq"], r["EHq"] - r["SG"])
        viol["EHqSN"] = max(viol["EHqSN"], r["EHq"] - r["SN"])
        viol["norm"] = max(viol["norm"], abs(r["norm"] - 1.0))
        viol["slater"] = max(viol["slater"], r["dev_slater"])
    P(f"  ({time.time() - t0:.0f}s)   t      S_full  S_A'(Phi)  S_G    S_mix   S_f    S_N   E_qH(p^q)")
    for tt in (0.5, 1.0, 2.0, 10.0, 100.0, 1e3, 1e5, 1e6):
        i = int(np.argmin(np.abs(times - tt)))
        P(f"  {times[i]:10.3g}  " + "  ".join(f"{res[k][i]:6.3f}" for k in keys))
    P(f"  prop:entropies checks on all {len(times) - 1} times t > 0 (all must be <= 0): "
      f"|S_full-S_A'|-ln2: {viol['cut']:.2e}; S_G-S_mix: {viol['Gmix']:.2e}; S_mix-l ln2: {viol['mixmax']:.2e}; "
      f"S_N-min(S_full,S_mix,S_f): {viol['SN']:.2e}; E_qH(p^q)-S_G: {viol['SGq']:.2e}; E_qH(p^q)-S_N: {viol['EHqSN']:.2e}")
    P(f"  norm error {viol['norm']:.1e}; lem:slater (many-body vs correlation matrix) max deviation {viol['slater']:.1e}")
    late = times >= 1e3
    for k in keys:
        P(f"  {k:6s}: mean over t in [1e3,1e6] = {res[k][late].mean():.3f}, min = {res[k][late].min():.3f}, "
          f"max = {res[k][late].max():.3f}")
    P(f"  maxima: 2 l ln2 = {2*ell*LN2:.3f} (S_full), l ln2 = {ell*LN2:.3f} (S_mix, S_f)")
    # sign changes of S_N - S_G on the grid, refined by root finding (S_G, S_N recomputed exactly)
    d = res["SN"] - res["SG"]
    roots = []
    for i in range(1, len(times) - 1):
        if d[i] * d[i + 1] < 0:
            f = lambda tt: (lambda a: a[1] - a[0])(pl.SG_SN(tt))
            roots.append(brentq(f, times[i], times[i + 1], xtol=1e-6))
    pos = times[(d > 0) & (times > 0)]
    P(f"  log grid: S_N - S_G > 0 for t in [{pos.min():.3g}, {pos.max():.3g}]; sign changes at t = "
      + ", ".join(f"{r:.4f}" for r in roots) + f";  max_t (S_N - S_G) = {d.max():.4f} at t = {times[np.argmax(d)]:.3g}")
    # fine linear grid on [0.05, 12] (S_G, S_N only)
    tf = np.arange(0.05, 12.0 + 1e-9, 0.005)
    df = np.array([(lambda a: a[1] - a[0])(pl.SG_SN(tt)) for tt in tf])
    ch = np.nonzero(np.sign(df[1:]) != np.sign(df[:-1]))[0]
    P(f"  fine grid (step 0.005 on [0.05, 12]): S_N - S_G > 0 for all t in [0.05, {tf[ch[0]]:.3f}]; sign changes near t = "
      + ", ".join(f"{tf[i]:.3f}" for i in ch) + f"; S_N - S_G < 0 on the whole grid for t > {tf[ch[-1] + 1]:.3f}: "
      + f"{bool(np.all(df[ch[-1] + 1:] < 0))};  min over t in [6, 1e6] of S_G - S_N on the log grid: "
      + f"{np.min(-d[times >= 6]):.3f}")
    np.savez_compressed(os.path.join(DATA, "numerics_entropies_small.npz"), times=times, N=N, J=J, h=h, ell=ell,
                        crossings=np.array(roots), tfine=tf, dfine=df, **res)


if __name__ == "__main__":
    main()
