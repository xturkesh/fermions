"""
Numerical sanity checks for Section sec:main ("Main results: disorder-free localization").

Every inequality / identity stated in sections/main.tex is tested here, either exactly (ED of the full
lattice gauge theory for small N, using the independent builders of model_duality_check.py) or on
random sectors of the one-body problem h_q = -J Delta_N + 2h diag(q).  Conventions as in BRIEF.md:
sites x = 1..N (0-based index x-1 in the code), CDW = even sites occupied, Gamma = diag((-1)^x).

Nothing printed here is used in a proof.  Quantities called "F_emp" are statistical estimates of the
function F of Hypothesis (EC) (sample means of the correlator, made non-increasing); bounds evaluated with
F_emp are illustrations of the size of the constants, not rigorous statements.

Usage: python3 main_checks.py > main_checks.out      (runtime: a few minutes on 4 cores)
"""
import itertools
import math
import os
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import model_duality_check as mdc   # noqa: E402  (ED builders: physical(), build_U(), ...)

LN2 = np.log(2.0)
RNG = np.random.default_rng(20261122)


def P(*a):
    print(*a, flush=True)


# --------------------------------------------------------------------------- one-body helpers
def eig_q(q, J, h):
    N = len(q)
    E, V = eigh_tridiagonal(2.0 * h * np.asarray(q, float), -J * np.ones(N - 1))
    return E, V


def Uof(E, V, t):
    return (V * np.exp(-1j * E * t)) @ V.T


def corr(V):
    A = np.abs(V)
    return A @ A.T                      # Q(x,y) = sum_a |phi_a(x) phi_a(y)|


def Rmat(V):
    B = V ** 2
    return B @ B.T                      # R(x,y) = sum_a phi_a(x)^2 phi_a(y)^2


def hfun(nu):
    nu = np.clip(nu, 1e-300, 1 - 1e-16)
    return -nu * np.log(nu) - (1 - nu) * np.log(1 - nu)


def shannon(p):
    p = np.asarray(p, float)
    p = p[p > 1e-300]
    return float(-(p * np.log(p)).sum())


def vn_entropy_from_matrix(M):
    s = np.linalg.svd(M, compute_uv=False)
    p = s ** 2
    p = p[p > 1e-15]
    return float(-(p * np.log(p)).sum())


def vn_entropy_rho(r):
    w = np.linalg.eigvalsh(r)
    w = w[w > 1e-15]
    return float(-(w * np.log(w)).sum())


def poibin(nu):
    """pmf of a sum of independent Bernoulli(nu_i)."""
    p = np.zeros(len(nu) + 1)
    p[0] = 1.0
    for k, v in enumerate(nu):
        p[1:k + 2] = p[1:k + 2] * (1 - v) + p[0:k + 1] * v
        p[0] *= (1 - v)
    return p


def F_envelope(Fd):
    """non-increasing majorant min(1, max_{d'>=d} Fd(d'))."""
    out = np.minimum(1.0, np.maximum.accumulate(Fd[::-1])[::-1])
    return out


def moments(F):
    d = np.arange(len(F))
    M = {p: float(((d + 1.0) ** p * F).sum()) for p in (0, 1, 2)}
    Mp = {p: float(((d + 1.0) ** p * np.sqrt(F)).sum()) for p in (0, 1, 2)}
    S1 = float((d[1:] * F[1:]).sum())
    S1p = float((d[1:] * np.sqrt(F[1:])).sum())
    Theta = float(2 * (d[1:] * (2 * d[1:] + 1) * F[1:]).sum())
    return M, Mp, S1, S1p, Theta


def exp_majorant(Fd, dmin=1):
    """smallest C with Fd(d) <= C e^{-eta d} for the eta of a least-squares fit (illustrative)."""
    d = np.arange(len(Fd))
    m = (d >= dmin) & (Fd > 0)
    eta = -np.polyfit(d[m], np.log(Fd[m]), 1)[0]
    C = float(np.max(Fd[m] * np.exp(eta * d[m])))
    return max(C, 1.0), float(eta)


TIMES = np.concatenate([[0.0], np.logspace(-1, 5, 61)])


# =========================================================================== A. transport, moments, F_emp
def part_A(N=60, J=1.0, h=1.0, nsamp=300):
    P(f"=== A. Theorem thm:transport (N={N}, J={J}, h={h}, lambda={2*h/J}, {nsamp} sectors, "
      f"{len(TIMES)} times in [0,1e5]) ===")
    S = np.array([1 if (x + 1) % 2 == 0 else 0 for x in range(N)])      # CDW, 1-based even sites
    Qsum = np.zeros((N, N))
    worst = dict(tr=0.0, D=0.0, mom=0.0, sym=0.0)
    sup_tr_mid = []
    sup_mom = {p: [] for p in (0, 1, 2)}
    y0 = N // 2
    dist = np.abs(np.arange(N)[:, None] - np.arange(N)[None, :])
    for _ in range(nsamp):
        q = RNG.choice([-1.0, 1.0], size=N)
        E, V = eig_q(q, J, h)
        Q = corr(V)
        Qsum += Q
        cutQ = np.array([Q[:m, m:].sum() for m in range(1, N)])       # sum_{x<=m<z} Q(x,z)
        smid = 0.0
        smom = {p: 0.0 for p in (0, 1, 2)}
        for t in TIMES:
            U = Uof(E, V, t)
            W = np.abs(U) ** 2
            worst['sym'] = max(worst['sym'], np.max(np.abs(U - U.T)))
            nx = W[:, S == 1].sum(1)
            for m in range(1, N):
                tr = nx[:m].sum() - S[:m].sum()
                D = W[:m, m:].sum()
                worst['tr'] = max(worst['tr'], abs(tr) - D)
                worst['D'] = max(worst['D'], D - cutQ[m - 1])
                if m == N // 2:
                    smid = max(smid, abs(tr))
            for p in (0, 1, 2):
                mom = ((dist[:, y0] + 1.0) ** p * W[:, y0]).sum()
                bnd = ((dist[:, y0] + 1.0) ** p * Q[:, y0]).sum()
                worst['mom'] = max(worst['mom'], mom - bnd)
                smom[p] = max(smom[p], mom)
        sup_tr_mid.append(smid)
        for p in (0, 1, 2):
            sup_mom[p].append(smom[p])
    Qm = Qsum / nsamp
    Fd = np.array([max(Qm[x, x + d] for x in range(N - d)) for d in range(N)])
    F = F_envelope(Fd)
    M, Mp, S1, S1p, Theta = moments(F)
    P(f"  max_t,m,q (|<N_<=m>_q - |S cap [1,m]|| - sum_{{x<=m<z}}|U_xz|^2) = {worst['tr']:.2e}  (<= 0 required)")
    P(f"  max (sum_{{x<=m<z}}|U_xz|^2 - sum_{{x<=m<z}} Q(x,z))               = {worst['D']:.2e}  (<= 0)")
    P(f"  max (sum_x (|x-y|+1)^p |U_xy|^2 - sum_x (|x-y|+1)^p Q(x,y)), p<=2   = {worst['mom']:.2e}  (<= 0)")
    P(f"  symmetry of U_q(t): max|U - U^T| = {worst['sym']:.1e}")
    P(f"  F_emp (envelope of max_x E Q(x,x+d)): F(1,5,10,20,40) = "
      + ", ".join(f"{F[d]:.3g}" for d in (1, 5, 10, 20, 40)))
    P(f"  moments of F_emp (truncated at d<N): M_0={M[0]:.3g} M_1={M[1]:.3g} M_2={M[2]:.3g} "
      f"M'_1={Mp[1]:.3g}  Sigma_1={S1:.3g}  Sigma'_1={S1p:.3g}  Theta={Theta:.3g}")
    P(f"  E sup_t |transport across the middle cut| = {np.mean(sup_tr_mid):.3f}   "
      f"(theorem: <= Sigma_1(F) <= M_1(F); with F_emp: {S1:.3g})")
    for p in (0, 1, 2):
        P(f"  E sup_t sum_x (|x-y|+1)^{p}|U_xy|^2 (y=N/2) = {np.mean(sup_mom[p]):.3f}   "
          f"(theorem: <= 2 M_{p}(F); with F_emp: {2*M[p]:.3g})")
    return F


# =========================================================================== B. ED of the LGT: anticommutator, norms
def lgt_setup(N, J, h, with_P=True):
    Pd = mdc.physical(N, J, h)
    H = Pd['H'].toarray()
    w, V = np.linalg.eigh(H)
    Id = Pd['Id']
    ct = []
    for x in range(N):
        S = Id
        for k in range(x, N):
            S = S @ Pd['sz'][k]
        ct.append((Pd['f'][x] @ S).toarray())
    G = [g.toarray() for g in Pd['G']]
    n = [m.toarray() for m in Pd['n']]
    D = H.shape[0]
    Pq = {}
    for q in (itertools.product([1, -1], repeat=N) if with_P else []):
        Pm = np.eye(D)
        for j in range(N):
            Pm = Pm @ (np.eye(D) + q[j] * G[j]) / 2
        Pq[q] = Pm
    return dict(H=H, w=w, V=V, ct=ct, G=G, n=n, Pq=Pq, D=D, nq=Pd['nq'])


def heis(L, X, t):
    Ut = (L['V'] * np.exp(-1j * L['w'] * t)) @ L['V'].conj().T
    return Ut.conj().T @ X @ Ut


def opnorm(X):
    return float(np.linalg.norm(X, 2))


def part_B(N=4, J=1.0, h=0.7):
    P(f"=== B. ED of the LGT (N={N}, J={J}, h={h}): propagator identity and operator norms ===")
    L = lgt_setup(N, J, h)
    err_anti = 0.0
    err_cn = 0.0
    err_nn = 0.0
    lr_viol = -np.inf
    st_viol = -np.inf
    qs = list(L['Pq'])
    # state Psi0 = |up..up>_z x CDW
    nocc = [1 if (x + 1) % 2 == 0 else 0 for x in range(N)]
    psi0 = np.zeros(L['D'], complex)
    psi0[mdc.index(nocc + [0] * N)] = 1.0
    A0 = L['ct'][0].conj().T @ L['ct'][1]
    A0 = A0 + A0.conj().T                      # c~_1^+ c~_2 + h.c.
    B0 = L['n'][N - 1]                         # n_N
    for t in (0.3, 1.0, 2.7, 7.0):
        Uq = {q: Uof(*eig_q(np.array(q, float), J, h), t) for q in qs}
        for x in range(N):
            cxt = heis(L, L['ct'][x], t)
            for y in range(N):
                anti = cxt @ L['ct'][y].conj().T + L['ct'][y].conj().T @ cxt
                pred = sum(Uq[q][x, y] * L['Pq'][q] for q in qs)
                err_anti = max(err_anti, np.max(np.abs(anti - pred)))
                cn = opnorm(cxt @ L['n'][y] - L['n'][y] @ cxt)
                err_cn = max(err_cn, abs(cn - max(abs(Uq[q][x, y]) for q in qs)))
                nxt = heis(L, L['n'][x], t)
                nn = opnorm(nxt @ L['n'][y] - L['n'][y] @ nxt)
                pr = max(abs(Uq[q][x, y]) * np.sqrt(max(0.0, 1 - abs(Uq[q][x, y]) ** 2)) for q in qs)
                err_nn = max(err_nn, abs(nn - pr))
        # sector-averaged LR bound for A = c~_1^+c~_2 + h.c., B = n_N (monomial form)
        At = heis(L, A0, t)
        Cm = At @ B0 - B0 @ At
        avg = np.mean([opnorm(L['Pq'][q] @ Cm) for q in qs])
        bnd = np.mean([2 * 1.0 * 2 * (abs(Uq[q][0, N - 1]) + abs(Uq[q][1, N - 1])) for q in qs])
        lr_viol = max(lr_viol, avg - bnd)
        st_viol = max(st_viol, abs(np.vdot(psi0, Cm @ psi0)) - avg)
    P(f"  max |{{c~_x(t),c~_y^+}} - sum_q U_q(t)_xy P_q|            = {err_anti:.1e}")
    P(f"  max | ||[c~_x(t),n_y]|| - max_q |U_q(t)_xy| |             = {err_cn:.1e}")
    P(f"  max | ||[n_x(t),n_y]|| - max_q |u|(1-|u|^2)^(1/2) |       = {err_nn:.1e}")
    P(f"  sector-averaged LR (A=c~_1^+c~_2+h.c., B=n_N): max(E_q||P_q[A(t),B]|| - bound) = {lr_viol:.2e} (<=0)")
    P(f"  state bound: max(|<Psi0,[A(t),B]Psi0>| - E_q||P_q[A(t),B]||) = {st_viol:.2e} (<=0)")
    # definite-charge states: S_full = S(rho^q_A)
    err_q = 0.0
    for NN in (4, 5):
        LL = lgt_setup(NN, J, h, with_P=False) if NN != N else L
        l = 2
        for trial in range(4):
            xi = RNG.choice([1, -1], size=NN)
            nn = [int(b) for b in RNG.integers(0, 2, size=NN)]
            gauge = np.array([1.0])
            for j in range(NN):
                gauge = np.kron(gauge, np.array([1.0, xi[j]]) / np.sqrt(2))
            matter = np.zeros(2 ** NN)
            matter[mdc.index(nn)] = 1
            psi = np.kron(matter, gauge).astype(complex)
            c0 = LL['V'].conj().T @ psi
            q = np.array([(xi[j - 1] if j > 0 else 1) * (-1) ** nn[j] * xi[j] for j in range(NN)], float)
            E, V = eig_q(q, J, h)
            occ = [j for j in range(NN) if nn[j]]
            for t in (0.5, 2.0, 9.0):
                pt = LL['V'] @ (np.exp(-1j * LL['w'] * t) * c0)
                T = pt.reshape([2] * (2 * NN))
                order = list(range(l)) + list(range(NN, NN + l)) + list(range(l, NN)) + list(range(NN + l, 2 * NN))
                M = T.transpose(order).reshape(4 ** l, -1)
                Sfull = vn_entropy_from_matrix(M)
                U = Uof(E, V, t)
                K = U[:, occ] @ U[:, occ].conj().T
                nu = np.linalg.eigvalsh(K[:l, :l]) if occ else np.zeros(l)
                SG = float(hfun(nu).sum())
                err_q = max(err_q, abs(Sfull - SG))
    P(f"  definite-charge states |xi>_sigma x |n> (N=4,5, l=2): max |S_full(t) - sum_i h(nu_i)| = {err_q:.1e}")


# =========================================================================== C. LR in free fermions, clean sector
def dGamma_norm_herm(W):
    w = np.linalg.eigvalsh(W)
    return max(w[w > 0].sum(), -w[w < 0].sum())


def part_C(N=60, J=1.0, h=1.0, nsamp=200):
    P(f"=== C. Proposition prop:LR in sectors (N={N}, lambda={2*h/J}, {nsamp} sectors) ===")
    x, y = 10, 40                                  # 0-based; A = hop(x,x+1), B = n_y
    a = np.zeros((N, N)); a[x, x + 1] = a[x + 1, x] = 1.0
    b = np.zeros((N, N)); b[y, y] = 1.0
    viol = -np.inf
    sup_avg = 0.0
    sup_bnd = 0.0
    for _ in range(nsamp):
        q = RNG.choice([-1.0, 1.0], size=N)
        E, V = eig_q(q, J, h)
        sv, sb = 0.0, 0.0
        for t in TIMES:
            U = Uof(E, V, t)
            At = U.conj().T @ a @ U
            Z = At @ b - b @ At                       # anti-Hermitian; [dGamma(At),dGamma(b)] = dGamma(Z)
            nrm = dGamma_norm_herm(1j * Z)
            bnd = 2 * 1.0 * 2 * (abs(U[x, y]) + abs(U[x + 1, y]))
            viol = max(viol, nrm - bnd)
            sv = max(sv, nrm); sb = max(sb, bnd)
        sup_avg += sv / nsamp; sup_bnd += sb / nsamp
    P(f"  A = c~_x^+c~_(x+1) + h.c. (x={x+1}), B = n_y (y={y+1}): max(||P_q[A(t),B]|| - monomial bound) = {viol:.2e} (<=0)")
    P(f"  E_q sup_t ||P_q[A(t),B]|| = {sup_avg:.3e};  E_q sup_t (monomial bound) = {sup_bnd:.3e}")
    # clean sector: R(1,N) = 3/(2(N+1)), sup_t |U(1,N)| >= sqrt(R)
    for NN in (20, 50, 100, 200):
        E, V = eig_q(np.ones(NN), J, h)
        R1N = Rmat(V)[0, NN - 1]
        ts = np.linspace(0, 4 * NN, 4000)
        supU = max(abs(Uof(E, V, t)[0, NN - 1]) for t in ts)
        P(f"  clean sector N={NN}: R(1,N) = {R1N:.6f} vs 3/(2(N+1)) = {1.5/(NN+1):.6f};  "
          f"max_t|U(1,N)| on grid = {supU:.4f} >= sqrt(R) = {np.sqrt(R1N):.4f}")


# =========================================================================== D. imbalance
def kappa_all_intervals(Rsum, N):
    """max over intervals B=[a,b] of sum_{x in B, y notin B} R(x,y) = |B| - sum_{x,y in B} R."""
    C = np.zeros((N + 1, N + 1))
    C[1:, 1:] = np.cumsum(np.cumsum(Rsum, 0), 1)
    best = 0.0
    for a_ in range(N):
        for b_ in range(a_, N):
            inner = C[b_ + 1, b_ + 1] - C[a_, b_ + 1] - C[b_ + 1, a_] + C[a_, a_]
            best = max(best, (b_ - a_ + 1) - inner)
    return best


def imbalance_lower_bounds(N, h, J, kap):
    glob = h * h / (N * (J + abs(h)) ** 2)
    blk = max(((N // L) / N) * h * h * max(0.0, np.sqrt(L / 2) - 2 * kap) ** 2 / (L * (J + abs(h)) ** 2)
              for L in range(1, N + 1))
    return glob, blk


def part_D(J=1.0):
    P("=== D. Theorem thm:imbalance ===")
    # D1 identity and pointwise steps
    err_id = 0.0
    viol_pt = -np.inf
    viol_cs = -np.inf
    for N, h in ((64, 0.6), (64, -1.3), (31, 0.4)):
        Gm = np.array([(-1.0) ** (x + 1) for x in range(N)])
        for _ in range(50):
            q = RNG.choice([-1.0, 1.0], size=N)
            E, V = eig_q(q, J, h)
            a = (V ** 2 * Gm[:, None]).sum(0)
            rhs = 2 * h * (V ** 2 * (q * Gm)[:, None]).sum(0)
            err_id = max(err_id, np.max(np.abs(E * a - rhs)))
            if N % 2 == 0:
                R = Rmat(V)
                r = q * Gm
                for L in (4, 9, 20):
                    for a0 in (0, N // 3, N - L):
                        B = np.zeros(N); B[a0:a0 + L] = 1
                        m = (V ** 2 * B[:, None]).sum(0)
                        XB = (E * a * m).sum()
                        XB2 = 2 * h * (r * (R @ B)).sum()
                        err_id = max(err_id, abs(XB - XB2))
                        cross = (R[B == 1][:, B == 0]).sum()
                        viol_pt = max(viol_pt, abs(r[B == 1].sum()) - 2 * cross - abs(XB) / (2 * abs(h)))
                        lhs = XB ** 2
                        rhs2 = (E ** 2 * m).sum() * (a ** 2 * m).sum()
                        viol_cs = max(viol_cs, lhs - rhs2)
                        viol_cs = max(viol_cs, (E ** 2 * m).sum() - L * (2 * J * J + 4 * h * h) + 1e-12 * 0)
    P(f"  identity E_a a_a = 2h<phi_a,V Gamma phi_a> and X_B = 2h sum r_x R(x,B): max err = {err_id:.1e}")
    P(f"  pointwise |sum_B r| - 2 sum_(B x B^c) R - |X_B|/(2|h|): max = {viol_pt:.1e} (<=0)")
    P(f"  Cauchy-Schwarz and sum_a E_a^2 m_a(B) <= L(2J^2+4h^2): max violation = {viol_cs:.1e} (<=0)")
    # D2 time average identity, exact over all sectors, N=8
    N, h = 8, 0.7
    Gm = np.array([(-1.0) ** (x + 1) for x in range(N)])
    S = np.array([1.0 if (x + 1) % 2 == 0 else 0.0 for x in range(N)])
    Tt = 1e7
    lhs, rhs = 0.0, 0.0
    I0 = 0.0
    errI = 0.0
    for q in itertools.product([1.0, -1.0], repeat=N):
        q = np.array(q)
        E, V = eig_q(q, J, h)
        Gab = V.T @ (Gm[:, None] * V)
        om = E[:, None] - E[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            avg = np.where(np.abs(om) > 0, (np.exp(1j * om * Tt) - 1) / (1j * om * Tt), 1.0)
        lhs += (np.abs(Gab) ** 2 * avg).sum().real / N / 2 ** N
        rhs += (np.diag(Gab) ** 2).sum() / N / 2 ** N
        U0 = Uof(E, V, 0.0)
        I0 += (2 / N) * (Gm * (np.abs(U0) ** 2 @ S)).sum() / 2 ** N
        U1 = Uof(E, V, 1.37)
        errI = max(errI, abs((2 / N) * (Gm * (np.abs(U1) ** 2 @ S)).sum() - np.trace(np.diag(Gm) @ U1.conj().T @ np.diag(Gm) @ U1).real / N))
    P(f"  N=8, h=0.7, all 256 sectors: (1/T)int_0^T I(t)dt (T=1e7, exact) = {lhs:.8f}; "
      f"E_q N^-1 sum_a a_a^2 = {rhs:.8f}; I(0) = {I0:.3f}; |I_q(t) - N^-1 Tr(Gamma U* Gamma U)| = {errI:.1e}")
    # h = 0: chiral obstruction
    E, V = eig_q(np.ones(10), J, 0.0)
    P(f"  h=0, N=10: max_a |a_a| = {np.max(np.abs((V**2*np.array([(-1.0)**(x+1) for x in range(10)])[:,None]).sum(0))):.1e} (Ibar = 0)")
    # D3 bounds vs actual, kappa as sup over ALL intervals
    P("  lower bounds (kappa_N = max over all intervals B of E sum_{x in B,y notin B} R; 300 sectors):")
    for N in (32, 64):
        for lam in (0.5, 1.0, 2.0, 4.0, -1.8):
            h = lam * J / 2
            Gm = np.array([(-1.0) ** (x + 1) for x in range(N)])
            Rs = np.zeros((N, N)); Ib = 0.0; ns = 300
            for _ in range(ns):
                q = RNG.choice([-1.0, 1.0], size=N)
                E, V = eig_q(q, J, h)
                Rs += Rmat(V) / ns
                Ib += ((V ** 2 * Gm[:, None]).sum(0) ** 2).sum() / N / ns
            kap = kappa_all_intervals(Rs, N)
            glob, blk = imbalance_lower_bounds(N, h, J, kap)
            P(f"    N={N:3d} lambda={lam:+.1f}: Ibar = {Ib:.4f}, kappa_N = {kap:.2f}, "
              f"global bound = {glob:.2e}, best block bound = {blk:.2e}, Ibar >= max: {Ib >= max(glob, blk)}")
    # wrong (J+h) version at h<0 (erratum)
    N, h = 40, -0.9
    wrong = h * h / (N * (J + h) ** 2)
    P(f"  erratum check: with (J+h) instead of (J+|h|), h=-0.9, N=40, the 'global bound' would be {wrong:.2f} > 1")
    # slow convergence in N at lambda = 0.5
    out = []
    for N in (50, 100, 200, 400):
        h = 0.25
        Gm = np.array([(-1.0) ** (x + 1) for x in range(N)])
        ns = 400 if N <= 200 else 200
        Ib = 0.0
        for _ in range(ns):
            q = RNG.choice([-1.0, 1.0], size=N)
            E, V = eig_q(q, J, h)
            Ib += ((V ** 2 * Gm[:, None]).sum(0) ** 2).sum() / N / ns
        out.append(f"N={N}: {Ib:.4f}")
    P("  lambda=0.5, E Ibar: " + ", ".join(out) + "  (slowly converging in N)")
    # Szarek constant and elementary 1/sqrt(3)
    mins, mins3 = np.inf, np.inf
    for L in range(1, 401):
        k = np.arange(L + 1)
        logp = np.array([math.lgamma(L + 1) - math.lgamma(i + 1) - math.lgamma(L - i + 1) for i in k]) - L * LN2
        EA = float((np.exp(logp) * np.abs(2 * k - L)).sum())
        mins = min(mins, EA / np.sqrt(L / 2)); mins3 = min(mins3, EA / np.sqrt(L / 3))
    P(f"  Khintchine: min_L<=400 E|sum_1^L r_i|/sqrt(L/2) = {mins:.6f} (>=1, Szarek), /sqrt(L/3) = {mins3:.4f}")


# =========================================================================== E. area law
def slater_vector(W, order=None):
    """many-body vector (2^n, JW with modes in 'order' first..) of the Slater determinant of columns of W."""
    n, M = W.shape
    if order is None:
        order = list(range(n))
    Wr = W[order, :]
    vec = np.zeros(2 ** n, complex)
    for occ in itertools.combinations(range(n), M):
        idx = sum(1 << (n - 1 - j) for j in occ)
        vec[idx] = np.linalg.det(Wr[list(occ), :])
    return vec


def part_E(F):
    P("=== E. Theorem thm:arealaw and the Slater lemma ===")
    s = np.linspace(1e-9, 1 - 1e-9, 2000001)
    diff = hfun(s) - 2 * LN2 * np.sqrt(s * (1 - s))
    P(f"  h(nu) <= 2 ln2 sqrt(nu(1-nu)) on 2e6 grid points: max(h - bound) = {diff.max():.2e} (<=0; equality at 1/2)")
    k = np.arange(1, 200001, dtype=float)
    logc = np.array([math.lgamma(2 * kk + 1) - 2 * math.lgamma(kk + 1) - kk * math.log(4) for kk in k]) - np.log(2 * k - 1)
    dk = 1 / (2 * k * (2 * k - 1)) - np.exp(logc) * LN2
    P(f"  d_1 = {dk[0]:.4f} > 0, max_(2<=k<=2e5) d_k = {dk[1:].max():.3e} (<0); "
      f"ratio e_k/c_k decreasing: {bool(np.all(np.diff(1/(2*k*(2*k-1))/np.exp(logc)) < 0))}")
    # Slater lemma: reduced state = rho[C_A], entropy = sum h, N_A ~ Poisson binomial
    err_sp, err_S, err_pN = 0.0, 0.0, 0.0
    for trial in range(10):
        n, M = 8, 4
        Z = RNG.normal(size=(n, M)) + 1j * RNG.normal(size=(n, M))
        W, _ = np.linalg.qr(Z)
        K = W @ W.conj().T
        for Aset in ([0, 1, 2], [2, 3, 4, 5], [1, 4, 6]):
            order = Aset + [j for j in range(n) if j not in Aset]
            vec = slater_vector(W, order)
            la = len(Aset)
            Mm = vec.reshape(2 ** la, -1)
            rho = Mm @ Mm.conj().T
            ev = np.sort(np.linalg.eigvalsh(rho))
            nu = np.linalg.eigvalsh(K[np.ix_(Aset, Aset)])
            prod = np.sort([np.prod([nu[i] if (T >> i) & 1 else 1 - nu[i] for i in range(la)]) for T in range(2 ** la)])
            err_sp = max(err_sp, np.max(np.abs(ev - prod)))
            err_S = max(err_S, abs(vn_entropy_rho(rho) - hfun(nu).sum()))
            NA = np.array([bin(i).count('1') for i in range(2 ** la)])
            pN = np.array([np.real(rho.diagonal()[NA == kk].sum()) for kk in range(la + 1)])
            err_pN = max(err_pN, np.max(np.abs(pN - poibin(nu))))
    P(f"  Slater lemma (n=8, M=4, left/interior/non-contiguous A): spectrum err {err_sp:.1e}, "
      f"entropy err {err_S:.1e}, N_A law err {err_pN:.1e}")
    # chain of inequalities along the dynamics
    N, J, h = 60, 1.0, 1.0
    l = N // 2
    S = [x for x in range(N) if (x + 1) % 2 == 0]
    v1 = v2 = v3 = -np.inf
    supS = []
    for _ in range(200):
        q = RNG.choice([-1.0, 1.0], size=N)
        E, V = eig_q(q, J, h)
        sS = 0.0
        for t in TIMES:
            U = Uof(E, V, t)
            K = U[:, S] @ U[:, S].conj().T
            nu = np.linalg.eigvalsh(K[:l, :l])
            Sq = float(hfun(nu).sum())
            tn = np.linalg.svd(K[l:, :l], compute_uv=False).sum()
            ab = np.abs(K[:l, l:]).sum()
            AU = np.abs(U)
            zb = (AU @ AU.T)[:l, l:]
            v1 = max(v1, Sq - 2 * LN2 * tn)
            v2 = max(v2, tn - ab)
            v3 = max(v3, np.max(np.abs(K[:l, l:]) - zb))
            sS = max(sS, Sq)
        supS.append(sS)
    P(f"  S <= 2ln2||P_B K P_A||_1: {v1:.1e};  ||.||_1 <= sum|K_xy|: {v2:.1e};  |K_xy| <= sum_z|U_xz||U_yz|: {v3:.1e} (all <=0)")
    viol = -np.inf
    dpad = np.arange(2000)
    for Ftest in (np.concatenate([F, np.zeros(2000 - len(F))]), np.minimum(1, 3 * np.exp(-0.2 * dpad))):
        tail = np.cumsum(Ftest[::-1])[::-1]
        for x0 in range(0, 60, 7):
            for dd in range(1, 50):
                zz = np.arange(x0 - 900, x0 + dd + 900)
                lhs = Ftest[np.maximum(np.abs(x0 - zz), np.abs(x0 + dd - zz))].sum()
                rhs = 2 * tail[int(math.ceil(dd / 2))]
                viol = max(viol, lhs - rhs)
    P(f"  sum_z F(max(|x-z|,|y-z|)) <= 2 sum_(m>=ceil(d/2)) F(m): max violation = {viol:.1e} (<=0, up to truncation)")
    M, Mp, S1, S1p, Theta = moments(F)
    P(f"  E sup_t S(rho^q_A(t)) (N=60, lambda=2, CDW, l=30) = {np.mean(supS):.3f};  theorem bound 2ln2*Theta(F_emp) = "
      f"{2*LN2*Theta:.3g} <= 8ln2*M_2(F_emp) = {8*LN2*M[2]:.3g}")


# =========================================================================== F. entropy taxonomy, N = 10
def part_F(N=10, J=1.0, h=1.0, l=5, times=(0.5, 1.0, 2.0, 10.0, 100.0, 1e5), parity_even=True):
    occ = [x for x in range(N) if ((x + 1) % 2 == 0) == parity_even]
    P(f"=== F. Proposition prop:entropies (N={N}, h=J={h}, l={l}, CDW on {'even' if parity_even else 'odd'} sites) ===")
    M = len(occ)
    confs = list(itertools.combinations(range(N), M))
    cidx = np.array([sum(1 << (N - 1 - j) for j in c) for c in confs])
    rows = np.array(confs)
    qs = np.array(list(itertools.product([1.0, -1.0], repeat=N)))        # index = bits (1 <-> q=-1), site 1 MSB
    eig = [eig_q(q, J, h) for q in qs]
    nidx = np.arange(2 ** N)
    pc = np.array([bin(i >> (N - l)).count('1') for i in range(2 ** N)])  # N_A of configuration
    viol = dict(cut=-np.inf, G_mix=-np.inf, mix_up=-np.inf, SN=-np.inf)
    P("     t    S_full  S_A'(Phi)  S_G    S_mix   S_f    S_N   #distinct rho^q_A  H_t")
    for t in times:
        Phi = np.zeros((2 ** N, 2 ** N), complex)        # [n, q]
        SG = 0.0
        rhos = []
        for k, q in enumerate(qs):
            E, V = eig[k]
            U = Uof(E, V, t)
            amp = np.linalg.det(U[rows][:, :, occ])
            v = np.zeros(2 ** N, complex)
            v[cidx] = amp
            Phi[:, k] = np.exp(1j * h * t * q.sum()) * v / 2 ** (N / 2)
            Mm = v.reshape(2 ** l, 2 ** (N - l))
            r = Mm @ Mm.conj().T
            rhos.append(r)
            SG += vn_entropy_rho(r) / 2 ** N
        rho_mix = sum(rhos) / 2 ** N
        Smix = vn_entropy_rho(rho_mix)
        keys = {}
        for r in rhos:
            key = tuple(np.round(np.concatenate([r.real.ravel(), r.imag.ravel()]), 9))
            keys[key] = keys.get(key, 0) + 1
        Ht = shannon(np.array(list(keys.values())) / 2 ** N)
        Sdual = vn_entropy_from_matrix(Phi.reshape(2 ** l, 2 ** (N - l), 2 ** l, 2 ** (N - l))
                                       .transpose(0, 2, 1, 3).reshape(4 ** l, -1))
        # original variables: Psi[n, x] = Phi[n, q(n,x)], q-bits = n ^ x ^ (x >> 1)
        X = nidx[None, :]
        Qi = (nidx[:, None] ^ X ^ (X >> 1))
        Psi = Phi[nidx[:, None], Qi]
        Sfull = vn_entropy_from_matrix(Psi.reshape(2 ** l, 2 ** (N - l), 2 ** l, 2 ** (N - l))
                                       .transpose(0, 2, 1, 3).reshape(4 ** l, -1))
        Mf = Psi.reshape(2 ** l, -1)
        rho_f = Mf @ Mf.conj().T
        Sf = vn_entropy_rho(rho_f)
        pN = np.array([(np.abs(Psi) ** 2)[pc == kk].sum() for kk in range(l + 1)])
        SN = shannon(pN)
        viol['cut'] = max(viol['cut'], abs(Sfull - Sdual) - LN2)
        viol['G_mix'] = max(viol['G_mix'], SG - Smix)
        viol['mix_up'] = max(viol['mix_up'], Smix - min(SG + Ht, l * LN2))
        viol['SN'] = max(viol['SN'], SN - min(Sfull, Smix, Sf))
        flag = "  <-- S_G < S_N" if SG < SN else ""
        P(f"  {t:8.1f}  {Sfull:6.3f}  {Sdual:6.3f}  {SG:6.3f}  {Smix:6.3f}  {Sf:6.3f}  {SN:6.3f}   "
          f"{len(keys):5d}             {Ht:5.3f}{flag}")
    P(f"  max(|S_full - S_A'| - ln2) = {viol['cut']:.3f} (<=0);  max(S_G - S_mix) = {viol['G_mix']:.3f} (<=0);  "
      f"max(S_mix - min(S_G+H_t, l ln2)) = {viol['mix_up']:.3f} (<=0);  max(S_N - min(S_full,S_mix,S_f)) = {viol['SN']:.3f} (<=0)")
    P(f"  maxima: 2 l ln2 = {2*l*LN2:.3f} (S_full), l ln2 = {l*LN2:.3f} (S_mix, S_f)")


def part_F_crosscheck(N=4, J=1.0, h=1.3, l=2, t=1.7):
    """same pipeline as part_F against direct ED of the full LGT (N=4)."""
    L = lgt_setup(N, J, h, with_P=False)
    occ = [x for x in range(N) if (x + 1) % 2 == 0]
    nocc = [1 if x in occ else 0 for x in range(N)]
    psi0 = np.zeros(L['D'], complex)
    psi0[mdc.index(nocc + [0] * N)] = 1
    pt = L['V'] @ (np.exp(-1j * L['w'] * t) * (L['V'].conj().T @ psi0))
    T = pt.reshape([2] * (2 * N))
    order = list(range(l)) + list(range(N, N + l)) + list(range(l, N)) + list(range(N + l, 2 * N))
    S_ed = vn_entropy_from_matrix(T.transpose(order).reshape(4 ** l, -1))
    Sf_ed = vn_entropy_rho((lambda Mx: Mx @ Mx.conj().T)(T.reshape(2 ** l, -1)))
    # pipeline
    M = len(occ)
    confs = list(itertools.combinations(range(N), M))
    cidx = np.array([sum(1 << (N - 1 - j) for j in c) for c in confs]); rows = np.array(confs)
    qs = np.array(list(itertools.product([1.0, -1.0], repeat=N)))
    Phi = np.zeros((2 ** N, 2 ** N), complex)
    for k, q in enumerate(qs):
        U = Uof(*eig_q(q, J, h), t)
        v = np.zeros(2 ** N, complex); v[cidx] = np.linalg.det(U[rows][:, :, occ])
        Phi[:, k] = np.exp(1j * h * t * q.sum()) * v / 2 ** (N / 2)
    nidx = np.arange(2 ** N)
    Psi = Phi[nidx[:, None], nidx[:, None] ^ nidx[None, :] ^ (nidx[None, :] >> 1)]
    S_pl = vn_entropy_from_matrix(Psi.reshape(2 ** l, 2 ** (N - l), 2 ** l, 2 ** (N - l)).transpose(0, 2, 1, 3).reshape(4 ** l, -1))
    Sf_pl = vn_entropy_rho((lambda Mx: Mx @ Mx.conj().T)(Psi.reshape(2 ** l, -1)))
    P(f"  cross-check of the N=10 pipeline against ED (N=4): |S_full diff| = {abs(S_ed-S_pl):.1e}, |S_f diff| = {abs(Sf_ed-Sf_pl):.1e}")


# =========================================================================== G. number entropy
def part_G(F, N=64, J=1.0, h=1.0, nsamp=400):
    P("=== G. Theorem thm:numberentropy ===")
    viol = -np.inf
    for _ in range(3000):
        m = RNG.integers(1, 30)
        p = RNG.random(m) ** RNG.integers(1, 6)
        p /= p.sum()
        k = np.arange(m) * RNG.integers(1, 3)
        var = (p * k ** 2).sum() - (p * k).sum() ** 2
        viol = max(viol, shannon(p) - 0.5 * np.log(2 * np.pi * np.e * (var + 1 / 12)))
    for nn in range(1, 60):
        for v in np.linspace(0.01, 0.99, 25):
            p = poibin(np.full(nn, v)); var = nn * v * (1 - v)
            viol = max(viol, shannon(p) - 0.5 * np.log(2 * np.pi * np.e * (var + 1 / 12)))
    P(f"  Massey bound H(X) <= 1/2 ln(2 pi e (Var X + 1/12)) on random/binomial laws: max(H - bound) = {viol:.2e} (<=0)")
    l = N // 2
    S = [x for x in range(N) if (x + 1) % 2 == 0]
    nT = len(TIMES)
    pmix = np.zeros((nT, l + 1))
    EVq = np.zeros(nT)
    supVar = 0.0
    supD2 = 0.0
    viol2 = -np.inf
    for _ in range(nsamp):
        q = RNG.choice([-1.0, 1.0], size=N)
        E, V = eig_q(q, J, h)
        Q = corr(V)
        Y = Q[:l, l:].sum()
        sv = 0.0; sd = 0.0
        for it, t in enumerate(TIMES):
            U = Uof(E, V, t)
            K = U[:, S] @ U[:, S].conj().T
            nu = np.clip(np.linalg.eigvalsh(K[:l, :l]), 0, 1)
            p = poibin(nu)
            pmix[it] += p / nsamp
            varq = (nu * (1 - nu)).sum()
            viol2 = max(viol2, abs(varq - (np.abs(K[:l, l:]) ** 2).sum()))
            EVq[it] += varq / nsamp
            dev = (p * np.arange(l + 1)).sum() - len([s for s in S if s < l])
            sv = max(sv, varq); sd = max(sd, dev ** 2)
            viol2 = max(viol2, abs(dev) - (np.abs(U[:l, l:]) ** 2).sum() - 1e-12)
            viol2 = max(viol2, (np.abs(U[:l, l:]) ** 2).sum() - Y - 1e-12)
        supVar += sv / nsamp; supD2 += sd / nsamp
    SN = np.array([shannon(p) for p in pmix])
    Var = np.array([(p * np.arange(l + 1) ** 2).sum() - (p * np.arange(l + 1)).sum() ** 2 for p in pmix])
    massey_meas = 0.5 * np.log(2 * np.pi * np.e * (Var + 1 / 12))
    M, Mp, S1, S1p, Theta = moments(F)
    K1, K2 = Theta, S1p ** 2
    thm = 0.5 * np.log(2 * np.pi * np.e * (K1 + K2 + 1 / 12))
    C, eta = exp_majorant(F[:40])
    d = np.arange(0, 2000)
    Fx = np.minimum(1.0, C * np.exp(-eta * d))
    _, _, _, S1px, Thx = moments(Fx)
    thmx = 0.5 * np.log(2 * np.pi * np.e * (Thx + S1px ** 2 + 1 / 12))
    P(f"  checks Var_q = sum|K_xy|^2, |<N_A>_q - |S cap A|| <= sum|U_xz|^2 <= sum Q: max violation = {viol2:.1e}")
    P(f"  N={N}, lambda={2*h/J}, CDW, A=[1,{l}], {nsamp} sectors, t in [0,1e5]:")
    P(f"    S_N(t): max_t = {SN.max():.3f}, S_N(1e5) = {SN[-1]:.3f};  Var_mix(N_A): max_t = {Var.max():.3f}")
    P(f"    Massey bound with the MEASURED variance (non-rigorous): max_t = {massey_meas.max():.3f}")
    P(f"    E sup_t Var_q = {supVar:.3f} (<= K_1);  E sup_t (<N_A>_q - |S cap A|)^2 = {supD2:.3f} (<= K_2)")
    P(f"    theorem constant with F_emp (truncated at d<N): K_1 = Theta = {K1:.3g}, K_2 = Sigma'_1^2 = {K2:.3g}, "
      f"bound = {thm:.3f}")
    P(f"    with an exponential majorant F(d) = min(1, {C:.2f} e^(-{eta:.3f} d)) of F_emp, untruncated: "
      f"K_1 = {Thx:.3g}, K_2 = {S1px**2:.3g}, bound = {thmx:.3f}")
    P(f"    per-sector bound sup_t E_q S_N^q <= 1/2 ln(2 pi e (K_1 + 1/12)) = {0.5*np.log(2*np.pi*np.e*(K1+1/12)):.3f}")


# =========================================================================== H. infinite volume: finite propagation
def part_H(J=1.0, h=1.0):
    P("=== H. Theorem thm:infvol: finite propagation speed and convergence M -> infinity ===")
    Mbig = 150
    viol1 = viol2 = -np.inf
    qZ = RNG.choice([-1.0, 1.0], size=2 * Mbig + 1)          # sites -Mbig..Mbig
    Eb, Vb = eig_q(qZ, J, h)
    nh = 2 * J + 2 * abs(h)
    rows = []
    for t in (0.5, 1.0, 2.0, 5.0):
        c = abs(t) * nh
        Ub = Uof(Eb, Vb, t)
        for M in (5, 10, 20, 40):
            qs = qZ[Mbig - M: Mbig + M + 1]
            E, V = eig_q(qs, J, h)
            U = Uof(E, V, t)
            for xi in range(2 * M + 1):
                for yi in range(2 * M + 1):
                    x, y = xi - M, yi - M
                    dlt = abs(U[xi, yi] - Ub[xi + Mbig - M, yi + Mbig - M])
                    n0 = (M + 1 - abs(x)) + (M + 1 - abs(y))
                    bnd = 2 * np.exp(c + n0 * np.log(c) - math.lgamma(n0 + 1))
                    viol1 = max(viol1, dlt - bnd)
                    dd = abs(x - y)
                    viol2 = max(viol2, abs(U[xi, yi]) - np.exp(c + dd * np.log(c) - math.lgamma(dd + 1)))
            rows.append((t, M, abs(U[M, M + 3] - Ub[Mbig, Mbig + 3])))
    P(f"  |e^-ith_M(x,y) - e^-ith_M'(x,y)| <= 2 e^c c^n0/n0!: max violation = {viol1:.1e};  "
      f"|e^-ith(x,y)| <= e^c c^d/d!: max violation = {viol2:.1e} (<=0)")
    P("  |U^(M)(t)_{0,3} - U^(Mbig)(t)_{0,3}| (one sample): " +
      ", ".join(f"t={t},M={M}: {e:.1e}" for (t, M, e) in rows if M in (5, 20)))
    # sector averages converge
    t = 3.0
    Ms = (4, 8, 16, 32, 64)
    acc = np.zeros(len(Ms))
    ns = 2000
    for _ in range(ns):
        qZ = RNG.choice([-1.0, 1.0], size=2 * Ms[-1] + 1)      # one sequence on [-64,64] for all M (coupling)
        for i, M in enumerate(Ms):
            E, V = eig_q(qZ[Ms[-1] - M: Ms[-1] + M + 1], J, h)
            acc[i] += abs(Uof(E, V, t)[M, M + 3]) / ns
    P(f"  E_q |U^(M)_q(t=3)_(0,3)| (coupled sectors, {ns} samples): " + ", ".join(f"M={M}: {a:.5f}" for M, a in zip(Ms, acc)))


# =========================================================================== I. Duhamel sensitivity lemma
def part_I(N=60, J=1.0, h=1.0):
    P("=== I. Lemma lem:duhamel ===")
    viol1 = viol2 = -np.inf
    for _ in range(40):
        q = RNG.choice([-1.0, 1.0], size=N)
        z = RNG.integers(0, N)
        D = RNG.choice(N, size=RNG.integers(1, 6), replace=False)
        qp = q.copy(); qp[D] *= -1
        E, V = eig_q(q, J, h); Ep, Vp = eig_q(qp, J, h)
        Qp = corr(Vp)
        for s in (0.5, 3.0, 20.0, 150.0):
            lhs = np.linalg.norm((Uof(E, V, s) - Uof(Ep, Vp, s))[:, z])
            rr = np.linspace(0, s, 6001)
            Ur = (Vp[D, :] * Vp[z, :][None, :]) @ np.exp(-1j * Ep[:, None] * rr[None, :])
            vals = np.linalg.norm(Ur, axis=0)
            integral = np.trapezoid(vals, rr) if hasattr(np, 'trapezoid') else np.trapz(vals, rr)
            viol1 = max(viol1, lhs - 4 * abs(h) * integral * (1 + 1e-3))
            viol2 = max(viol2, lhs - 4 * abs(h) * s * Qp[D, z].sum())
    P(f"  ||(U_q(s)-U_q'(s))delta_z|| <= 4|h| int_0^s ||1_D U_q'(r) delta_z|| dr: max violation = {viol1:.1e}")
    P(f"  ... <= 4|h| s sum_(x in D) Q_q'(x,z): max violation = {viol2:.1e} (<=0)")


if __name__ == '__main__':
    F = part_A()
    part_B()
    part_C()
    part_D()
    part_E(F)
    part_F_crosscheck()
    part_F(parity_even=True)
    part_F(parity_even=False, times=(1.0,))
    part_G(F)
    part_H()
    part_I()
    P("done")
