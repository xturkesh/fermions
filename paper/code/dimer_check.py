"""
Checks for Section sec:dimer ("Disorder-free localization depends on the initial state").

Conventions: BRIEF.md section 2 / Section sec:model (imported from model_duality_check.py, which builds
H, G_j, f_j, the link Pauli matrices and the duality unitary U from scratch).

Part A  (exact diagonalization of the full LGT, N = 2, 4, 6; adapted from notes/.../code/dimer_check.py)
  A1  dimer product states Psi^dim(n) = |n> x chi^dim(n): odd links |up>_z, even links b_{2k} in the
      sigma^x-eigenstate with eigenvalue xi_{2k} = (-1)^{N_{<=2k}(n)}.  Sector weights from the PHYSICAL
      projectors P_q = prod_j (1 + q_j G_j)/2 (no use of U): weight 2^{-M} on dimer configurations
      (q_{2k-1} = q_{2k}), 0 elsewhere.  And U^+ Psi^dim(n) = 2^{-M/2} sum_eta |n> x |q(eta)>.
  A2  Psi^dim(n) = ct_{s_1}^+ ... ct_{s_m}^+ Psi^dim(0)   (ct = dressed fermion f_j prod_{k>=j} sz_{b_k}).
  A3  general even-link sigma^x eigenvalues xi and Fock n:  q_{2k} = d_k q_{2k-1},
      d_k = (-1)^{n_{2k-1}+n_{2k}} xi_{2k-2} xi_{2k};  q_{2k-1} = s_k y_k with s_k = xi_{2k-2}(-1)^{n_{2k-1}}.
      In particular CDW with all even links |+> gives the ANTI-dimer law (d_k = -1).
  A4  G_{2k-1} G_{2k} Psi^dim(n) = Psi^dim(n)   (every dimer is Z2-neutral).
  A5  dynamics: full ED of the LGT vs the dimer-law average of the one-body formulas
      (single dressed fermion and CDW dimer states; densities n_y, imbalance), h/J = 0.25, 0.75.
      Control: the same with the iid law (should fail).
Part B  (one-body transfer matrices; JSBS conventions T_{v-E,1} = [[v-E,-1],[1,0]])
  B1  criticality (JSBS Def. 1) iff E = +-lambda and 0 < lambda < 1; M, eta_+-, |<e^{2 i eta}>| = |1-2 lambda^2|.
  B2  transfer-matrix bound ||T^z(n,m)|| <= K1 exp(K2 |z-lambda| (n-m)) on random dimer configurations.
Part C  imbalance identity (lambda^2 - E_a^2) a_a = phi_a(N)^2 - phi_a(1)^2 for dimer configurations;
        time-averaged CDW imbalance Ibar_N (dimer law -> 0, iid law > 0).
Part D  deterministic lower bound M_2(T) >= c T (Prop. prop:dimerlower): Plancherel identity and the
        ratio M_2(T)/T for several dimer configurations (random, clean, periodic) on a finite chain.
"""
import itertools
import os
import sys

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh_tridiagonal, expm
from scipy.sparse.linalg import expm_multiply

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from model_duality_check import (physical, build_U, gauge_product_state, index, bits)  # noqa: E402

OUT = open(os.path.join(HERE, 'dimer_check.out'), 'w')


def pr(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    OUT.write(s + '\n')
    OUT.flush()


def dimer_kinds(N, n, xi=None):
    """gauge product state: odd links |up>_z; even link b_{2k} sigma^x eigenvalue xi_{2k}
    (default: xi_{2k} = (-1)^{N_{<=2k}(n)}, the dimer choice).  n is 1-indexed list n[0]=n_1."""
    kinds, xis = [], {}
    for j in range(1, N + 1):
        if j % 2 == 1:
            kinds.append('zu')
        else:
            v = (-1) ** sum(n[:j]) if xi is None else xi[j // 2 - 1]
            xis[j] = v
            kinds.append('x+' if v == 1 else 'x-')
    return kinds, xis


def fock(N, n):
    v = np.zeros(2 ** N)
    v[index(list(n))] = 1.
    return v


def dual_vec(N, n, q):
    vq = np.zeros(2 ** N)
    vq[index([(1 - t) // 2 for t in q])] = 1.
    return np.kron(fock(N, n), vq)


def q_of_eta(eta):
    return [eta[j // 2] for j in range(2 * len(eta))]


def hq(N, J, h, q):
    A = np.zeros((N, N))
    for j in range(N - 1):
        A[j, j + 1] = A[j + 1, j] = -J
    return A + np.diag(2 * h * np.asarray(q, float))


def sector_weights_physical(P, psi, N):
    """||P_q psi||^2 with P_q = prod_j (1 + q_j G_j)/2 built from the physical G_j."""
    Gd = [P['G'][j].toarray() if N <= 4 else P['G'][j] for j in range(N)]
    out = {}
    for q in itertools.product([1, -1], repeat=N):
        v = psi.copy()
        for j in range(N):
            v = 0.5 * (v + q[j] * (Gd[j] @ v))
        out[q] = float(np.vdot(v, v).real)
    return out


def partA(rng):
    pr('=== Part A: ED of the LGT for dimer product states ===')
    J = 1.0
    for N in [2, 4, 6]:
        M = N // 2
        P = physical(N, J, 0.3)
        U = build_U(N)
        Ud = U.conj().T
        dims = []
        # list of matter Fock states: vacuum, CDW (even sites), each single fermion, 3 random
        tests = [('vacuum', [0] * N), ('CDW', [1 if (j + 1) % 2 == 0 else 0 for j in range(N)])]
        tests += [(f'single@{x}', [1 if j == x - 1 else 0 for j in range(N)]) for x in range(1, N + 1)]
        tests += [(f'random{r}', list(rng.integers(0, 2, size=N))) for r in range(3)]
        err_w = err_U = err_neutral = 0.
        for name, n in tests:
            kinds, _ = dimer_kinds(N, n)
            psi = np.kron(fock(N, n), gauge_product_state(N, kinds))
            # A1: physical sector weights
            w = sector_weights_physical(P, psi, N)
            for q, val in w.items():
                isdimer = all(q[2 * k] == q[2 * k + 1] for k in range(M))
                err_w = max(err_w, abs(val - (2.0 ** (-M) if isdimer else 0.)))
            # A1: dual vector
            Phi = Ud @ psi
            pred = sum(2 ** (-M / 2) * dual_vec(N, n, q_of_eta(eta))
                       for eta in itertools.product([1, -1], repeat=M))
            err_U = max(err_U, float(np.max(np.abs(Phi - pred))))
            # A4: dimer neutrality
            for k in range(M):
                GG = P['G'][2 * k] @ P['G'][2 * k + 1]
                err_neutral = max(err_neutral, float(np.max(np.abs(GG @ psi - psi))))
        pr(f'N={N}: A1 max|weight - 2^-M 1[dimer]| = {err_w:.1e} (physical P_q);  '
           f'max|U^+Psi - 2^(-M/2) sum_eta |n,q(eta)>| = {err_U:.1e};  A4 max|G_(2k-1)G_(2k)Psi - Psi| = {err_neutral:.1e}')
        # A2: dressed-fermion creation over the dimer vacuum
        ct = []
        for x in range(N):
            S = P['Id']
            for k in range(x, N):
                S = S @ P['sz'][k]
            ct.append((P['f'][x] @ S).tocsr())
        kinds0, _ = dimer_kinds(N, [0] * N)
        vac = np.kron(fock(N, [0] * N), gauge_product_state(N, kinds0))
        err_c = 0.
        for nidx in range(2 ** N):
            n = bits(nidx, N)
            kinds, _ = dimer_kinds(N, n)
            target = np.kron(fock(N, n), gauge_product_state(N, kinds))
            v = vac.copy()
            for s in reversed([j for j in range(N) if n[j] == 1]):
                v = ct[s].T @ v
            err_c = max(err_c, float(np.max(np.abs(v - target))))
        pr(f'N={N}: A2 max over all n of |ct^+_(s1)...ct^+_(sm) Psi^dim(0) - Psi^dim(n)| = {err_c:.1e}')
        # A3: general xi
        err_d = 0.
        anti = None
        for trial in range(6):
            n = list(rng.integers(0, 2, size=N))
            xi = list(rng.choice([1, -1], size=M))
            if trial == 0:   # CDW with all even links |+>: anti-dimer law expected
                n = [1 if (j + 1) % 2 == 0 else 0 for j in range(N)]
                xi = [1] * M
            kinds, xis = dimer_kinds(N, n, xi)
            psi = np.kron(fock(N, n), gauge_product_state(N, kinds))
            w = sector_weights_physical(P, psi, N)
            xi_full = {0: 1}
            xi_full.update({2 * (k + 1): xi[k] for k in range(M)})
            d = [(-1) ** (n[2 * k] + n[2 * k + 1]) * xi_full[2 * k] * xi_full[2 * k + 2] for k in range(M)]
            for q, val in w.items():
                ok = all(q[2 * k + 1] == d[k] * q[2 * k] for k in range(M))
                err_d = max(err_d, abs(val - (2.0 ** (-M) if ok else 0.)))
            if trial == 0:
                anti = d
        pr(f'N={N}: A3 max|weight - 2^-M 1[q_(2k)=d_k q_(2k-1)]| over 6 (n,xi) = {err_d:.1e};  '
           f'CDW with all even links |+>: d = {anti} (anti-dimers)')
    # A5: dynamics
    pr('--- A5: dynamics, full ED of the LGT vs dimer-law average of one-body formulas ---')
    for N in [4, 6]:
        M = N // 2
        for h in [0.25, 0.75]:
            P = physical(N, J, h)
            Hs = P['H'].tocsc()
            nd = [P['n'][x].diagonal() for x in range(N)]
            G = np.array([(-1.) ** (x + 1) for x in range(N)])
            res = []
            for name, n in [('single@%d' % (N // 2), [1 if j == N // 2 - 1 else 0 for j in range(N)]),
                            ('CDW', [1 if (j + 1) % 2 == 0 else 0 for j in range(N)])]:
                kinds, _ = dimer_kinds(N, n)
                psi0 = np.kron(fock(N, n), gauge_product_state(N, kinds)).astype(complex)
                occ = [j for j in range(N) if n[j] == 1]
                times = [0.5, 1.7, 4.0, 11.0]
                err_dim = err_iid = 0.
                for t in times:
                    psi = expm_multiply(-1j * t * Hs, psi0)
                    ed = np.array([np.vdot(psi, nd[x] * psi).real for x in range(N)])
                    dim = np.zeros(N)
                    for eta in itertools.product([1, -1], repeat=M):
                        Ut = expm(-1j * t * hq(N, J, h, q_of_eta(eta)))
                        dim += (np.abs(Ut[:, occ]) ** 2).sum(1) / 2 ** M
                    iid = np.zeros(N)
                    for q in itertools.product([1, -1], repeat=N):
                        Ut = expm(-1j * t * hq(N, J, h, q))
                        iid += (np.abs(Ut[:, occ]) ** 2).sum(1) / 2 ** N
                    err_dim = max(err_dim, float(np.max(np.abs(ed - dim))))
                    err_iid = max(err_iid, float(np.max(np.abs(ed - iid))))
                    if name == 'CDW':
                        imb_ed = 2 / N * float(G @ ed)
                        imb_dim = 2 / N * float(G @ dim)
                        err_dim = max(err_dim, abs(imb_ed - imb_dim))
                res.append(f'{name}: dimer-law err = {err_dim:.1e}, iid-law control err = {err_iid:.2f}')
            pr(f'N={N}, h/J={h}: ' + ';  '.join(res))


def partB():
    pr('=== Part B: transfer matrices of the dimer model (JSBS conventions, t = 1) ===')

    def A(a):
        return np.array([[a, -1.], [1., 0.]], dtype=complex)

    worst_crit = 0.
    lines = []
    for lam in [0.1, 0.25, 0.3, 0.5, 1 / np.sqrt(2), 0.8, 0.95, 0.999, 1.0, 1.2, np.sqrt(2), 2.0, 3.0]:
        Tp = A(lam - lam) @ A(lam - lam)        # (+,+) dimer at E = lambda
        Tm = A(-lam - lam) @ A(-lam - lam)      # (-,-) dimer at E = lambda
        comm = np.max(np.abs(Tp @ Tm - Tm @ Tp))
        trm = np.trace(Tm).real
        crit = (np.max(np.abs(Tp + np.eye(2))) < 1e-14) and (abs(trm) < 2) and comm < 1e-14
        if lam < 1:
            th = np.arccos(-lam)
            S = np.array([[np.sin(th), np.cos(th)], [0., 1.]])
            Mm = np.linalg.inv(S)
            R = lambda a: np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
            e1 = np.max(np.abs(Mm @ A(-2 * lam) @ S - R(th)))
            e2 = np.max(np.abs(Mm @ Tm @ S - R(2 * th)))
            e3 = np.max(np.abs(Mm @ Tp @ S - R(np.pi)))
            kappa = np.linalg.cond(S)
            # independent: eigenvalues of Tm are e^{+- i eta_-}
            ev = np.linalg.eigvals(Tm)
            avg = 0.5 * (1 + ev[0] ** 2)        # <e^{2 i eta}> with eta_+ = pi, e^{i eta_-} = ev[0]
            avg2 = 0.5 * (1 + ev[1] ** 2)
            worst_crit = max(worst_crit, e1, e2, e3)
            lines.append(f'  lambda={lam:.4f}: critical={crit}, tr T_-={trm:+.4f}, [T+,T-]={comm:.0e}, '
                         f'|<e^(2i eta)>| = {abs(avg):.4f} / {abs(avg2):.4f} (formula |1-2l^2| = {abs(1 - 2 * lam ** 2):.4f}), '
                         f'|<e^(4i eta)> - 1| = {abs(0.5 * (1 + ev[0] ** 4) - 1):.2e}, kappa(M) = {kappa:.4f} '
                         f'(formula {np.sqrt((1 + lam) / (1 - lam)):.4f})')
        else:
            lines.append(f'  lambda={lam:.4f}: critical={crit}, tr T_-={trm:+.4f} ({"parabolic" if abs(abs(trm) - 2) < 1e-12 else "hyperbolic"})')
    for l in lines:
        pr(l)
    pr(f'  max error of M A(-2l) M^-1 = R_theta, M T_- M^-1 = R_2theta, M T_+ M^-1 = R_pi: {worst_crit:.1e}')
    # B1': for E not in {+-lambda}, the dimer matrices do not commute
    mn = np.inf
    for lam in [0.3, 0.7, 1.5]:
        for E in np.linspace(-4, 4, 8001):
            if min(abs(E - lam), abs(E + lam)) < 1e-3:
                continue
            Tp = A(lam - E) @ A(lam - E)
            Tm = A(-lam - E) @ A(-lam - E)
            c = np.max(np.abs(Tp @ Tm - Tm @ Tp))
            # exact: [A(a)^2, A(b)^2] = -ab(a-b) sigma_x  with a = lam-E, b = -lam-E
            mn = min(mn, c / abs((lam - E) * (-lam - E) * 2 * lam))
    pr(f'  B1: min over E (|E -+ lambda| >= 1e-3) of |[T+,T-]| / |ab(a-b)| = {mn:.6f} (exact value 1)')
    # anti-dimer polymers (lam,-lam) and (-lam,lam) (sector law with d_k = -1): commute at no energy
    mnA = np.inf
    for lam in [0.3, 0.7, 1.5]:
        for E in np.linspace(-4, 4, 8001):
            T1 = A(-lam - E) @ A(lam - E)
            T2 = A(lam - E) @ A(-lam - E)
            mnA = min(mnA, float(np.max(np.abs(T1 @ T2 - T2 @ T1))))
    pr(f'  B1 (anti-dimer law): min over lambda in {{0.3,0.7,1.5}}, E in [-4,4] of max|[T_(l,-l), T_(-l,l)]| = {mnA:.3e} (> 0)')
    # B2: transfer-matrix bound on random configurations
    rng = np.random.default_rng(7)
    worst = 0.
    for lam in [0.25, 0.5, 0.9]:
        th = np.arccos(-lam)
        kappa = np.sqrt((1 + lam) / (1 - lam))
        K1 = kappa * (2 + 2 * lam) ** 2
        K2 = kappa * (3 + 4 * lam) / 2
        for trial in range(40):
            L = 400
            r = rng.integers(0, 2)
            eta = rng.choice([-1, 1], size=L // 2 + 2)
            v = np.array([lam * eta[(n + r) // 2] for n in range(L)])
            dz = rng.uniform(0, 1) * np.exp(2j * np.pi * rng.uniform())
            dz *= 10 ** rng.uniform(-3, 0)
            z = lam + dz
            m = rng.integers(0, L // 2)
            T = np.eye(2, dtype=complex)
            for n in range(m, L):
                T = A(v[n] - z) @ T
                logratio = np.log(np.linalg.norm(T, 2)) - np.log(K1) - K2 * abs(dz) * (n + 1 - m)
                worst = max(worst, float(np.exp(min(logratio, 50.))))
    pr(f'  B2: max over 120 random (config, alignment, z, m, n) of ||T^z(n,m)|| / (K1 e^(K2|z-l|(n-m))) = {worst:.4f} (must be <= 1)')


def partC(rng):
    pr('=== Part C: imbalance identity for dimer configurations ===')
    err = 0.
    for lam in [0.3, 0.5, 1.0, 1.5, 3.0]:
        for trial in range(20):
            N = 200
            q = np.repeat(rng.choice([-1, 1], size=N // 2), 2).astype(float)
            E, V = eigh_tridiagonal(lam * q, -np.ones(N - 1))
            G = (-1.) ** np.arange(1, N + 1)
            a = (G[:, None] * V ** 2).sum(0)
            lhs = (lam ** 2 - E ** 2) * a
            rhs = V[-1] ** 2 - V[0] ** 2
            err = max(err, float(np.max(np.abs(lhs - rhs))))
    pr(f'  C1: max |(l^2 - E_a^2) a_a - (phi_a(N)^2 - phi_a(1)^2)| over 100 configurations, N=200: {err:.1e}')
    errM = 0.
    for _ in range(2000):
        w, E = 3 * rng.normal(size=2)
        A1 = np.array([[w, -1.], [1., 0.]])
        A2 = A1 @ A1
        B = np.array([[1., E], [E, 1.]])
        Sw = np.array([[w * w - 1., -w], [-w, 1.]])
        errM = max(errM, float(np.max(np.abs(A2.T @ B @ A2 - B - w * (w + 2 * E) * Sw))) / (1 + w ** 4 + E ** 2))
    pr(f'  C0: matrix identity (A(w)^2)^T B_E A(w)^2 - B_E = w(w+2E) S(w), 2000 random (w,E): max rel. error {errM:.1e}')
    for lam in [0.5, 1.5]:
        row = []
        for N in [100, 200, 400, 800]:
            S = 400 if N <= 200 else 100
            vd, vi = [], []
            for s in range(S):
                for law, store in [('dimer', vd), ('iid', vi)]:
                    if law == 'dimer':
                        q = np.repeat(rng.choice([-1, 1], size=N // 2), 2).astype(float)
                    else:
                        q = rng.choice([-1, 1], size=N).astype(float)
                    E, V = eigh_tridiagonal(lam * q, -np.ones(N - 1))
                    G = (-1.) ** np.arange(1, N + 1)
                    a = (G[:, None] * V ** 2).sum(0)
                    store.append((a ** 2).sum() / N)
            row.append(f'N={N}: dimer {np.mean(vd):.2e}+-{np.std(vd) / np.sqrt(S):.0e} (N*Ibar={N * np.mean(vd):.3f}), '
                       f'iid {np.mean(vi):.4f}+-{np.std(vi) / np.sqrt(S):.4f}')
        pr(f'  C2: lambda={lam}: time-averaged CDW imbalance Ibar_N = E (1/N) sum_a a_a^2:')
        for r_ in row:
            pr('      ' + r_)


def partD():
    pr('=== Part D: deterministic lower bound M_2(T) >= c T (finite chain, x in the middle) ===')
    # Plancherel identity: M_p(T) = (1/(2 pi T)) sum_n |n-x|^p int |G^{E+i/(2T)}(n)|^2 dE, small chain
    lam, N, x, T, p = 0.5, 40, 20, 7.0, 2.0
    rng = np.random.default_rng(11)
    q = np.repeat(rng.choice([-1, 1], size=N // 2), 2).astype(float)
    E, V = eigh_tridiagonal(lam * q, -np.ones(N - 1))
    c = V[x - 1]
    w = (np.arange(1, N + 1) - x) ** p
    K = 1.0 / (1.0 + (E[:, None] - E[None, :]) ** 2 * T ** 2)
    A_ = V * c[None, :]
    MT_exact = float(np.einsum('x,xa,ab,xb->', w, A_, K, A_))
    # quadrature over the whole real line with E = 3 tan(u)
    u = np.linspace(-np.pi / 2, np.pi / 2, 400001)[1:-1]
    Es = 3 * np.tan(u)
    jac = 3 / np.cos(u) ** 2
    eps = 1 / (2 * T)
    Gn = V @ (c[:, None] / (E[:, None] - (Es[None, :] + 1j * eps)))     # G^{E+i eps}(n), shape (N, nE)
    integrand = (w[:, None] * np.abs(Gn) ** 2).sum(0) * jac
    MT_res = float(np.trapezoid(integrand, u)) / (2 * np.pi * T)
    pr(f'  D0: Plancherel identity, N={N}, T={T}: exact eigen-formula {MT_exact:.8f}, resolvent quadrature {MT_res:.8f}')
    # lower bound constant
    for lam in [0.25, 0.5]:
        kappa = np.sqrt((1 + lam) / (1 - lam))
        K1 = kappa * (2 + 2 * lam) ** 2
        K2 = kappa * (3 + 4 * lam) / 2
        c0 = 1 / (3 + 2 * lam)
        c1 = c0 / (K1 * np.exp(K2))
        cl = c1 ** 2 / (2 ** (p + 3) * np.pi * (p + 1))
        N = 600
        x = N // 2 + 1  # odd site
        configs = {}
        configs['clean (+)'] = np.ones(N)
        configs['alternating dimers'] = np.array([(-1) ** ((j // 2) % 2) for j in range(N)], float)
        configs['period-6 (++,++,--)'] = np.array([1 if (j // 2) % 3 < 2 else -1 for j in range(N)], float)
        for s in range(3):
            configs[f'random {s}'] = np.repeat(rng.choice([-1, 1], size=N // 2), 2).astype(float)
        Ts = [1, 3, 10, 30, 100, 250]
        pr(f'  lambda={lam}: explicit constant c_(lambda,2) = {cl:.3e};  ratio M_2(T)/T for T = {Ts}:')
        for name, qq in configs.items():
            E, V = eigh_tridiagonal(lam * qq, -np.ones(N - 1))
            c = V[x - 1]
            wv = (np.arange(1, N + 1) - x) ** 2.
            A_ = V * c[None, :]
            B = (A_ * wv[:, None]).T @ A_          # B_ab = sum_x w_x A_xa A_xb
            rat = []
            for T in Ts:
                K = 1.0 / (1.0 + (E[:, None] - E[None, :]) ** 2 * T ** 2)
                rat.append(float(np.sum(B * K)) / T)
            pr(f'      {name:22s}: ' + ' '.join(f'{r_:.3e}' for r_ in rat) + f'   (min ratio / c = {min(rat) / cl:.1e})')


if __name__ == '__main__':
    rng = np.random.default_rng(2026)
    partA(rng)
    partB()
    partC(rng)
    partD()
    OUT.close()
