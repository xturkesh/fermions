"""
Numerical checks for Section "Number entropy in U(1) gauge theories: kinematic bounds" (sections/u1.tex).

Conventions (as in the text): sites j=1..L, staggered fermions, rho_j = n_j - (1-(-1)^j)/2, links b_j=(j,j+1),
j=1..L (b_L carries no hopping), fixed boundary flux E_0 = e0, Gauss law in the sector Q:
        E_j - E_{j-1} = rho_j + Q_j       (j = 1..L),
so that in a fixed sector the electric fields are functions of the occupations ("gauge fields integrated out").
Quantum link: |E_j| <= S (configurations violating it are absent); U_j = S^+_j raises E_j by one.
Wilson link: E_j in Z, unbounded; U_j = e^{i phi_j}.

A. Prop. prop:gausscap: the identity N_A = E_b - E_{a-1} + #odd(A) - sum_A Q on every basis state of a
   sector, the number of values of N_A (<= 2S+1 at the edge a=1, <= 4S+1 in the bulk), attainment of
   ln(#values), and the bound along gauge-invariant dynamics with random couplings and interactions.
A'. Superposition over background-charge sectors (initial state of Brenes et al. type: Neel/bare vacuum
   fermions, every link in an equal superposition of its electric eigenstates):
   S_N(t) <= ln(2S+1) + H(sum_A Q) (edge) and <= ln(4S+1) + H(sum_A Q) (bulk).
B. Prop. prop:wilson: Wilson links, single sector, S_N([1,l]) = H(E_l); energy bound on the cut average of
   <E_l^2>_t and of S_N([1,l])(t), uniformly in t.
C. Remark rem:randomwalk: standard (single-link) Z2 LGT -> one fermion in V(x) = 2h sum_{b>=x} Q_b
   (sector-resolved ED), rank-alternating potential for two fermions (not free), and eigenfunction
   correlators: random-walk potential vs i.i.d. Bernoulli.
"""
import itertools, time
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh, expm, eigh_tridiagonal

rng = np.random.default_rng(11)
T0 = time.time()


def H_shannon(p):
    p = np.asarray(p, float); p = p[p > 1e-15]
    return float(-(p * np.log(p)).sum())


def odd(j):
    return j % 2 == 1


# ----------------------------------------------------------------------------------------------
# QLM / Wilson sector: enumerate allowed occupation configurations
# ----------------------------------------------------------------------------------------------
def fields(n, Q, e0):
    """n, Q: arrays indexed 1..L stored at 0..L-1. returns E_1..E_L"""
    L = len(n)
    rho = np.array([n[j - 1] - (1 if odd(j) else 0) for j in range(1, L + 1)])
    return e0 + np.cumsum(rho + np.asarray(Q))


def sector_basis(L, Q, e0, S=None, Nf=None):
    confs = []
    for n in itertools.product([0, 1], repeat=L):
        if Nf is not None and sum(n) != Nf:
            continue
        E = fields(n, Q, e0)
        if S is not None and np.any(np.abs(E) > S + 1e-9):
            continue
        confs.append((n, E))
    return confs


def sector_H(confs, L, w, m, g2, S=None, theta=0.0, vrand=None, Vnn=0.0):
    """H = -w sum_{j<L} (f_j^dag U_j f_{j+1} + h.c.) + m sum (-1)^j n_j + (g2/2) sum (E_j+theta)^2
       + sum_j v_j n_j + Vnn sum n_j n_{j+1}   (all gauge invariant) in the integrated representation.
       QLM: f_j^dag S^+_j f_{j+1} has the matrix element sqrt(S(S+1)-E(E+1)) (E = E_j before the hop),
       normalised by sqrt(S(S+1)) so that every bond term has norm <= w.  Wilson: matrix element 1.
       Jordan-Wigner signs are trivial (nearest-neighbour hops between consecutive modes)."""
    idx = {c[0]: k for k, c in enumerate(confs)}
    D = len(confs); Hm = np.zeros((D, D))
    for k, (n, E) in enumerate(confs):
        diag = sum(m * (-1) ** j * n[j - 1] for j in range(1, L + 1)) + 0.5 * g2 * np.sum((E + theta) ** 2)
        if vrand is not None:
            diag += float(np.dot(vrand, n))
        diag += Vnn * sum(n[j] * n[j + 1] for j in range(L - 1))
        Hm[k, k] = diag
        for j in range(1, L):                      # bond (j, j+1), link b_j between them
            if n[j] == 1 and n[j - 1] == 0:        # f_j^dag U_j f_{j+1}: particle j+1 -> j, E_j -> E_j + 1
                nn = list(n); nn[j] = 0; nn[j - 1] = 1; nn = tuple(nn)
                if nn in idx:
                    Ej = E[j - 1]
                    amp = 1.0 if S is None else np.sqrt(S * (S + 1) - Ej * (Ej + 1)) / np.sqrt(S * (S + 1))
                    k2 = idx[nn]
                    Hm[k2, k] += -w * amp; Hm[k, k2] += -w * amp
    return Hm


print('=' * 100)
print('A. Quantum link models: Gauss-law caps on the number entropy (Prop. prop:gausscap)')
print('=' * 100)
L = 10
for S in (0.5, 1.0, 1.5):
    e0 = -0.5 if S == 0.5 else (0.0 if S == 1.0 else 0.5)
    for trial in range(3):
        Q = np.zeros(L, int) if trial == 0 else rng.choice([-1, 0, 1], size=L)
        confs = sector_basis(L, Q, e0, S=S)
        if len(confs) == 0:
            continue
        for (a, b) in [(1, 5), (4, 7)]:
            A = list(range(a, b + 1))
            # identity on every basis state
            err = 0
            vals = set()
            for (n, E) in confs:
                NA = sum(n[j - 1] for j in A)
                Eam1 = e0 if a == 1 else E[a - 2]
                rhs = E[b - 1] - Eam1 + sum(1 for j in A if odd(j)) - sum(Q[j - 1] for j in A)
                err = max(err, abs(NA - rhs)); vals.add(NA)
            cap = 2 * S + 1 if a == 1 else 4 * S + 1
            # random state
            psi = rng.normal(size=len(confs)) + 1j * rng.normal(size=len(confs)); psi /= np.linalg.norm(psi)
            NAs = np.array([sum(n[j - 1] for j in A) for (n, E) in confs])
            p = [np.sum(np.abs(psi[NAs == v]) ** 2) for v in sorted(vals)]
            # maximal-entropy superposition: one basis state per value, equal weights
            pmax = [1.0 / len(vals)] * len(vals)
            print(f'S={S} Q={"0" if trial == 0 else "random"} A=[{a},{b}]: sector dim={len(confs)}, identity err={err:.0e}, '
                  f'#values of N_A={len(vals)} (cap {int(cap)}), S_N(random state)={H_shannon(p):.3f}, '
                  f'S_N(equal-weight superposition)={H_shannon(pmax):.3f} <= ln cap={np.log(cap):.3f}')

# dynamics under random gauge-invariant Hamiltonians (with interactions), S = 1, bulk and edge intervals
print('Dynamics (S=1, L=10, Q=0, e0=0, random on-site potential and n n interaction, Fock initial state):')
S = 1.0; e0 = 0.0; Q = np.zeros(L, int)
confs = sector_basis(L, Q, e0, S=S, Nf=L // 2)
Hm = sector_H(confs, L, w=1.0, m=0.3, g2=0.5, S=S, vrand=rng.uniform(-1, 1, size=L), Vnn=0.7)
E_, V_ = np.linalg.eigh(Hm)
k0 = [k for k, (n, E) in enumerate(confs) if n == tuple(1 if odd(j) else 0 for j in range(1, L + 1))][0]
c0 = V_[k0, :]
for (a, b) in [(1, 5), (4, 7)]:
    A = list(range(a, b + 1)); NAs = np.array([sum(n[j - 1] for j in A) for (n, E) in confs])
    cap = 2 * S + 1 if a == 1 else 4 * S + 1
    mx = 0
    for t in np.concatenate([np.linspace(0, 20, 41), np.logspace(1.5, 4, 20)]):
        psi = V_ @ (np.exp(-1j * t * E_) * c0)
        p = [np.sum(np.abs(psi[NAs == v]) ** 2) for v in np.unique(NAs)]
        mx = max(mx, H_shannon(p))
    print(f'   A=[{a},{b}]: sup_t S_N(t) = {mx:.3f} <= ln({int(cap)}) = {np.log(cap):.3f}')
print(f'[time {time.time() - T0:.1f}s]')

# ----------------------------------------------------------------------------------------------
# A'. superpositions over background-charge sectors
# ----------------------------------------------------------------------------------------------
print()
print("A'. Superposition over background-charge sectors (bare-vacuum fermions, links in equal superposition)")
for (S, L, intervals) in [(0.5, 8, [(1, 4), (3, 6)]), (1.0, 6, [(1, 3), (2, 5)])]:
    e0 = 0.5 if S == 0.5 else 0.0
    n0 = tuple(1 if odd(j) else 0 for j in range(1, L + 1))       # bare vacuum: rho = 0 everywhere
    Evals = np.arange(-S, S + 1)
    for (a, b) in intervals:
        A = list(range(a, b + 1))
        cap = 2 * S + 1 if a == 1 else 4 * S + 1
        # each link configuration E in {-S..S}^L defines one sector Q (Q_j = E_j - E_{j-1} - rho_j(n0))
        sectors = {}
        for E in itertools.product(Evals, repeat=L):
            E = np.array(E); Eprev = np.concatenate([[e0], E[:-1]])
            Qv = tuple(int(round(x)) for x in (E - Eprev))         # rho(n0) = 0
            sectors[Qv] = sectors.get(Qv, 0) + 1
        weight = {k: v / (len(Evals) ** L) for k, v in sectors.items()}
        assert all(v == 1 for v in sectors.values())
        HQ = H_shannon([sum(wq for Qv, wq in weight.items() if sum(Qv[j - 1] for j in A) == s)
                        for s in set(sum(Qv[j - 1] for j in A) for Qv in weight)])
        ts = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0]
        pt = [dict() for _ in ts]
        for Qv, wq in weight.items():
            confs = sector_basis(L, Qv, e0, S=S, Nf=sum(n0))
            k0 = [k for k, (n, E) in enumerate(confs) if n == n0][0]
            Hm = sector_H(confs, L, w=1.0, m=0.0, g2=0.4, S=S)
            Ee, Vv = np.linalg.eigh(Hm); c0 = Vv[k0, :]
            NAs = np.array([sum(n[j - 1] for j in A) for (n, E) in confs])
            for it, t in enumerate(ts):
                psi = Vv @ (np.exp(-1j * t * Ee) * c0)
                for v in np.unique(NAs):
                    pt[it][v] = pt[it].get(v, 0.0) + wq * np.sum(np.abs(psi[NAs == v]) ** 2)
        SN = [H_shannon(list(d.values())) for d in pt]
        print(f'S={S} A=[{a},{b}]: {len(weight)} sectors; H(sum_A Q)={HQ:.3f}; S_N(t) for t={ts}: '
              + ' '.join(f'{x:.3f}' for x in SN) + f'  <= ln({int(cap)}) + H(sum_A Q) = {np.log(cap) + HQ:.3f}')
print(f'[time {time.time() - T0:.1f}s]')

# ----------------------------------------------------------------------------------------------
# B. Wilson links
# ----------------------------------------------------------------------------------------------
print()
print('=' * 100)
print('B. Wilson links (lattice Schwinger model, gauge field integrated out in the sector Q=0, e0=0)')
print('=' * 100)
L = 12
for (w, m, g2) in [(1.0, 0.0, 0.25), (1.0, 0.5, 0.5), (1.0, 0.0, 2.0)]:
    Q = np.zeros(L, int); e0 = 0
    n0 = tuple(1 if odd(j) else 0 for j in range(1, L + 1))
    confs = sector_basis(L, Q, e0, S=None, Nf=sum(n0))
    Hm = sector_H(confs, L, w=w, m=m, g2=g2, S=None)
    Ee, Vv = np.linalg.eigh(Hm)
    k0 = [k for k, (n, E) in enumerate(confs) if n == n0][0]; c0 = Vv[k0, :]
    Earr = np.array([E for (n, E) in confs])                           # D x L, E_1..E_L
    H0 = Hm[k0, k0]
    bound_E2 = (2.0 / g2) * (H0 / L + w + abs(m))
    bound_SN = 0.5 * np.log(2 * np.pi * np.e * (bound_E2 + 1.0 / 12))
    worst_E2 = 0; worst_SN = 0; worst_single = 0; ident = 0
    for t in np.concatenate([np.linspace(0, 30, 61), np.logspace(1.6, 4, 15)]):
        psi = Vv @ (np.exp(-1j * t * Ee) * c0); P = np.abs(psi) ** 2
        E2 = (P @ Earr ** 2).mean()
        SNs = []
        for l in range(1, L + 1):
            vals, inv = np.unique(Earr[:, l - 1], return_inverse=True)
            pE = np.bincount(inv, weights=P)
            Nl = np.array([sum(n[:l]) for (n, E) in confs])
            vN, invN = np.unique(Nl, return_inverse=True)
            pN = np.bincount(invN, weights=P)
            ident = max(ident, abs(H_shannon(pE) - H_shannon(pN)))
            SNs.append(H_shannon(pN))
        worst_E2 = max(worst_E2, E2); worst_SN = max(worst_SN, np.mean(SNs)); worst_single = max(worst_single, max(SNs))
    print(f'w={w} m={m} g^2={g2}: |S_N([1,l]) - H(E_l)| <= {ident:.1e}; sup_t (1/L)sum_l <E_l^2>_t = {worst_E2:.3f} '
          f'<= {bound_E2:.3f}; sup_t (1/L)sum_l S_N([1,l])(t) = {worst_SN:.3f} <= {bound_SN:.3f} '
          f'(largest single cut {worst_single:.3f})')
print(f'[time {time.time() - T0:.1f}s]')

# ----------------------------------------------------------------------------------------------
# C. random-walk potential (standard single-link Z2 LGT)
# ----------------------------------------------------------------------------------------------
print()
print('=' * 100)
print('C. Standard Z2 LGT (single-link electric term, boundary link): random-walk potential (Remark rem:randomwalk)')
print('=' * 100)
I2 = sp.identity(2, format='csr'); Zm = sp.csr_matrix(np.diag([1., -1.]))
Xm = sp.csr_matrix(np.array([[0., 1.], [1., 0.]])); am = sp.csr_matrix(np.array([[0., 1.], [0., 0.]]))


def op_on(ops, nq):
    out = sp.identity(1, format='csr')
    for k in range(nq):
        out = sp.kron(out, ops.get(k, I2), format='csr')
    return out


def build_std(N, J, h):
    nq = 2 * N
    f = []
    for j in range(N):
        d = {k: Zm for k in range(j)}; d[j] = am; f.append(op_on(d, nq))
    n = [f[j].T @ f[j] for j in range(N)]
    sz = [op_on({N + b: Zm}, nq) for b in range(N)]; sx = [op_on({N + b: Xm}, nq) for b in range(N)]
    H = sp.csr_matrix((2 ** nq, 2 ** nq))
    for b in range(N - 1):
        hop = f[b].T @ f[b + 1]; H = H - J * (sz[b] @ (hop + hop.T))
    for b in range(N):
        H = H - h * sx[b]
    Id = sp.identity(2 ** nq, format='csr')
    G = []
    for j in range(N):
        par = Id - 2 * n[j]
        left = sx[j - 1] if j >= 1 else Id
        G.append(left @ par @ sx[j])
    return H.tocsr(), G, sum(n)


def V_rw(N, h, q):
    Qb = np.cumprod(q)
    return np.array([-h * Qb.sum() + 2 * h * Qb[x:].sum() for x in range(N)])


J, h = 1.0, 0.7
for N in (3, 4, 5):
    H, G, Ntot = build_std(N, J, h)
    d = Ntot.diagonal(); ids = np.where(np.abs(d - 1) < 1e-9)[0]
    Hs = H[ids][:, ids].toarray(); Gs = [g[ids][:, ids].toarray() for g in G]
    X = sum(rng.normal() * g for g in Gs); ev, V = np.linalg.eigh(X)
    lab = np.sign(np.array([np.einsum('ij,ij->j', V, g @ V) for g in Gs]).T).astype(int)
    err = 0; nsec = 0
    for q in itertools.product([1, -1], repeat=N):
        cols = np.where(np.all(lab == np.array(q)[None, :], axis=1))[0]
        B = V[:, cols]; Es = np.sort(np.linalg.eigvalsh(B.T @ Hs @ B))
        A = -J * (np.eye(N, k=1) + np.eye(N, k=-1)) + np.diag(V_rw(N, h, np.array(q)))
        pred = np.sort(np.linalg.eigvalsh(A))
        err = max(err, np.abs(Es - pred).max()); nsec += 1
    print(f'N={N}: one-fermion sector, {nsec} charge sectors, max|dE| (sector by sector) = {err:.2e}')

# two fermions: rank-alternating potential, and failure of any free-fermion description
N = 5
H, G, Ntot = build_std(N, J, h)
d = Ntot.diagonal()
for q in [(1, 1, 1, 1, 1), (1, -1, -1, 1, -1)]:
    specs = {}
    for k in (0, 1, 2):
        ids = np.where(np.abs(d - k) < 1e-9)[0]
        Hs = H[ids][:, ids].toarray(); Gs = [g[ids][:, ids].toarray() for g in G]
        P = np.eye(len(ids))
        for j in range(N):
            P = P @ (np.eye(len(ids)) + q[j] * Gs[j]) / 2
        wv, vv = np.linalg.eigh(P); B = vv[:, wv > 0.5]
        specs[k] = np.sort(np.linalg.eigvalsh(B.T @ Hs @ B))
    Qb = np.cumprod(q); Sx = np.array([Qb[x:].sum() for x in range(N)])
    pairs = [(x, y) for x in range(N) for y in range(x + 1, N)]; pid = {p: i for i, p in enumerate(pairs)}
    H2 = np.zeros((len(pairs), len(pairs)))
    for (x, y), i in pid.items():
        H2[i, i] = -h * Qb.sum() + 2 * h * (Sx[x] - Sx[y])
        for (x2, y2) in [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]:
            if (x2, y2) in pid:
                H2[i, pid[(x2, y2)]] = -J
    e_rank = np.abs(np.sort(np.linalg.eigvalsh(H2)) - specs[2]).max()
    e1, e0 = specs[1], specs[0]
    free = np.sort([e1[a] + e1[b] - e0[0] for a in range(N) for b in range(a + 1, N)])
    print(f'N={N}, sector q={q}: 2-fermion spectrum vs rank-alternating formula: max|dE| = {e_rank:.1e}; '
          f'vs best free-fermion guess e1_a+e1_b-e0: max|dE| = {np.abs(free - specs[2]).max():.3f} (nonzero => not free)')
print(f'[time {time.time() - T0:.1f}s]')

# eigenfunction correlators: random walk vs iid Bernoulli
print('Eigenfunction correlator E Q(x0,x0+d) (N=600, x0=300, 200 samples):')
for h in (0.25, 0.5):
    N = 600; x0 = 300; Sm = 200
    Qrw = np.zeros(N); Qbe = np.zeros(N)
    for s in range(Sm):
        Qb = rng.choice([-1, 1], size=N)
        V = 2 * h * np.cumsum(Qb[::-1])[::-1]; V = V - V[x0]
        E, U = eigh_tridiagonal(V, -J * np.ones(N - 1)); A = np.abs(U); Qrw += A @ A[x0] / Sm
        q = rng.choice([-1, 1], size=N)
        E2, U2 = eigh_tridiagonal(2 * h * q.astype(float), -J * np.ones(N - 1)); A2 = np.abs(U2); Qbe += A2 @ A2[x0] / Sm
    ds = [1, 4, 16, 32, 64, 128]
    print(f'  h/J={h}: d = {ds}')
    print('     random-walk V  : ' + ' '.join(f'{Qrw[x0 + d]:.1e}' for d in ds))
    print('     iid Bernoulli  : ' + ' '.join(f'{Qbe[x0 + d]:.1e}' for d in ds))
print(f'[total time {time.time() - T0:.1f}s]')
