"""
Exact-diagonalization check of Section 2 (sec:model) of the paper, in the paper's exact conventions.

Conventions (BRIEF.md, section 2):
  sites j = 1..N (spinless fermions f_j), links b_j = (j, j+1), j = 1..N, b_N = boundary link (no hopping);
  H = -J sum_{j=1}^{N-1} sz_{b_j} (f_j^+ f_{j+1} + h.c.) - h sum_{j=1}^N sx_{b_{j-1}} sx_{b_j},  sx_{b_0} := 1;
  G_j = sx_{b_{j-1}} (-1)^{n_j} sx_{b_j};
  U(|n> x |q>) = |n> x |x(n,q)>,  x_{b_j}(n,q) = prod_{k<=j} q_k (-1)^{n_k}  (sigma^x eigenbasis);
  U^+ H U = sum_q (dGamma(h_q) - h sum_j q_j) x |q><q|,   h_q = -J Delta_N + 2h diag(q).

Qubit layout.  Physical space H_N: qubits 0..N-1 = fermion modes f_1..f_N (Jordan-Wigner, |n> =
(f_1^+)^{n_1}...(f_N^+)^{n_N}|0>), qubits N..2N-1 = links b_1..b_N in the sigma^z basis (bit 0 = up).
Dual space K_N: qubits 0..N-1 = c-modes (same JW convention), qubits N..2N-1 = tau-spins in the tau^x
eigenbasis (bit 0 <-> q_j = +1).  Hence tau^x_j = diag(1,-1) and tau^z_j = bit flip on that qubit.
U is built explicitly as (1 x Hadamard^{(x)N}) . Permutation; nothing about the claimed dual operators is
used in the construction of H or U.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.linalg import expm, eigh

np.set_printoptions(linewidth=150)
I2 = sp.identity(2, format='csr')
Zm = sp.csr_matrix(np.diag([1., -1.]))
Xm = sp.csr_matrix(np.array([[0., 1.], [1., 0.]]))
Am = sp.csr_matrix(np.array([[0., 1.], [0., 0.]]))   # annihilator: |1> -> |0>
Had = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
LN2 = np.log(2.)


def op_on(ops, nq):
    out = sp.identity(1, format='csr')
    for k in range(nq):
        out = sp.kron(out, ops.get(k, I2), format='csr')
    return out


def ferm_ops(N, nq):
    """JW annihilators on qubits 0..N-1 of an nq-qubit register."""
    f = []
    for j in range(N):
        d = {k: Zm for k in range(j)}
        d[j] = Am
        f.append(op_on(d, nq))
    return f


def maxabs(A):
    A = sp.csr_matrix(A)
    return 0.0 if A.nnz == 0 else float(np.max(np.abs(A.data)))


def bits(idx, nq):
    return [(idx >> (nq - 1 - k)) & 1 for k in range(nq)]


def index(bl):
    v = 0
    for b in bl:
        v = 2 * v + b
    return v


def build_U(N):
    """U: K_N -> H_N.  Dual tau-qubit bit 0 <-> q=+1; link qubit after Hadamard: bit 0 <-> x=+1."""
    nq = 2 * N
    D = 2 ** nq
    rows, cols = [], []
    for col in range(D):
        b = bits(col, nq)
        n, qb = b[:N], b[N:]
        q = [1 - 2 * t for t in qb]
        x, run = [], 1
        for k in range(N):
            run *= q[k] * (-1) ** n[k]
            x.append(run)
        xb = [(1 - t) // 2 for t in x]
        rows.append(index(n + xb))
        cols.append(col)
    P = sp.csr_matrix((np.ones(D), (rows, cols)), shape=(D, D))
    HN = np.array([[1.]])
    for _ in range(N):
        HN = np.kron(HN, Had)
    return sp.kron(sp.identity(2 ** N), sp.csr_matrix(HN), format='csr') @ P


def physical(N, J, h, variant='bl'):
    """variant 'bl': boundary link (paper's model); 'open': N-1 links; 'ring': N links, bond (N,1)."""
    M = N - 1 if variant == 'open' else N
    nq = N + M
    f = ferm_ops(N, nq)
    n = [f[j].T @ f[j] for j in range(N)]
    Id = sp.identity(2 ** nq, format='csr')
    par = [Id - 2 * n[j] for j in range(N)]
    sz = [op_on({N + b: Zm}, nq) for b in range(M)]
    sx = [op_on({N + b: Xm}, nq) for b in range(M)]

    def sxl(b):  # b = link index 0..M-1 (link b_{b+1}); b=-1 is b_0
        if variant == 'ring':
            return sx[b % N]
        if b < 0 or b >= M:
            return Id
        return sx[b]
    H = sp.csr_matrix((2 ** nq, 2 ** nq))
    nb = N if variant == 'ring' else N - 1
    for b in range(nb):
        i, j = b, (b + 1) % N
        hop = f[i].T @ f[j]
        H = H - J * (sz[b] @ (hop + hop.T))
    G = []
    for j in range(N):
        H = H - h * (sxl(j - 1) @ sxl(j))
        G.append(sxl(j - 1) @ par[j] @ sxl(j))
    return dict(H=H.tocsr(), f=f, n=n, par=par, sz=sz, sx=sx, G=G, nq=nq, Id=Id)


def dual(N, J, h):
    nq = 2 * N
    c = ferm_ops(N, nq)
    n = [c[j].T @ c[j] for j in range(N)]
    Id = sp.identity(2 ** nq, format='csr')
    par = [Id - 2 * n[j] for j in range(N)]
    tx = [op_on({N + j: Zm}, nq) for j in range(N)]
    tz = [op_on({N + j: Xm}, nq) for j in range(N)]
    H = sp.csr_matrix((2 ** nq, 2 ** nq))
    for j in range(N - 1):
        hop = c[j].T @ c[j + 1]
        H = H - J * (hop + hop.T)
    for j in range(N):
        H = H - h * (tx[j] @ par[j])
    return dict(H=H.tocsr(), c=c, n=n, par=par, tx=tx, tz=tz, Id=Id, nq=nq)


def hq(N, J, h, q, ring=False, W=1):
    A = np.zeros((N, N))
    for j in range(N - 1):
        A[j, j + 1] = A[j + 1, j] = -J
    if ring:
        A[N - 1, 0] = A[0, N - 1] = -J * W
    return A + np.diag(2 * h * np.asarray(q, float))


def mb_energies(eps, parity=None):
    out = []
    for S in itertools.product([0, 1], repeat=len(eps)):
        if parity is not None and sum(S) % 2 != parity:
            continue
        out.append(float(np.dot(S, eps)))
    return out


# ---------------------------------------------------------------------------------------------
def check_operator_dictionary(N, J, h):
    P = physical(N, J, h)
    Dl = dual(N, J, h)
    U = build_U(N)
    Ud = U.T.conj()
    e = {}
    e['unitarity'] = maxabs(Ud @ U - P['Id'])
    e['H'] = maxabs(Ud @ P['H'] @ U - Dl['H'])
    e['G'] = max(maxabs(Ud @ P['G'][j] @ U - Dl['tx'][j]) for j in range(N))
    e['f'] = max(maxabs(Ud @ P['f'][j] @ U - Dl['c'][j] @ Dl['tz'][j]) for j in range(N))
    e['n'] = max(maxabs(Ud @ P['n'][j] @ U - Dl['n'][j]) for j in range(N))
    err = 0.
    for j in range(N):
        tgt = Dl['tz'][j] @ (Dl['tz'][j + 1] if j + 1 < N else Dl['Id'])
        err = max(err, maxabs(Ud @ P['sz'][j] @ U - tgt))
    e['sz'] = err
    err = 0.
    for j in range(N):
        tgt = Dl['Id']
        for k in range(j + 1):
            tgt = tgt @ Dl['tx'][k] @ Dl['par'][k]
        err = max(err, maxabs(Ud @ P['sx'][j] @ U - tgt))
    e['sx'] = err
    # dressed fermions c~_j = U c_j U^+ = f_j prod_{k>=j} sz_{b_k}; gauge invariance; CAR
    err = errG = 0.
    ct = []
    for j in range(N):
        S = P['Id']
        for k in range(j, N):
            S = S @ P['sz'][k]
        ctj = P['f'][j] @ S
        ct.append(ctj)
        err = max(err, maxabs(U @ Dl['c'][j] @ Ud - ctj))
        for k in range(N):
            errG = max(errG, maxabs(ctj @ P['G'][k] - P['G'][k] @ ctj))
    e['dressed'] = err
    e['dressed_gauge_inv'] = errG
    errC = 0.
    for i in range(N):
        for j in range(N):
            errC = max(errC, maxabs(ct[i] @ ct[j].T + ct[j].T @ ct[i] - (P['Id'] if i == j else 0 * P['Id'])))
    e['dressed_CAR'] = errC
    # parity, 1D Gauss law sx_{b_j} = (-1)^{N_{<=j}} prod_{k<=j} G_k, H in dressed variables
    parN = P['Id']
    for j in range(N):
        parN = parN @ P['par'][j]
    parNd = Dl['Id']
    for j in range(N):
        parNd = parNd @ Dl['par'][j]
    e['parity'] = maxabs(Ud @ parN @ U - parNd)
    err = 0.
    for j in range(N):
        tgt = P['Id']
        for k in range(j + 1):
            tgt = tgt @ P['par'][k] @ P['G'][k]
        err = max(err, maxabs(P['sx'][j] - tgt))
    e['gauss_law_1d'] = err
    Hd = sp.csr_matrix(P['H'].shape)
    for j in range(N - 1):
        hop = ct[j].T @ ct[j + 1]
        Hd = Hd - J * (hop + hop.T)
    for j in range(N):
        Hd = Hd - h * (P['G'][j] @ P['par'][j])
    e['H_dressed_form'] = maxabs(P['H'] - Hd)
    prodG = P['Id']
    for j in range(N):
        prodG = prodG @ P['G'][j]
    e['prodG=(-1)^Nf sx_bN'] = maxabs(prodG - parN @ P['sx'][N - 1])
    return e


def check_spectra(N, J, h):
    out = {}
    # boundary-link model
    P = physical(N, J, h, 'bl')
    Ev = np.sort(np.linalg.eigvalsh(P['H'].toarray()))
    pred = []
    for q in itertools.product([1, -1], repeat=N):
        pred += [e - h * sum(q) for e in mb_energies(np.linalg.eigvalsh(hq(N, J, h, q)))]
    out['bl'] = (len(Ev), float(np.max(np.abs(Ev - np.sort(pred)))))
    # open: parity constraint prod q = (-1)^{N_f}
    P = physical(N, J, h, 'open')
    Ev = np.sort(np.linalg.eigvalsh(P['H'].toarray()))
    pred = []
    for q in itertools.product([1, -1], repeat=N):
        p = 0 if np.prod(q) == 1 else 1
        pred += [e - h * sum(q) for e in mb_energies(np.linalg.eigvalsh(hq(N, J, h, q)), p)]
    out['open'] = (len(Ev), float(np.max(np.abs(Ev - np.sort(pred)))))
    if N >= 3:
        P = physical(N, J, h, 'ring')
        Ev = np.sort(np.linalg.eigvalsh(P['H'].toarray()))
        pred = []
        for q in itertools.product([1, -1], repeat=N):
            p = 0 if np.prod(q) == 1 else 1
            for W in (1, -1):
                pred += [e - h * sum(q) for e in mb_energies(np.linalg.eigvalsh(hq(N, J, h, q, True, W)), p)]
        out['ring'] = (len(Ev), float(np.max(np.abs(Ev - np.sort(pred)))))
    return out


def check_open_ring_operator_identities(N, J, h):
    """open: H_bl restricted to sx_{b_N}=+1 equals H_open; U^+ P^+_{b_N} U = projector on prod q = (-1)^{N_c}.
       ring: U^+ H_ring U = -J sum_{j<N} hop - J tau^z_1 (c_N^+ c_1 + h.c.) - h (-1)^{N_c} prod_{k>=2} tau^x_k (-1)^{n_1}
             - h sum_{j>=2} tau^x_j (-1)^{n_j};  U^+ W U = tau^z_1;  U^+ G^ring_1 U = (-1)^{N_c} prod_{k>=2} tau^x_k."""
    e = {}
    P = physical(N, J, h, 'bl')
    Po = physical(N, J, h, 'open')
    Dl = dual(N, J, h)
    U = build_U(N)
    Ud = U.T.conj()
    plus = np.array([[1.], [1.]]) / np.sqrt(2)
    V = sp.kron(sp.identity(2 ** (2 * N - 1)), sp.csr_matrix(plus), format='csr')  # embed |+> on b_N
    e['open: [H,sx_bN]'] = maxabs(P['H'] @ P['sx'][N - 1] - P['sx'][N - 1] @ P['H'])
    e['open: V^+ H V = H_open'] = maxabs(V.T @ P['H'] @ V - Po['H'])
    parNd = Dl['Id']
    prodtx = Dl['Id']
    for j in range(N):
        parNd = parNd @ Dl['par'][j]
        prodtx = prodtx @ Dl['tx'][j]
    Pp = (P['Id'] + P['sx'][N - 1]) / 2
    e['open: U^+ P+ U'] = maxabs(Ud @ Pp @ U - (Dl['Id'] + prodtx @ parNd) / 2)
    if N >= 3:
        Pr = physical(N, J, h, 'ring')
        T = Dl['Id']
        for k in range(1, N):
            T = T @ Dl['tx'][k]
        Hr = sp.csr_matrix(Dl['H'].shape)
        for j in range(N - 1):
            hop = Dl['c'][j].T @ Dl['c'][j + 1]
            Hr = Hr - J * (hop + hop.T)
        hop = Dl['c'][N - 1].T @ Dl['c'][0]
        Hr = Hr - J * (Dl['tz'][0] @ (hop + hop.T))
        Hr = Hr - h * (parNd @ T @ Dl['par'][0])
        for j in range(1, N):
            Hr = Hr - h * (Dl['tx'][j] @ Dl['par'][j])
        e['ring: U^+ H U'] = maxabs(Ud @ Pr['H'] @ U - Hr)
        Wop = Pr['Id']
        for b in range(N):
            Wop = Wop @ Pr['sz'][b]
        e['ring: U^+ W U = tau^z_1'] = maxabs(Ud @ Wop @ U - Dl['tz'][0])
        e['ring: U^+ G_1 U'] = maxabs(Ud @ Pr['G'][0] @ U - parNd @ T)
        e['ring: U^+ G_j U (j>=2)'] = max(maxabs(Ud @ Pr['G'][j] @ U - Dl['tx'][j]) for j in range(1, N))
    return e


def check_ring_sectorwise(N, J, h):
    """Project the ring Hamiltonian onto joint eigenspaces of (G_j, W); compare with twist W and with the
    (wrong) JW-type twist W(-1)^{N_f+1}."""
    Pr = physical(N, J, h, 'ring')
    D = 2 ** Pr['nq']
    H = Pr['H'].toarray()
    Wop = np.eye(D)
    for b in range(N):
        Wop = Wop @ Pr['sz'][b].toarray()
    Gs = [g.toarray() for g in Pr['G']]
    errA = errB = 0.
    tot = 0
    for q in itertools.product([1, -1], repeat=N):
        Pq = np.eye(D)
        for j in range(N):
            Pq = Pq @ (np.eye(D) + q[j] * Gs[j]) / 2
        for W in (1, -1):
            PW = Pq @ (np.eye(D) + W * Wop) / 2
            w, v = np.linalg.eigh(PW)
            B = v[:, w > 0.5]
            if B.shape[1] == 0:
                continue
            tot += B.shape[1]
            Es = np.sort(np.linalg.eigvalsh(B.T @ H @ B))
            p = 0 if np.prod(q) == 1 else 1
            pa = np.sort([e - h * sum(q) for e in mb_energies(np.linalg.eigvalsh(hq(N, J, h, q, True, W)), p)])
            pb = np.sort([e - h * sum(q) for e in mb_energies(np.linalg.eigvalsh(hq(N, J, h, q, True, W * (-1) ** (p + 1))), p)])
            errA = max(errA, np.abs(Es - pa).max() if len(pa) == len(Es) else np.inf)
            errB = max(errB, np.abs(Es - pb).max() if len(pb) == len(Es) else np.inf)
    return tot, D, errA, errB


# ---------------------------------------------------------------------------------------------
def split_perm(N, l):
    """permutation matrix Pi with (Pi v)[split index] = v[standard index]; standard layout
    [f_1..f_N, b_1..b_N], split layout [f_1..f_l, b_1..b_l, f_{l+1}..f_N, b_{l+1}..b_N]."""
    nq = 2 * N
    order = list(range(l)) + list(range(N, N + l)) + list(range(l, N)) + list(range(N + l, 2 * N))
    D = 2 ** nq
    rows, cols = [], []
    for s in range(D):
        b = bits(s, nq)
        rows.append(index([b[k] for k in order]))
        cols.append(s)
    return sp.csr_matrix((np.ones(D), (rows, cols)), shape=(D, D))


def build_U_right(N, l):
    """U_{>l} on 2(N-l) qubits [f_{l+1}..f_N, b_{l+1}..b_N]: x_{b_j} = prod_{l<k<=j} q_k (-1)^{n_k}."""
    return build_U(N - l)


def check_cut_locality(N, J, h):
    U = build_U(N)
    P = physical(N, J, h)
    err = 0.
    for l in range(1, N):
        Pi = split_perm(N, l)
        Ul = build_U(l)
        Ur = build_U_right(N, l)
        Uprod = Pi.T @ sp.kron(Ul, Ur, format='csr') @ Pi
        V = P['Id']
        for k in range(l, N):
            V = V @ P['sz'][k]
        Pp = (P['Id'] + P['sx'][l - 1]) / 2
        Pm = (P['Id'] - P['sx'][l - 1]) / 2
        C = Pp + Pm @ V
        err = max(err, maxabs(U - C @ Uprod))
    return err


def entropy_A(vec, N, l):
    """von Neumann entropy of qubits {0..l-1} u {N..N+l-1} (left segment: modes 1..l and links/tau 1..l)."""
    nq = 2 * N
    psi = vec.reshape([2] * nq)
    A = list(range(l)) + list(range(N, N + l))
    B = [k for k in range(nq) if k not in A]
    M = np.transpose(psi, A + B).reshape(2 ** len(A), -1)
    s = np.linalg.svd(M, compute_uv=False)
    p = s ** 2
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log(p)))


def check_entropy_defect(N, rng, ntrial=200):
    U = build_U(N).toarray()
    D = 2 ** (2 * N)
    worst = 0.
    for _ in range(ntrial):
        phi = rng.normal(size=D) + 1j * rng.normal(size=D)
        phi /= np.linalg.norm(phi)
        for l in range(1, N):
            worst = max(worst, abs(entropy_A(U @ phi, N, l) - entropy_A(phi, N, l)))
    # definite charge, definite parity: equality; also fixed (parity, prod_{k<=l} q_k) superpositions
    worst_eq = 0.
    worst_eq2 = 0.
    for _ in range(ntrial // 4):
        for pi in (0, 1):
            phi_m = rng.normal(size=2 ** N) + 1j * rng.normal(size=2 ** N)
            for s in range(2 ** N):
                if sum(bits(s, N)) % 2 != pi:
                    phi_m[s] = 0
            phi_m /= np.linalg.norm(phi_m)
            q = rng.integers(0, 2, size=N)
            vq = np.zeros(2 ** N)
            vq[index(list(q))] = 1
            phi = np.kron(phi_m, vq)
            for l in range(1, N):
                worst_eq = max(worst_eq, abs(entropy_A(U @ phi, N, l) - entropy_A(phi, N, l)))
            # superposition over q with fixed Q_l = prod_{k<=l} q_k, fixed parity
            for l in range(1, N):
                Ql = rng.choice([1, -1])
                phi2 = np.zeros(D, complex)
                for s in range(D):
                    b = bits(s, 2 * N)
                    n, qb = b[:N], b[N:]
                    if sum(n) % 2 != pi:
                        continue
                    if np.prod([1 - 2 * t for t in qb[:l]]) != Ql:
                        continue
                    phi2[s] = rng.normal() + 1j * rng.normal()
                phi2 /= np.linalg.norm(phi2)
                worst_eq2 = max(worst_eq2, abs(entropy_A(U @ phi2, N, l) - entropy_A(phi2, N, l)))
    return worst, worst_eq, worst_eq2


def sharpness_examples():
    """N=2, l=1.  (i) Phi = |00> x (|q1=+>+|q1=->)/sqrt2 x |q2=+>: S_A'(Phi)=0, S_A(U Phi)=ln2.
                  (ii) Psi = |00> x |up>_{b1} x |+>_{b2}: S_A(Psi)=0, S_A'(U^+ Psi)=ln2."""
    N = 2
    U = build_U(N).toarray()
    e0 = np.zeros(4)
    e0[0] = 1
    pl = np.array([1, 1.]) / np.sqrt(2)
    phi = np.kron(e0, np.kron(pl, np.array([1, 0.])))
    a = (entropy_A(phi, N, 1), entropy_A(U @ phi, N, 1))
    psi = np.kron(e0, np.kron(np.array([1, 0.]), pl))
    b = (entropy_A(psi, N, 1), entropy_A(U.T @ psi, N, 1))
    return a, b


# ---------------------------------------------------------------------------------------------
def dGamma(N, A):
    """second quantization of a one-body N x N matrix on the 2^N-dim JW Fock space."""
    c = ferm_ops(N, N)
    out = sp.csr_matrix((2 ** N, 2 ** N), dtype=complex)
    for x in range(N):
        for y in range(N):
            if A[x, y] != 0:
                out = out + A[x, y] * (c[x].T @ c[y])
    return out.toarray()


def check_dynamics(N, J, h, times, occ, rng, general_psi=False):
    """Psi0 = |up..up>_z x psi.  Compare ED of the full LGT with the sector-average formulas of Cor. cor:sector:
    n_x, dressed bilinears f_x^+ (prod_{x<=k<y} sz) f_y, sx strings, N_A distribution, and the Loschmidt
    formula for sz_{b_j} (with the phase e^{iht(sum q - sum q')})."""
    P = physical(N, J, h)
    D = 2 ** P['nq']
    if general_psi:
        vf = rng.normal(size=2 ** N) + 1j * rng.normal(size=2 ** N)
        vf /= np.linalg.norm(vf)
    else:
        vf = np.zeros(2 ** N, complex)
        vf[index([1 if j in occ else 0 for j in range(N)])] = 1
    vs = np.zeros(2 ** N)
    vs[0] = 1
    psi0 = np.kron(vf, vs)
    Hd = P['H'].toarray()
    w, V = np.linalg.eigh(Hd)
    c0 = V.conj().T @ psi0
    # dressed bilinear operators and N_A projectors
    bil = {}
    for x in range(N):
        for y in range(N):
            if x < y:
                S = P['Id']
                for k in range(x, y):
                    S = S @ P['sz'][k]
                bil[(x, y)] = (P['f'][x].T @ S @ P['f'][y]).tocsr()
    l = N // 2
    NA = sum(P['n'][x] for x in range(l)).diagonal()
    sxd = [P['sx'][b] for b in range(N)]
    szd = [P['sz'][b] for b in range(N)]
    nd = [P['n'][x].diagonal() for x in range(N)]
    # dual side, Fock space operators
    cF = ferm_ops(N, N)
    nF = [(cF[x].T @ cF[x]).diagonal() for x in range(N)]
    NAF = sum(nF[x] for x in range(l))
    errs = []
    qs = list(itertools.product([1, -1], repeat=N))
    for t in times:
        psi = V @ (np.exp(-1j * w * t) * c0)
        ED = []
        ED += [np.vdot(psi, nd[x] * psi).real for x in range(N)]
        ED += [np.vdot(psi, bil[k] @ psi) for k in sorted(bil)]
        ED += [np.vdot(psi, sxd[b] @ psi).real for b in range(N)]
        ED += [np.vdot(psi, szd[b] @ psi).real for b in range(N)]
        ED += [np.sum(np.abs(psi[NA == k]) ** 2) for k in range(l + 1)]
        ED = np.array(ED)
        # sector formulas
        psiq = {}
        for q in qs:
            psiq[q] = expm(-1j * t * dGamma(N, hq(N, J, h, q))) @ vf
        PR = np.zeros(len(ED), complex)
        for q in qs:
            ph = psiq[q]
            row = []
            row += [np.vdot(ph, nF[x] * ph).real for x in range(N)]
            for (x, y) in sorted(bil):
                row.append(np.vdot(ph, (cF[x].T @ cF[y]) @ ph))
            for b in range(N):
                Qb = np.prod(q[:b + 1])
                strg = np.ones(2 ** N)
                for k in range(b + 1):
                    strg = strg * (1 - 2 * nF[k])
                row.append(Qb * np.vdot(ph, strg * ph).real)
            for b in range(N):
                qq = list(q)
                qq[b] *= -1
                if b + 1 < N:
                    qq[b + 1] *= -1
                qq = tuple(qq)
                phase = np.exp(1j * h * t * (sum(q) - sum(qq)))
                row.append(phase * np.vdot(psiq[qq], ph))
            row += [np.sum(np.abs(ph[NAF == k]) ** 2) for k in range(l + 1)]
            PR += np.array(row) / 2 ** N
        errs.append(float(np.max(np.abs(ED - PR))))
    return max(errs)


def check_slater_formulas(N, J, h, t, occ):
    """Fock psi: <psi_q'(t), psi_q(t)> = det(Phi_q'^* Phi_q), <(-1)^{N_<=j}> = det(Phi^* R_j Phi)."""
    vf = np.zeros(2 ** N, complex)
    vf[index([1 if j in occ else 0 for j in range(N)])] = 1
    cF = ferm_ops(N, N)
    nF = [(cF[x].T @ cF[x]).diagonal() for x in range(N)]
    occl = sorted(occ)
    err = 0.
    for q in itertools.product([1, -1], repeat=N):
        qq = list(q)
        qq[0] *= -1
        qq[1] *= -1
        a = expm(-1j * t * dGamma(N, hq(N, J, h, q))) @ vf
        b = expm(-1j * t * dGamma(N, hq(N, J, h, qq))) @ vf
        Ua = expm(-1j * t * hq(N, J, h, q))[:, occl]
        Ub = expm(-1j * t * hq(N, J, h, qq))[:, occl]
        err = max(err, abs(np.vdot(b, a) - np.linalg.det(Ub.conj().T @ Ua)))
        for j in range(N):
            strg = np.ones(2 ** N)
            for k in range(j + 1):
                strg = strg * (1 - 2 * nF[k])
            R = np.diag([-1. if k <= j else 1. for k in range(N)])
            err = max(err, abs(np.vdot(a, strg * a) - np.linalg.det(Ua.conj().T @ R @ Ua)))
    return err


# ---------------------------------------------------------------------------------------------
def gauge_product_state(N, kinds):
    """kinds[j] in {'zu','zd','x+','x-'} or a 2-vector (amplitudes in the sigma^z basis)."""
    vecs = {'zu': np.array([1, 0.]), 'zd': np.array([0, 1.]),
            'x+': np.array([1, 1.]) / np.sqrt(2), 'x-': np.array([1, -1.]) / np.sqrt(2)}
    g = np.array([1.])
    for k in kinds:
        g = np.kron(g, vecs[k] if isinstance(k, str) else np.asarray(k))
    return g


def x_of(n, q):
    x, run = [], 1
    for k in range(len(n)):
        run *= q[k] * (-1) ** n[k]
        x.append(run)
    return x


def check_product_laws(N, rng):
    U = build_U(N).toarray()
    out = {}
    # (c) z-product with down set D, general psi: U^+ Psi0 = 2^{-N/2} sum_q eps_D(q) (W_D psi) x |q>
    s = rng.choice(['zu', 'zd'], size=N)
    Dset = [j for j in range(N) if s[j] == 'zd']
    psi = rng.normal(size=2 ** N) + 1j * rng.normal(size=2 ** N)
    psi /= np.linalg.norm(psi)
    Phi = U.T @ np.kron(psi, gauge_product_state(N, s))
    WD = np.ones(2 ** N)
    for nidx in range(2 ** N):
        nb = bits(nidx, N)
        for j in Dset:
            WD[nidx] *= (-1) ** sum(nb[:j + 1])
    pred = np.zeros_like(Phi)
    for q in itertools.product([1, -1], repeat=N):
        eps = np.prod([np.prod(q[:j + 1]) for j in Dset]) if Dset else 1
        vq = np.zeros(2 ** N)
        vq[index([(1 - t) // 2 for t in q])] = 1
        pred += 2 ** (-N / 2) * eps * np.kron(WD * psi, vq)
    out['z-product (general psi)'] = float(np.max(np.abs(Phi - pred)))
    # (b) general product gauge state, Fock matter: mu(q) = |chi(x(n,q))|^2, matter component = chi |n>
    chis = []
    for _ in range(N):
        v = rng.normal(size=2) + 1j * rng.normal(size=2)
        chis.append(v / np.linalg.norm(v))
    n = list(rng.integers(0, 2, size=N))
    vf = np.zeros(2 ** N)
    vf[index(n)] = 1
    Phi = U.T @ np.kron(vf, gauge_product_state(N, chis))
    pred = np.zeros_like(Phi)
    plus, minus = np.array([1, 1.]) / np.sqrt(2), np.array([1, -1.]) / np.sqrt(2)
    for q in itertools.product([1, -1], repeat=N):
        x = x_of(n, q)
        amp = np.prod([np.vdot(plus if x[j] == 1 else minus, chis[j]) for j in range(N)])
        vq = np.zeros(2 ** N)
        vq[index([(1 - t) // 2 for t in q])] = 1
        pred += amp * np.kron(vf, vq)
    out['general product, Fock'] = float(np.max(np.abs(Phi - pred)))
    # (d) sigma^x product, Fock: single sector q_j = (-1)^{n_j} xi_{j-1} xi_j
    xi = list(rng.choice([1, -1], size=N))
    Phi = U.T @ np.kron(vf, gauge_product_state(N, ['x+' if v == 1 else 'x-' for v in xi]))
    qd = [(-1) ** n[j] * (xi[j - 1] if j > 0 else 1) * xi[j] for j in range(N)]
    vq = np.zeros(2 ** N)
    vq[index([(1 - t) // 2 for t in qd])] = 1
    out['sigma^x product, Fock'] = float(np.max(np.abs(Phi - np.kron(vf, vq))))
    # (e) dimer state (N even): odd links |up>_z, even links b_{2k} sigma^x eigenvalue (-1)^k, CDW on even sites
    if N % 2 == 0:
        Mh = N // 2
        kinds = []
        for j in range(1, N + 1):
            if j % 2 == 1:
                kinds.append('zu')
            else:
                kinds.append('x+' if (-1) ** (j // 2) == 1 else 'x-')
        ncdw = [1 if (j + 1) % 2 == 0 else 0 for j in range(N)]  # 1-indexed even sites occupied
        vf2 = np.zeros(2 ** N)
        vf2[index(ncdw)] = 1
        Phi = U.T @ np.kron(vf2, gauge_product_state(N, kinds))
        pred = np.zeros_like(Phi)
        for eta in itertools.product([1, -1], repeat=Mh):
            q = [eta[j // 2] for j in range(N)]
            vq = np.zeros(2 ** N)
            vq[index([(1 - t) // 2 for t in q])] = 1
            pred += 2 ** (-Mh / 2) * np.kron(vf2, vq)
        out['dimer law'] = float(np.max(np.abs(Phi - pred)))
    return out


def check_gauss_projection():
    """q = (1,...,1): a_alpha = <phi_alpha, Gamma phi_alpha> = 0 for N even (so the infinite-time CDW imbalance
    vanishes); ballistic second moment <X^2>_t = 2 J^2 t^2 on a long chain."""
    J, h = 1.0, 0.7
    worst = 0.
    for N in range(2, 201, 2):
        E, V = np.linalg.eigh(hq(N, J, h, np.ones(N)))
        Gam = np.array([(-1.) ** (x + 1) for x in range(N)])
        worst = max(worst, np.max(np.abs((V ** 2 * Gam[:, None]).sum(0))))
    N = 1201
    x0 = N // 2
    E, V = np.linalg.eigh(hq(N, J, h, np.ones(N)))
    rel = 0.
    for t in [10., 50., 100.]:
        u = V @ (np.exp(-1j * E * t) * V[x0])
        m2 = ((np.arange(N) - x0) ** 2 * np.abs(u) ** 2).sum()
        rel = max(rel, abs(m2 / (2 * J ** 2 * t ** 2) - 1))
    return worst, rel


def check_clean_correlator():
    """q=(1..1): eigenfunction correlator Q^{[1,N]}(1,N) = sum_alpha |phi_alpha(1) phi_alpha(N)| = 1 for all N."""
    worst = 0.
    for N in range(2, 301):
        E, V = np.linalg.eigh(hq(N, 1.0, 0.7, np.ones(N)))
        worst = max(worst, abs(np.sum(np.abs(V[0] * V[N - 1])) - 1))
    return worst


def dynamical_defect(N, J, h, l, times):
    """|S_A(Psi(t)) - S_A'(U^+ Psi(t))| along the dynamics from |up..up>_z x CDW."""
    P = physical(N, J, h)
    U = build_U(N).toarray()
    vf = np.zeros(2 ** N)
    vf[index([1 if (j + 1) % 2 == 0 else 0 for j in range(N)])] = 1
    vs = np.zeros(2 ** N)
    vs[0] = 1
    psi0 = np.kron(vf, vs)
    w, V = np.linalg.eigh(P['H'].toarray())
    c0 = V.T @ psi0
    out = []
    for t in times:
        psi = V @ (np.exp(-1j * w * t) * c0)
        out.append((t, entropy_A(psi, N, l), entropy_A(U.T @ psi, N, l)))
    return out


def check_loschmidt_nophase(N, J, h, times):
    """<sz_b(t)>: ED vs. eq. (loschmidt) WITHOUT the phase e^{iht(sum q - sum q')} (should fail at O(1))."""
    P = physical(N, J, h)
    vf = np.zeros(2 ** N, complex)
    vf[index([1 if (j + 1) % 2 == 0 else 0 for j in range(N)])] = 1
    vs = np.zeros(2 ** N)
    vs[0] = 1
    psi0 = np.kron(vf, vs)
    w, V = np.linalg.eigh(P['H'].toarray())
    c0 = V.conj().T @ psi0
    qs = list(itertools.product([1, -1], repeat=N))
    err = 0.
    for t in times:
        psi = V @ (np.exp(-1j * w * t) * c0)
        psiq = {q: expm(-1j * t * dGamma(N, hq(N, J, h, q))) @ vf for q in qs}
        for b in range(N):
            ed = np.vdot(psi, P['sz'][b] @ psi).real
            pr = 0.
            for q in qs:
                qq = list(q)
                qq[b] *= -1
                if b + 1 < N:
                    qq[b + 1] *= -1
                pr += np.vdot(psiq[tuple(qq)], psiq[q]).real / 2 ** N
            err = max(err, abs(ed - pr))
    return err


if __name__ == '__main__':
    J, h = 1.0, 0.7
    rng = np.random.default_rng(20261122)
    print('=== (1) operator dictionary of Thm thm:duality (explicit U), J=1, h=0.7 ===')
    for N in [2, 3, 4, 5, 6]:
        e = check_operator_dictionary(N, J, h)
        print(f'N={N}: ' + ', '.join(f'{k}={v:.1e}' for k, v in e.items()))
    print('=== (2) spectra: full LGT vs union of free-fermion sectors ===')
    for N in [2, 3, 4, 5, 6]:
        o = check_spectra(N, J, h)
        print(f'N={N}: ' + ', '.join(f'{k}: dim={v[0]}, max|dE|={v[1]:.1e}' for k, v in o.items()))
    print('=== (3) open / ring variants: operator identities ===')
    for N in [3, 4, 5, 6]:
        e = check_open_ring_operator_identities(N, J, h)
        print(f'N={N}: ' + ', '.join(f'{k}={v:.1e}' for k, v in e.items()))
    for N in [3, 4]:
        tot, D, eA, eB = check_ring_sectorwise(N, J, h)
        print(f'ring N={N}, sector by sector ({tot}/{D}): twist=W err={eA:.1e}; twist=W(-1)^(Nf+1) err={eB:.2f}')
    print('=== (4) cut locality U = C_l (U_<=l x U_>l), all l ===')
    for N in [2, 3, 4, 5, 6]:
        print(f'N={N}: max error = {check_cut_locality(N, J, h):.1e}')
    print('=== (5) entanglement of a left segment under the duality ===')
    for N in [3, 4]:
        w, weq, weq2 = check_entropy_defect(N, rng)
        print(f'N={N}: random states max|S_A(U Phi)-S_A\'(Phi)| = {w:.4f} (ln2 = {LN2:.4f}); '
              f'definite charge+parity: {weq:.1e}; fixed (parity, prod_(k<=l) q_k): {weq2:.1e}')
    a, b = sharpness_examples()
    print(f'sharpness (N=2,l=1): (S_A\'(Phi), S_A(U Phi)) = ({a[0]:.4f}, {a[1]:.4f});  '
          f'(S_A(Psi), S_A\'(U^+ Psi)) = ({b[0]:.4f}, {b[1]:.4f})')
    print('=== (6) dynamics: ED of the LGT vs the sector-average formulas of Cor. cor:sector ===')
    times = np.linspace(0, 10.0, 11)
    for N in [4, 5, 6]:
        for hh in [0.4, 1.3]:
            occ = set(range(1, N, 2))       # 1-indexed even sites (CDW)
            e1 = check_dynamics(N, J, hh, times, occ, rng)
            e2 = check_dynamics(N, J, hh, times[:6], occ, rng, general_psi=True)
            print(f'N={N}, h={hh}: CDW max err = {e1:.1e};  random (non-Fock) psi max err = {e2:.1e}')
    print(f'control: Loschmidt formula WITHOUT the sector phase, N=4, h=1.3: max err = {check_loschmidt_nophase(4, J, 1.3, times):.2f}')
    print(f'Slater/Loschmidt determinant formulas (N=5, t=2.3): err = {check_slater_formulas(5, J, h, 2.3, {0, 2, 4}):.1e}')
    print('=== (7) sector laws of product initial states (Lemma lem:productlaw) ===')
    for N in [2, 3, 4, 5, 6]:
        o = check_product_laws(N, rng)
        print(f'N={N}: ' + ', '.join(f'{k}: {v:.1e}' for k, v in o.items()))
    print('=== (8) Gauss-law projection (Remark rem:gaussproj) ===')
    worst, rel = check_gauss_projection()
    print(f'q=(1..1), N even 2..200: max_alpha |<phi_a, Gamma phi_a>| = {worst:.1e};  '
          f'|<X^2>_t/(2J^2t^2) - 1| (N=1201, t<=100) = {rel:.1e}')
    print(f'q=(1..1): max_N |Q^[1,N](1,N) - 1| (N=2..300) = {check_clean_correlator():.1e}')
    print('=== (9) duality defect along the dynamics, |up..up>_z x CDW, N=6, l=3, h=1 ===')
    for t, sA, sAp in dynamical_defect(6, 1.0, 1.0, 3, [0., 0.5, 1., 2., 5., 10., 50.]):
        print(f'   t={t:5.1f}: S_A(Psi)={sA:.4f}, S_A\'(U^+Psi)={sAp:.4f}, |diff|={abs(sA - sAp):.4f}')
