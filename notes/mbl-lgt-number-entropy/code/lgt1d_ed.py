"""
Exact diagonalization check of the 1D Z2 LGT (SKKM-type) -> free fermions with binary potential.

Full Hilbert space: N fermion modes (Jordan-Wigner qubits 0..N-1) + M link spins (qubits N..N+M-1).
Spin basis: standard sigma^z basis (|up>=(1,0)).  Nothing about the mapping is built in.

Variants:
  'open_bl' : open chain, links (j,j+1) for j=1..N, the last one a static boundary link (N,N+1);
              electric term  -h sum_{j=1}^N sx_{j-1,j} sx_{j,j+1},  sx_{0,1}:=1.
              Claim: H = (+)_q [ dGamma(h_q) - h sum q ],  h_q = -J Lap_open + 2h diag(q), q in {+-1}^N free.
  'open'    : open chain, links (j,j+1) j=1..N-1 only; electric -h sx_{12} - h sum_{j=2}^{N-1} sx sx - h sx_{N-1,N}.
              Claim: same, but with parity constraint prod q = (-1)^{N_f}.
  'ring'    : ring of N sites, N links; claim: sectors (q,W), prod q=(-1)^{N_f}, h_{q,W} = ring hopping with
              factor W on bond (N,1).
  'std_bl'  : standard single-link electric term -h sum_b sx_b (with boundary link), one-particle sector:
              claim: h = -J Lap + V_q(x),  V_q(x) = -h sum_b Q_b + 2h sum_{b>=x} Q_b,  Q_b = prod_{k<=b} q_k.
"""
import numpy as np, itertools, sys
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import eigh, expm

I2 = sp.identity(2, format='csr')
Zm = sp.csr_matrix(np.diag([1., -1.]))
Xm = sp.csr_matrix(np.array([[0., 1.], [1., 0.]]))
am = sp.csr_matrix(np.array([[0., 1.], [0., 0.]]))   # |0>=empty,|1>=occ ; a|1>=|0>
nm = sp.csr_matrix(np.diag([0., 1.]))

def op_on(ops, nq):
    """ops: dict qubit->2x2 sparse; returns kron product"""
    out = sp.identity(1, format='csr')
    for k in range(nq):
        out = sp.kron(out, ops.get(k, I2), format='csr')
    return out

def build(N, J, h, variant):
    if variant in ('open_bl', 'std_bl'):
        M = N
    elif variant == 'open':
        M = N - 1
    elif variant == 'ring':
        M = N
    nq = N + M
    f = []
    for j in range(N):
        d = {k: Zm for k in range(j)}; d[j] = am
        f.append(op_on(d, nq))
    n = [f[j].T @ f[j] for j in range(N)]
    sz = [op_on({N + b: Zm}, nq) for b in range(M)]
    sx = [op_on({N + b: Xm}, nq) for b in range(M)]
    Id = sp.identity(2 ** nq, format='csr')
    H = sp.csr_matrix((2 ** nq, 2 ** nq))
    # hopping bonds: link b connects sites b and b+1 (0-indexed), ring: link N-1 connects N-1 and 0
    nb = N - 1 if variant != 'ring' else N
    for b in range(nb):
        i, j = b, (b + 1) % N
        hop = f[i].T @ f[j]
        H = H - J * (sz[b] @ (hop + hop.T))
    # electric
    def sxl(b):  # link left of site j is b=j-1, right is b=j
        if variant == 'ring':
            return sx[b % N]
        if b < 0 or b >= M:
            return Id
        return sx[b]
    if variant == 'std_bl':
        for b in range(M):
            H = H - h * sx[b]
    else:
        for j in range(N):
            H = H - h * (sxl(j - 1) @ sxl(j))
    Q = []
    for j in range(N):
        par = Id - 2 * n[j]
        Q.append(sxl(j - 1) @ par @ sxl(j))
    return dict(H=H.tocsr(), f=f, n=n, sz=sz, sx=sx, Q=Q, nq=nq, M=M)

def hq(N, J, h, q, ring=False, W=1):
    A = np.zeros((N, N))
    for j in range(N - 1):
        A[j, j + 1] = A[j + 1, j] = -J
    if ring:
        A[N - 1, 0] = A[0, N - 1] = -J * W
    A += np.diag(2 * h * np.asarray(q, float))
    return A

def mb_energies(eps, parity=None):
    N = len(eps); out = []
    for S in itertools.product([0, 1], repeat=N):
        if parity is not None and sum(S) % 2 != parity:
            continue
        out.append(float(np.dot(S, eps)))
    return out

def check_spectrum(N, J, h, variant):
    m = build(N, J, h, variant)
    Ev = np.linalg.eigvalsh(m['H'].toarray())
    pred = []
    for q in itertools.product([1, -1], repeat=N):
        const = -h * sum(q)
        if variant == 'open_bl':
            eps = np.linalg.eigvalsh(hq(N, J, h, q))
            pred += [e + const for e in mb_energies(eps)]
        elif variant == 'open':
            p = 0 if np.prod(q) == 1 else 1
            eps = np.linalg.eigvalsh(hq(N, J, h, q))
            pred += [e + const for e in mb_energies(eps, p)]
        elif variant == 'ring':
            p = 0 if np.prod(q) == 1 else 1
            for W in (1, -1):
                eps = np.linalg.eigvalsh(hq(N, J, h, q, ring=True, W=W))
                pred += [e + const for e in mb_energies(eps, p)]
    pred = np.sort(pred)
    return len(Ev), len(pred), np.max(np.abs(np.sort(Ev) - pred)) if len(pred) == len(Ev) else None

def check_ring_alt(N, J, h):
    """alternative hypothesis: twist W*(-1)^{N_f+1} (JW-like) -- should FAIL"""
    m = build(N, J, h, 'ring')
    Ev = np.sort(np.linalg.eigvalsh(m['H'].toarray()))
    pred = []
    for q in itertools.product([1, -1], repeat=N):
        const = -h * sum(q); p = 0 if np.prod(q) == 1 else 1
        for W in (1, -1):
            Wt = W * (-1) ** (p + 1)
            eps = np.linalg.eigvalsh(hq(N, J, h, q, ring=True, W=Wt))
            pred += [e + const for e in mb_energies(eps, p)]
    pred = np.sort(pred)
    return np.max(np.abs(Ev - pred))

def slater_overlap(Phi1, Phi2):
    return np.linalg.det(Phi1.conj().T @ Phi2)

def check_dynamics(N, J, h, times, occ):
    """open_bl. Psi0 = |up..up>_z (x) Fock(occ). Compare <n_x(t)>, <sx_b(t)>, <sz_b(t)>, dressed current."""
    m = build(N, J, h, 'open_bl'); nq = m['nq']
    # initial state
    vf = np.zeros(2 ** N); idx = int(''.join('1' if j in occ else '0' for j in range(N)), 2); vf[idx] = 1
    vs = np.zeros(2 ** N); vs[0] = 1     # all up
    psi0 = np.kron(vf, vs).astype(complex)
    psis = expm_multiply(-1j * m['H'], psi0, start=times[0], stop=times[-1], num=len(times), endpoint=True)
    occl = sorted(occ)
    errs = []
    rows = []
    for it, t in enumerate(times):
        psi = psis[it]
        ED_n = np.array([np.vdot(psi, m['n'][x] @ psi).real for x in range(N)])
        ED_sx = np.array([np.vdot(psi, m['sx'][b] @ psi).real for b in range(N)])
        ED_sz = np.array([np.vdot(psi, m['sz'][b] @ psi).real for b in range(N)])
        # dressed current  j_x = <i (f_x^dag sz f_{x+1} - h.c.)>
        ED_cur = []
        for x in range(N - 1):
            o = m['f'][x].T @ m['sz'][x] @ m['f'][x + 1]
            ED_cur.append((1j * (np.vdot(psi, o @ psi) - np.vdot(psi, o.T @ psi))).real)
        ED_cur = np.array(ED_cur)
        # mapped
        Pn = np.zeros(N); Psx = np.zeros(N); Psz = np.zeros(N); Pcur = np.zeros(N - 1)
        Phis = {}
        for q in itertools.product([1, -1], repeat=N):
            U = expm(-1j * t * hq(N, J, h, q))
            Phis[q] = U[:, occl]
        for q in itertools.product([1, -1], repeat=N):
            Phi = Phis[q]; C = Phi @ Phi.conj().T   # C_xy = <c_y^dag c_x>
            Pn += np.real(np.diag(C)) / 2 ** N
            for x in range(N - 1):
                # <c_x^dag c_{x+1}> = C_{x+1,x}
                v = C[x + 1, x]
                Pcur[x] += (1j * (v - np.conj(v))).real / 2 ** N
            for b in range(N):
                Wd = np.diag([(-1.) if k <= b else 1. for k in range(N)])
                sgn = np.prod(q[:b + 1])
                Psx[b] += sgn * np.real(np.linalg.det(Phi.conj().T @ Wd @ Phi)) / 2 ** N
                qq = list(q); qq[b] *= -1
                if b + 1 < N: qq[b + 1] *= -1
                ph = np.exp(1j * t * h * (sum(q) - sum(qq)))  # constants -h*sum(q) differ between sectors
                Psz[b] += np.real(ph * slater_overlap(Phis[tuple(qq)], Phi)) / 2 ** N
        e = max(np.abs(ED_n - Pn).max(), np.abs(ED_sx - Psx).max(), np.abs(ED_sz - Psz).max(), np.abs(ED_cur - Pcur).max())
        errs.append(e)
        rows.append((t, ED_n.round(4), ED_sx.round(4), ED_sz.round(4)))
    return errs, rows

def check_std_oneparticle(N, J, h):
    m = build(N, J, h, 'std_bl')
    H = m['H']
    Ntot = sum(m['n'])
    # restrict to N_f=1
    d = Ntot.diagonal()
    ids = np.where(np.abs(d - 1) < 1e-9)[0]
    Hs = H[ids][:, ids].toarray()
    Ev = np.sort(np.linalg.eigvalsh(Hs))
    pred = []
    for q in itertools.product([1, -1], repeat=N):
        Qb = np.cumprod(q)
        V = np.array([-h * Qb.sum() + 2 * h * Qb[x:].sum() for x in range(N)])
        A = hq(N, J, 0., q) + np.diag(V)
        pred += list(np.linalg.eigvalsh(A))
    pred = np.sort(pred)
    return len(Ev), len(pred), np.max(np.abs(Ev - pred))

def is_free_spectrum_std(N, J, h, q):
    """standard single-link model, 2 particles in sector q: compare with best free-fermion guess
    (antisymmetrized single-particle energies from the 1-particle problem in the same sector) -> should fail"""
    m = build(N, J, h, 'std_bl')
    H = m['H'].toarray()
    Qops = [m['Q'][j].toarray() for j in range(N)]
    Ntot = sum(m['n']).toarray()
    # projector onto sector q and N_f=2 via diagonal? Q not diagonal in sz basis -> build projector
    P = np.eye(H.shape[0])
    for j in range(N):
        P = P @ (np.eye(H.shape[0]) + q[j] * Qops[j]) / 2
    evN, vN = np.linalg.eigh(Ntot)
    P2 = vN[:, np.abs(evN - 2) < 1e-9]; P2 = P2 @ P2.T
    P1 = vN[:, np.abs(evN - 1) < 1e-9]; P1 = P1 @ P1.T
    def sec_spec(Pn):
        PP = Pn @ P @ Pn
        w, v = np.linalg.eigh(PP); B = v[:, w > 0.5]
        return np.sort(np.linalg.eigvalsh(B.T @ H @ B))
    e1 = sec_spec(P1)
    P0 = vN[:, np.abs(evN) < 1e-9]; P0 = P0 @ P0.T
    e0 = sec_spec(P0)
    e2 = sec_spec(P2)
    # free guess: E2 = e0 + (e1_a - e0) + (e1_b - e0), a<b
    guess = np.sort([e1[a] + e1[b] - e0[0] for a in range(len(e1)) for b in range(a + 1, len(e1))])
    return e2, guess

if __name__ == '__main__':
    J, h = 1.0, 0.7
    print('--- spectra: full LGT vs union of mapped free-fermion sectors ---')
    for var, Ns in (('open_bl', [2, 3, 4, 5, 6]), ('open', [2, 3, 4, 5, 6]), ('ring', [3, 4, 5, 6])):
        for N in Ns:
            d1, d2, err = check_spectrum(N, J, h, var)
            print(f'{var:8s} N={N}: dim={d1}, predicted={d2}, max|dE|={err:.2e}')
    for N in [4, 5]:
        print(f'ring N={N}: alternative (JW-type) twist hypothesis error = {check_ring_alt(N, J, h):.3f}  (should be O(1))')
    print('--- dynamics, open_bl, CDW initial state, all observables incl. string/Loschmidt ones ---')
    for N in [4, 5, 6]:
        occ = set(range(0, N, 2))
        times = np.linspace(0, 6.0, 13)
        errs, rows = check_dynamics(N, J, h, times, occ)
        print(f'N={N}: max error over t in [0,6] = {max(errs):.2e}')
        if N == 6:
            for r in rows[::4]:
                print('   t=%.1f n=%s sx=%s sz=%s' % r)
    print('--- standard single-link Z2 model, one-particle sector: random-walk potential ---')
    for N in [3, 4, 5, 6]:
        d1, d2, err = check_std_oneparticle(N, J, h)
        print(f'std_bl N={N}: dim(N_f=1)={d1}, predicted={d2}, max|dE|={err:.2e}')
    e2, guess = is_free_spectrum_std(4, J, h, (1, 1, 1, 1))
    print('std_bl N=4, sector q=+1, N_f=2: exact', e2.round(3))
    print('                 free-fermion guess', guess.round(3), ' (mismatch => interacting)')
