"""
Numerical checks for Section "Two dimensions" (sections/twod.tex).

Model (SKKM 2018, arXiv:1803.06574, Eqs. (1)-(2), with J_jk = J, h_j = h):
    H = -J sum_{b={i,j}} sz_b (f_i^dag f_j + h.c.) - h sum_i A_i,     A_i = prod_{b incident to i} sx_b,
    G_i = (-1)^{n_i} A_i,   W_gamma = prod_{b in gamma} sz_b  (gamma a cycle).
Everything is built from scratch in the occupation (x) sigma^z basis; nothing of the mapping is built in.

A. Theorem thm:2dmap, SECTOR BY SECTOR: H restricted to each joint eigenspace of all G_i and of the
   fundamental-cycle Wilson loops is compared with dGamma(h_{q,s}) - h sum q on the parity sector
   prod q = (-1)^{N_f}.  Graphs: planar 2x3 box, 2x2 torus (multigraph; plaquettes + the two
   non-contractible loops), K4.  Also the sector dimensions and plaquette/loop relations.
B. Proposition prop:2dlaws (sector laws + dynamics): <n_x(t)> by sparse time evolution of the full LGT
   versus the sector-averaged free-fermion prediction, for
     (i)  sigma^z-polarized links x Fock matter  -> q uniform on the parity hyperplane, zero flux;
     (ii) sigma^x-polarized links x checkerboard -> q = (-1)^{n^0} deterministic, fluxes uniform.
   Graphs: 2x4 box, 3x3 box, 2x2 torus.  Sector weights ||P_{q,phi} Psi_0||^2 checked directly.
C. Checkerboard identity h^2 = K^2 + 4h^2, gap |E| >= 2|h|, propagator formula, and the K^2 entries
   J^2 s s'(1 + W_p) on diagonal next-nearest neighbours (random pi-flux on a 16x16 box).
D. Proposition prop:cluster: ||e^{-ith_q} - e^{-ith^0_q}|| versus the bound 8J/g + 48 t J^2/g,
   and the probability to leave the equal-charge cluster, versus t (24x24 box).
"""
import itertools, time, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm, eigh

rng = np.random.default_rng(20261122)
T0 = time.time()


# ----------------------------------------------------------------------------------------------
# graphs
# ----------------------------------------------------------------------------------------------
def box(L1, L2):
    sites = [(x, y) for y in range(L2) for x in range(L1)]
    sid = {s: i for i, s in enumerate(sites)}
    edges = []
    for y in range(L2):
        for x in range(L1 - 1):
            edges.append((sid[(x, y)], sid[(x + 1, y)]))
    for y in range(L2 - 1):
        for x in range(L1):
            edges.append((sid[(x, y)], sid[(x, y + 1)]))
    plaq = []
    for y in range(L2 - 1):
        for x in range(L1 - 1):
            a, b, c, d = sid[(x, y)], sid[(x + 1, y)], sid[(x + 1, y + 1)], sid[(x, y + 1)]
            plaq.append([edges.index((a, b)), edges.index((b, c)), edges.index((d, c)), edges.index((a, d))])
    eps = [(-1) ** (x + y) for (x, y) in sites]
    return dict(name=f'box {L1}x{L2}', NV=len(sites), edges=edges, plaq=plaq, loops=[], eps=eps, sites=sites)


def torus(L1, L2):
    """(Z/L1) x (Z/L2) with edges h_{x,y}: (x,y)-(x+1,y) and v_{x,y}: (x,y)-(x,y+1); multigraph if L=2."""
    sites = [(x, y) for y in range(L2) for x in range(L1)]
    sid = {s: i for i, s in enumerate(sites)}
    edges, hid, vid = [], {}, {}
    for y in range(L2):
        for x in range(L1):
            hid[(x, y)] = len(edges); edges.append((sid[(x, y)], sid[((x + 1) % L1, y)]))
    for y in range(L2):
        for x in range(L1):
            vid[(x, y)] = len(edges); edges.append((sid[(x, y)], sid[(x, (y + 1) % L2)]))
    plaq = [[hid[(x, y)], vid[((x + 1) % L1, y)], hid[(x, (y + 1) % L2)], vid[(x, y)]]
            for y in range(L2) for x in range(L1)]
    loops = [[hid[(x, 0)] for x in range(L1)], [vid[(0, y)] for y in range(L2)]]
    eps = [(-1) ** (x + y) for (x, y) in sites]
    return dict(name=f'torus {L1}x{L2}', NV=len(sites), edges=edges, plaq=plaq, loops=loops, eps=eps, sites=sites)


def K4():
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    return dict(name='K4', NV=4, edges=edges, plaq=[], loops=[], eps=None, sites=list(range(4)))


def fundamental_cycles(NV, edges):
    """spanning tree by BFS from 0; returns tree mask and, for each non-tree edge, its cycle (edge list)."""
    adj = {i: [] for i in range(NV)}
    for b, (i, j) in enumerate(edges):
        adj[i].append((j, b)); adj[j].append((i, b))
    parent = {0: (None, None)}; order = [0]; k = 0
    while k < len(order):
        v = order[k]; k += 1
        for (u, b) in adj[v]:
            if u not in parent:
                parent[u] = (v, b); order.append(u)
    assert len(parent) == NV, 'graph not connected'
    tree = set(b for (v, b) in parent.values() if b is not None)

    def path_to_root(v):
        out = []
        while parent[v][0] is not None:
            out.append(parent[v][1]); v = parent[v][0]
        return out
    cycles = []
    for b, (i, j) in enumerate(edges):
        if b in tree:
            continue
        pi, pj = path_to_root(i), path_to_root(j)
        cyc = set(pi) ^ set(pj) ^ {b}
        cycles.append(sorted(cyc))
    return tree, cycles


def signs_for_flux(NE, tree, cycles, w):
    """s = +1 on tree edges, s_b = w(C_b) on the non-tree edge b of fundamental cycle C_b."""
    s = np.ones(NE)
    k = 0
    for b in range(NE):
        if b not in tree:
            s[b] = w[k]; k += 1
    return s


def hop_matrix(NV, edges, s, J):
    K = np.zeros((NV, NV))
    for b, (i, j) in enumerate(edges):
        K[i, j] += -J * s[b]; K[j, i] += -J * s[b]
    return K


# ----------------------------------------------------------------------------------------------
# many-body LGT in the occupation (x) sigma^z basis, restricted to fixed N_f
# basis state = (n, s): n in {0,1}^NV (bit i <-> site i), s in {0,1}^NE (bit b = 1 <-> spin down)
# ----------------------------------------------------------------------------------------------
def basis_Nf(NV, NE, Nf):
    nconf = [sum(1 << i for i in occ) for occ in itertools.combinations(range(NV), Nf)]
    nconf = np.array(sorted(nconf), dtype=np.int64)
    sconf = np.arange(1 << NE, dtype=np.int64)
    states = (nconf[:, None] << NE | sconf[None, :]).ravel()
    return np.sort(states)


def popcount(a):
    a = a.copy(); c = np.zeros_like(a)
    while np.any(a):
        c += a & 1; a >>= 1
    return c


def build_ops(G, J, h, Nf, want_G=True):
    NV, edges = G['NV'], G['edges']; NE = len(edges)
    st = basis_Nf(NV, NE, Nf); D = len(st)
    idx = lambda arr: np.searchsorted(st, arr)
    n = st >> NE; s = st & ((1 << NE) - 1)
    sz = 1 - 2 * ((s[:, None] >> np.arange(NE)[None, :]) & 1)       # D x NE, sigma^z eigenvalues
    rows, cols, vals = [], [], []
    for b, (i, j) in enumerate(edges):
        for (a, c) in ((i, j), (j, i)):                              # f_a^dag f_c
            m = (((n >> c) & 1) == 1) & (((n >> a) & 1) == 0)
            src = np.where(m)[0]
            nn = n[src]
            sign1 = (-1) ** popcount(nn & ((1 << c) - 1))              # f_c: (-1)^{N_<c(n)}
            n1 = nn & ~(1 << c)
            sign2 = (-1) ** popcount(n1 & ((1 << a) - 1))              # f_a^dag: (-1)^{N_<a(n1)}
            n2 = n1 | (1 << a)
            tgt = idx((n2 << NE) | s[src])
            rows.append(tgt); cols.append(src); vals.append(-J * sz[src, b] * sign1 * sign2)
    inc = [[b for b, e in enumerate(edges) if i in e] for i in range(NV)]
    Amask = [sum(1 << b for b in inc[i]) for i in range(NV)]
    for i in range(NV):
        tgt = idx((n << NE) | (s ^ Amask[i]))
        rows.append(tgt); cols.append(np.arange(D)); vals.append(-h * np.ones(D))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))
    out = dict(H=H, st=st, n=n, s=s, sz=sz, D=D, NE=NE)
    if want_G:
        Gs = []
        for i in range(NV):
            tgt = idx((n << NE) | (s ^ Amask[i]))
            par = 1 - 2 * ((n >> i) & 1)
            Gs.append(sp.csr_matrix((par.astype(float), (tgt, np.arange(D))), shape=(D, D)))
        out['G'] = Gs
    return out


def Wdiag(sz, cyc):
    return np.prod(sz[:, cyc], axis=1)


# ----------------------------------------------------------------------------------------------
# A. sector-resolved spectra
# ----------------------------------------------------------------------------------------------
def check_sectors(G, J, h):
    NV, edges = G['NV'], G['edges']; NE = len(edges)
    tree, cycles = fundamental_cycles(NV, edges)
    beta = NE - NV + 1
    assert len(cycles) == beta
    maxerr = 0.0; nsec = 0; dimerr = 0; total = 0
    cG = rng.normal(size=NV); cW = rng.normal(size=beta)
    allkeys = set(); countok = True
    for Nf in range(NV + 1):
        ops = build_ops(G, J, h, Nf)
        D = ops['D']; total += D
        Wd = [Wdiag(ops['sz'], cyc) for cyc in cycles]
        X = sum(cG[i] * ops['G'][i] for i in range(NV)).toarray() + np.diag(sum(cW[k] * Wd[k] for k in range(beta)))
        ev, V = eigh(X)
        Hd = ops['H'].toarray()
        # labels of each eigenvector
        qlab = np.array([np.einsum('ij,ij->j', V, ops['G'][i] @ V) for i in range(NV)]).T
        wlab = np.array([np.einsum('ij,ij->j', V, Wd[k][:, None] * V) for k in range(beta)]).T
        lab = np.hstack([qlab, wlab])
        assert np.abs(np.abs(lab) - 1).max() < 1e-8
        lab = np.sign(lab).astype(int)
        keys = {}
        for col, key in enumerate(map(tuple, lab)):
            keys.setdefault(key, []).append(col)
        allkeys |= set(keys)
        countok &= (len(keys) == 2 ** (NV - 1 + beta))
        for key, cols in keys.items():
            q = np.array(key[:NV]); w = np.array(key[NV:])
            B = V[:, cols]
            Es = np.sort(eigh(B.T @ Hd @ B, eigvals_only=True))
            s = signs_for_flux(NE, tree, cycles, w)
            eps = np.linalg.eigvalsh(hop_matrix(NV, edges, s, J) + np.diag(2 * h * q))
            pred = np.sort([eps[list(c)].sum() for c in itertools.combinations(range(NV), Nf)]) - h * q.sum()
            if np.prod(q) != (-1) ** Nf:
                dimerr += 1
                continue
            if len(pred) != len(Es):
                dimerr += 1
                continue
            maxerr = max(maxerr, np.abs(Es - pred).max()); nsec += 1
    # relations among plaquettes/loops in terms of the fundamental cycles: check that the plaquette and
    # loop edge sets are cycles (even degree), and the rank of plaquettes (+ loops) over F2
    def is_cycle(cyc):
        deg = np.zeros(NV, int)
        for b in cyc:
            i, j = edges[b]; deg[i] += 1; deg[j] += 1
        return np.all(deg % 2 == 0)
    allc = G['plaq'] + G['loops']
    ok = all(is_cycle(c) for c in allc)
    M = np.zeros((len(allc), NE), dtype=np.uint8)
    for r, c in enumerate(allc):
        M[r, c] = 1
    rank = f2rank(M)
    rank_pl = f2rank(M[:len(G['plaq'])]) if G['plaq'] else 0
    return dict(total=total, nsec=nsec, maxerr=maxerr, dimerr=dimerr, beta=beta, cycles_ok=ok,
                rank_all=rank, rank_plaq=rank_pl, nplaq=len(G['plaq']), nloops=len(G['loops']),
                nkeys=len(allkeys), countok=countok)


def f2rank(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = [i for i in range(r, rows) if M[i, c]]
        if not piv:
            continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


print('=' * 100)
print('A. Theorem thm:2dmap, sector by sector (J=1, h=0.8)')
print('=' * 100)
J, hh = 1.0, 0.8
for G in [box(2, 3), torus(2, 2), K4()]:
    r = check_sectors(G, J, hh)
    print(f"{G['name']:12s}: dim={r['total']}, beta=dim Z1={r['beta']}, distinct sectors (q,phi) found={r['nkeys']} "
          f"(predicted 2^(NV+beta)={2 ** (G['NV'] + r['beta'])}); in every N_f block exactly 2^(NV-1+beta) sectors, "
          f"all of parity prod q=(-1)^N_f and dimension C(NV,N_f): {r['countok'] and r['dimerr'] == 0}; "
          f"max|dE| over {r['nsec']} (sector,N_f) blocks = {r['maxerr']:.2e}")
    if G['plaq'] or G['loops']:
        print(f"              plaquettes={r['nplaq']}, loops={r['nloops']}: all are cycles: {r['cycles_ok']}; "
              f"F2-rank(plaquettes)={r['rank_plaq']}, F2-rank(plaquettes+loops)={r['rank_all']}")
print(f'[time {time.time() - T0:.1f}s]')


# ----------------------------------------------------------------------------------------------
# B. sector laws and dynamics
# ----------------------------------------------------------------------------------------------
def free_density(NV, edges, J, h, q, s, occ, times):
    Hq = hop_matrix(NV, edges, s, J) + np.diag(2 * h * np.asarray(q, float))
    E, U = np.linalg.eigh(Hq)
    out = []
    for t in times:
        Ut = (U * np.exp(-1j * t * E)) @ U.T
        out.append((np.abs(Ut[:, occ]) ** 2).sum(axis=1))
    return np.array(out)


def check_dynamics(G, J, h, occ, gauge, times):
    NV, edges = G['NV'], G['edges']; NE = len(edges)
    Nf = len(occ)
    ops = build_ops(G, J, h, Nf)
    st, D = ops['st'], ops['D']
    n0 = sum(1 << i for i in occ)
    psi0 = np.zeros(D, complex)
    if gauge == 'z':
        psi0[np.searchsorted(st, (n0 << NE) | 0)] = 1.0
    else:
        allS = (n0 << NE) | np.arange(1 << NE)
        psi0[np.searchsorted(st, allS)] = 2.0 ** (-NE / 2)
    psis = expm_multiply(-1j * ops['H'], psi0, start=times[0], stop=times[-1], num=len(times), endpoint=True)
    nocc = ((ops['n'][:, None] >> np.arange(NV)[None, :]) & 1).astype(float)
    ED = np.array([(np.abs(p) ** 2) @ nocc for p in psis])
    tree, cycles = fundamental_cycles(NV, edges)
    beta = len(cycles)
    pred = np.zeros_like(ED)
    # sector weights directly from projectors
    Wd = [Wdiag(ops['sz'], cyc) for cyc in cycles]
    werr = 0.0
    if gauge == 'z':
        par = (-1) ** Nf
        qs = [q for q in itertools.product([1, -1], repeat=NV) if np.prod(q) == par]
        for q in qs:
            pred += free_density(NV, edges, J, h, q, np.ones(NE), occ, times) / len(qs)
        # weights: W_gamma Psi0 = Psi0 and ||P_q Psi0||^2 = 2^{-(NV-1)} on the hyperplane, 0 off it
        werr = max(werr, max(np.abs(np.vdot(psi0, Wd[k] * psi0) - 1) for k in range(beta)) if beta else 0.0)
        for q in itertools.product([1, -1], repeat=NV):
            v = psi0.copy()
            for i in range(NV):
                v = 0.5 * (v + q[i] * (ops['G'][i] @ v))
            target = 2.0 ** (-(NV - 1)) if np.prod(q) == par else 0.0
            werr = max(werr, abs(np.vdot(v, v).real - target))
    else:
        q0 = [(-1) ** (1 if i in occ else 0) for i in range(NV)]
        for w in itertools.product([1, -1], repeat=beta):
            s = signs_for_flux(NE, tree, cycles, w)
            pred += free_density(NV, edges, J, h, q0, s, occ, times) / 2 ** beta
        # weights: G_i Psi0 = q0_i Psi0 ; ||P_phi Psi0||^2 = 2^{-beta} for every flux phi
        werr = max(np.abs(np.vdot(psi0, ops['G'][i] @ psi0) - q0[i]) for i in range(NV))
        lab = np.array([Wd[k] for k in range(beta)]).T if beta else np.zeros((D, 0))
        for w in itertools.product([1, -1], repeat=beta):
            m = np.all(lab == np.array(w)[None, :], axis=1)
            werr = max(werr, abs((np.abs(psi0[m]) ** 2).sum() - 2.0 ** (-beta)))
    return np.abs(ED - pred).max(), werr, ED, pred


print()
print('=' * 100)
print('B. Sector laws and dynamics: full LGT time evolution vs sector-averaged free fermions')
print('=' * 100)
times = np.linspace(0.0, 8.0, 17)
for (G, occ, label) in [(box(2, 4), None, 'checkerboard'), (box(3, 3), None, 'checkerboard'),
                        (torus(2, 2), None, 'checkerboard'), (box(2, 4), [0, 1, 6], 'Fock {0,1,6}')]:
    if occ is None:
        occ = [i for i in range(G['NV']) if G['eps'][i] == 1]
    for gauge in ('z', 'x'):
        for h in (0.4, 1.3):
            t1 = time.time()
            err, werr, ED, pred = check_dynamics(G, 1.0, h, occ, gauge, times)
            print(f"{G['name']:10s} matter={label:13s} links={'sigma^z-pol' if gauge == 'z' else 'sigma^x-pol'} "
                  f"h={h}: max_t,x |<n_x>_ED - <n_x>_pred| = {err:.2e}; sector-weight error = {werr:.2e} "
                  f"[{time.time() - t1:.1f}s]")
print(f'[time {time.time() - T0:.1f}s]')


# ----------------------------------------------------------------------------------------------
# C. checkerboard identity
# ----------------------------------------------------------------------------------------------
print()
print('=' * 100)
print('C. Checkerboard state in a random pi-flux sector: h_{q,s} = K_s - 2h Gamma')
print('=' * 100)
L = 16; G = box(L, L); NV = G['NV']; edges = G['edges']; NE = len(edges)
eps = np.array(G['eps'], float); Gam = np.diag(eps)
sid = {s: i for i, s in enumerate(G['sites'])}
eid = {e: b for b, e in enumerate(edges)}
for h in (0.3, 1.0, 3.0):
    s = rng.choice([-1.0, 1.0], size=NE)          # uniform signs = uniform fluxes (Lemma)
    K = hop_matrix(NV, edges, s, 1.0)
    q0 = -eps                                     # occupied on eps=+1 sites: q = (-1)^n = -eps
    Hq = K + np.diag(2 * h * q0)
    e1 = np.abs(Hq @ Hq - K @ K - 4 * h * h * np.eye(NV)).max()
    E = np.linalg.eigvalsh(Hq); eK = np.linalg.eigvalsh(K)
    gap = np.abs(E).min()
    mult = np.abs(np.sort(np.abs(E)) - np.sort(np.sqrt(eK ** 2 + 4 * h * h))).max()
    t = 2.7
    Ut = expm(-1j * t * Hq)
    w, Vk = np.linalg.eigh(K @ K + 4 * h * h * np.eye(NV)); R = np.sqrt(w)
    cosM = (Vk * np.cos(t * R)) @ Vk.T; sincM = (Vk * (np.sin(t * R) / R)) @ Vk.T
    e3 = np.abs(Ut - (cosM - 1j * Hq @ sincM)).max()
    # K^2 entries
    K2 = K @ K; e4 = 0.0
    for y in range(L - 1):
        for x in range(L - 1):
            i, j1, k, j2 = sid[(x, y)], sid[(x + 1, y)], sid[(x + 1, y + 1)], sid[(x, y + 1)]
            sa = s[eid[(i, j1)]] * s[eid[(j1, k)]]; sb = s[eid[(i, j2)]] * s[eid[(j2, k)]]
            Wp = sa * sb
            e4 = max(e4, abs(K2[i, k] - sa * (1 + Wp)))
    print(f'h={h}: |h^2-K^2-4h^2|={e1:.1e}; min|E|={gap:.4f} >= 2h={2 * h}; |sorted|E| - sqrt(eps_K^2+4h^2)|={mult:.1e}; '
          f'propagator formula err={e3:.1e}; K^2 diag-NNN entries vs J^2 s s\'(1+W_p) err={e4:.1e}')
print(f'[time {time.time() - T0:.1f}s]')


# ----------------------------------------------------------------------------------------------
# D. cluster confinement
# ----------------------------------------------------------------------------------------------
print()
print('=' * 100)
print('D. Strong-coupling cluster confinement on a 24x24 box (J=1), q i.i.d. Bernoulli(1/2), zero flux')
print('=' * 100)
L = 24; G = box(L, L); NV = G['NV']; edges = G['edges']; NE = len(edges)
K = hop_matrix(NV, edges, np.ones(NE), 1.0)
tlist = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]
coords = np.array(G['sites'])
bulk = [i for i, (x, y) in enumerate(G['sites']) if 6 <= x < L - 6 and 6 <= y < L - 6]


def clusters(q):
    lab = -np.ones(NV, int); c = 0
    nbr = {i: [] for i in range(NV)}
    for (i, j) in edges:
        if q[i] == q[j]:
            nbr[i].append(j); nbr[j].append(i)
    for i in range(NV):
        if lab[i] >= 0:
            continue
        stack = [i]; lab[i] = c
        while stack:
            v = stack.pop()
            for u in nbr[v]:
                if lab[u] < 0:
                    lab[u] = c; stack.append(u)
        c += 1
    return lab


for h in (3.0, 4.0, 6.0):
    g = 4 * h - 8.0
    nsamp = 6
    normdiff = np.zeros(len(tlist)); pout = np.zeros(len(tlist)); pout_x = np.zeros(len(tlist))
    for smp in range(nsamp):
        q = rng.choice([-1, 1], size=NV)
        Hq = K + np.diag(2 * h * q)
        mask = (q[:, None] == q[None, :]).astype(float)
        H0 = Hq * mask
        E, U = np.linalg.eigh(Hq); E0, U0 = np.linalg.eigh(H0)
        lab = clusters(q)
        for it, t in enumerate(tlist):
            Ut = (U * np.exp(-1j * t * E)) @ U.T
            U0t = (U0 * np.exp(-1j * t * E0)) @ U0.T
            normdiff[it] = max(normdiff[it], np.linalg.norm(Ut - U0t, 2))
            # probability to be outside the own cluster, averaged over bulk starting sites
            P = np.abs(Ut[:, bulk]) ** 2
            inside = (lab[:, None] == lab[bulk][None, :])
            pout[it] += (P * (~inside)).sum(axis=0).mean() / nsamp
    print(f'h={h} (g=4h-8J={g}):')
    for it, t in enumerate(tlist):
        bnd = 8.0 / g + 48.0 * t / g
        print(f'   t={t:7.1f}: max_q ||e^(-ith) - e^(-ith0)|| = {normdiff[it]:.3f}  (bound {bnd:9.3f}, trivial bound 2);'
              f'  E P(leave own cluster) = {pout[it]:.4f}')
    print(f'   time scale g/J^2 = {g:.0f}')
print(f'[time {time.time() - T0:.1f}s]')

# cluster-size statistics at p = 1/2 (subcritical site percolation): P(cluster radius >= R)
print()
Lp = 200; nrep = 20; Rs = [2, 4, 6, 8, 10, 12, 14]
cnt = np.zeros(len(Rs)); tot = 0
for rep in range(nrep):
    q = rng.choice([-1, 1], size=(Lp, Lp))
    # BFS from 200 random bulk sites
    for _ in range(200):
        x0, y0 = rng.integers(40, Lp - 40, size=2)
        seen = {(x0, y0)}; stack = [(x0, y0)]; R = 0
        while stack:
            x, y = stack.pop()
            R = max(R, abs(x - x0), abs(y - y0))
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                u = (x + dx, y + dy)
                if not (0 <= u[0] < Lp and 0 <= u[1] < Lp):
                    continue      # a cluster reaching the boundary has radius >= 40 > max(Rs) anyway
                if u not in seen and q[u] == q[x0, y0]:
                    seen.add(u); stack.append(u)
        tot += 1
        for k, Rr in enumerate(Rs):
            cnt[k] += (R >= Rr)
print('Empirical P(sup-norm radius of the equal-charge cluster of x >= R), q i.i.d. Bernoulli(1/2) on Z^2:')
print('   ' + '  '.join(f'R={Rr}: {c / tot:.2e}' for Rr, c in zip(Rs, cnt)))
print(f'[total time {time.time() - T0:.1f}s]')
