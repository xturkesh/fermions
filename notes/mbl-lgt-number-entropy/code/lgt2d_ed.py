"""ED check of the 2D Z2 LGT with star electric term on a 2x3 open lattice (6 sites, 7 links, 2 plaquettes):
H = -J sum_<ij> sz_ij (f_i^dag f_j + h.c.) - h sum_i prod_{j~i} sx_ij
Claim: H = (+)_{q, W} [dGamma(h_{q,W}) on parity sector prod q = (-1)^Nf]  - h sum q,
h_{q,W} = -J sum s_ij (|i><j| + h.c.) + 2h diag(q), link signs s realizing plaquette fluxes W_p."""
import numpy as np, itertools
import scipy.sparse as sp
from lgt1d_ed import op_on, Zm, Xm, am
J, h = 1.0, 0.8
sites = [(r, c) for r in range(2) for c in range(3)]; sid = {s: i for i, s in enumerate(sites)}; N = 6
links = [((r, c), (r, c + 1)) for r in range(2) for c in range(2)] + [((0, c), (1, c)) for c in range(3)]
M = len(links); nq = N + M
f = []
for j in range(N):
    d = {k: Zm for k in range(j)}; d[j] = am; f.append(op_on(d, nq))
n = [f[j].T @ f[j] for j in range(N)]
sz = [op_on({N + b: Zm}, nq) for b in range(M)]; sx = [op_on({N + b: Xm}, nq) for b in range(M)]
D = 2 ** nq; Id = sp.identity(D, format='csr')
H = sp.csr_matrix((D, D))
for b, (a, c) in enumerate(links):
    i, j = sid[a], sid[c]; hop = f[i].T @ f[j]
    H = H - J * (sz[b] @ (hop + hop.T))
for i, s in enumerate(sites):
    P = Id
    for b, (a, c) in enumerate(links):
        if s in (a, c): P = P @ sx[b]
    H = H - h * P
Ntot = sum(n).diagonal()
Ev = []
for k in range(N + 1):
    ids = np.where(np.abs(Ntot - k) < 1e-9)[0]
    Ev += list(np.linalg.eigvalsh(H[ids][:, ids].toarray()))
Ev = np.sort(Ev)
# prediction
signlink = {links.index(((0, 0), (1, 0))): 'W1', links.index(((0, 2), (1, 2))): 'W2'}
pred = []
for q in itertools.product([1, -1], repeat=N):
    p = 0 if np.prod(q) == 1 else 1
    for W1, W2 in itertools.product([1, -1], repeat=2):
        A = np.diag(2 * h * np.array(q, float))
        for b, (a, c) in enumerate(links):
            s = 1.0
            if signlink.get(b) == 'W1': s = W1
            if signlink.get(b) == 'W2': s = W2
            A[sid[a], sid[c]] = A[sid[c], sid[a]] = -J * s
        eps = np.linalg.eigvalsh(A)
        for S in itertools.product([0, 1], repeat=N):
            if sum(S) % 2 == p:
                pred.append(float(np.dot(S, eps)) - h * sum(q))
pred = np.sort(pred)
print(f'2D 2x3: dim={len(Ev)}, predicted={len(pred)}, max|dE|={np.abs(Ev - pred).max():.2e}')
