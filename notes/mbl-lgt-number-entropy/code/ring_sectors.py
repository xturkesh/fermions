"""Sector-resolved check for the ring: project onto joint eigenspaces of Gauss operators Q_j and the
Wilson loop W = prod sz, compare with dGamma(h_{q,Wt}) on parity sector prod q = (-1)^{N_f},
for the two hypotheses Wt = W and Wt = W*(-1)^{N_f+1}."""
import numpy as np, itertools
from lgt1d_ed import build, hq, mb_energies
J, h = 1.0, 0.7
for N in [3, 4, 5]:
    m = build(N, J, h, 'ring'); D = 2 ** m['nq']
    H = m['H'].toarray()
    Wop = np.eye(D)
    for b in range(N): Wop = Wop @ m['sz'][b].toarray()
    Qs = [m['Q'][j].toarray() for j in range(N)]
    errA = errB = 0.0; tot = 0
    for q in itertools.product([1, -1], repeat=N):
        P = np.eye(D)
        for j in range(N): P = P @ (np.eye(D) + q[j] * Qs[j]) / 2
        for W in (1, -1):
            PW = P @ (np.eye(D) + W * Wop) / 2
            w, v = np.linalg.eigh(PW); B = v[:, w > 0.5]
            if B.shape[1] == 0: continue
            tot += B.shape[1]
            Es = np.sort(np.linalg.eigvalsh(B.T @ H @ B))
            p = 0 if np.prod(q) == 1 else 1
            const = -h * sum(q)
            pa = np.sort([e + const for e in mb_energies(np.linalg.eigvalsh(hq(N, J, h, q, True, W)), p)])
            pb = np.sort([e + const for e in mb_energies(np.linalg.eigvalsh(hq(N, J, h, q, True, W * (-1) ** (p + 1))), p)])
            errA = max(errA, np.abs(Es - pa).max() if len(pa) == len(Es) else np.inf)
            errB = max(errB, np.abs(Es - pb).max() if len(pb) == len(Es) else np.inf)
    print(f'ring N={N}: total dim covered {tot}/{D}; hypothesis twist=W: err={errA:.2e}; twist=W(-1)^(Nf+1): err={errB:.2e}')
