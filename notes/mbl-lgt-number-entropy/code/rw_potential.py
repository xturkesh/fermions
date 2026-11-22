"""One matter fermion in the STANDARD 1D Z2 LGT (single-link electric term -h sum sx_b), uniform random sector:
exact reduction (verified by ED in lgt1d_ed.py) to  h = -J Lap + V,  V(x) = 2h sum_{b>=x} Q_b, Q_b iid +-1.
Compare eigenfunction correlator decay E Q(x0,x0+d) with the SKKM (iid Bernoulli) model.  Heuristic expectation:
superexponential decay ~ exp(-c d log d) away from returns of the walk."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
rng = np.random.default_rng(7)
J = 1.0
for h in [0.1, 0.25, 0.5]:
    for N in [400, 800]:
        x0 = N // 2; S = 300
        Qrw = np.zeros(N); Qbe = np.zeros(N); ipr_rw = 0; ipr_be = 0
        for s in range(S):
            Qb = rng.choice([-1, 1], size=N)
            V = 2 * h * np.cumsum(Qb[::-1])[::-1]
            V = V - V[x0]                       # energy reference irrelevant
            E, U = eigh_tridiagonal(V, -J * np.ones(N - 1))
            A = np.abs(U); Qrw += A @ A[x0] / S
            ipr_rw += (U[x0] ** 2 * (U ** 4).sum(0)).sum() / S   # spectral-weight-averaged IPR seen from x0
            q = rng.choice([-1, 1], size=N)
            E2, U2 = eigh_tridiagonal(2 * h * q.astype(float), -J * np.ones(N - 1))
            A2 = np.abs(U2); Qbe += A2 @ A2[x0] / S
            ipr_be += (U2[x0] ** 2 * (U2 ** 4).sum(0)).sum() / S
        ds = [1, 2, 4, 8, 16, 32, 64]
        print(f'h/J={h} N={N}: E Q(x0,x0+d), d={ds}')
        print('   random-walk V :', ' '.join(f'{Qrw[x0+d]:.1e}' for d in ds), f'  <IPR>={ipr_rw:.3f}')
        print('   iid Bernoulli :', ' '.join(f'{Qbe[x0+d]:.1e}' for d in ds), f'  <IPR>={ipr_be:.3f}')
