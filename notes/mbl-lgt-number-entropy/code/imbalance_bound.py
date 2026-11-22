"""Imbalance theorem, numerically: kappa = E sum_{x in B, y notin B} R(x,y) for an interior interval B (two cuts),
R(x,y) = sum_a phi_a(x)^2 phi_a(y)^2 (time-averaged transition probability), and the resulting rigorous lower bound
E Ibar >= h^2 / (16 l (J+h)^2), l = max(2, ceil(32 kappa^2)); compare with actual E Ibar. Also weak disorder, large N."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
rng = np.random.default_rng(11)
J = 1.0
for h in [0.1, 0.25, 0.5, 1.0]:
    for N in ([400, 1600] if h >= 0.25 else [800, 3200]):
        S = 40 if N <= 800 else 12
        kap = 0.0; Ib = 0.0
        a0, b0 = N // 4, 3 * N // 4
        for s in range(S):
            q = rng.choice([-1, 1], size=N)
            E, V = eigh_tridiagonal(2 * h * q.astype(float), -J * np.ones(N - 1))
            P = V ** 2
            inB = np.zeros(N, bool); inB[a0:b0] = True
            # sum_{x in B, y notin B} R = sum_a m_a(B) (1 - m_a(B))
            mB = P[inB].sum(0)
            kap += (mB * (1 - mB)).sum() / S
            st = (-1.0) ** np.arange(N)
            Ib += ((st @ P) ** 2).sum() / N / S
        l = max(2, int(np.ceil(32 * kap ** 2)))
        print(f'lambda=2h/J={2*h:4.2f} N={N:5d}: kappa(2 cuts)={kap:7.3f}  l={l:8d}  rigorous bound={h**2/(16*l*(J+h)**2):.2e}  actual E Ibar={Ib:.4f}')
