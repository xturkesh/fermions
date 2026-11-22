"""T1 tails: Gaussian U(1) dressing u = exp(i theta), theta_xy = g |x-y|^{-alpha} (x != y),
open chain of L sites, Neel state, cut in the middle.  Compare
  (i)  Var_{psi}(U^dag N_A U)  (static charge transferred by the dressing),
  (ii) dephased (infinite-time) Var for H = U D U^dag, D = sum eps_i n_i (random phases),
  (iii) ||Delta|| = ||U^dag N_A U - N_A|| = max(Tr K_+, Tr K_-), K = u^dag P_A u - P_A,
  (iv)  first-moment quantity sum_{x in A, y notin A} |theta_xy|  (T1 sufficient condition).
Expectation: (i),(ii) bounded iff ~ sum_r r theta_r^2 < inf (alpha > 1);
(iii) roughly iff alpha > 1; (iv) iff alpha > 2.
"""
import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(3)
g = 0.3
print(" alpha    L    Var_static  Var_dephased  ||Delta||   sum|theta_AB|")
for alpha in [0.75, 1.0, 1.5, 2.0, 3.0]:
    for L in [64, 128, 256, 512, 1024]:
        x = np.arange(L)
        R = np.abs(x[:, None] - x[None, :]).astype(float)
        theta = np.where(R > 0, g*np.where(R > 0, R, 1.0)**(-alpha), 0.0)
        u = expm(1j*theta)
        n = (x % 2 == 0).astype(float)
        A = x < L//2
        C = (u*n) @ u.conj().T                     # u C0 u^dag
        var_static = (np.abs(C[np.ix_(A, ~A)])**2).sum()
        vd = []
        for _ in range(8):
            ph = np.exp(-1j*rng.uniform(0, 2*np.pi, L))
            ut = (u*ph) @ u.conj().T
            Ct = (ut*n) @ ut.conj().T
            vd.append((np.abs(Ct[np.ix_(A, ~A)])**2).sum())
        PA = np.diag(A.astype(float))
        K = u.conj().T @ PA @ u - PA
        ev = np.linalg.eigvalsh((K + K.conj().T)/2)
        nD = max(ev[ev > 0].sum(), -ev[ev < 0].sum())
        fm = theta[np.ix_(A, ~A)].sum()
        print(f"{alpha:5.2f} {L:5d}   {var_static:9.4f}   {np.mean(vd):9.4f}    {nD:9.4f}   {fm:9.4f}")
