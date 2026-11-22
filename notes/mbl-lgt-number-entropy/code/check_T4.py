"""T4 sanity check: order-n quasi-conserved dressing of N_A for the random t-V (XXZ) chain.
H = H0 + J T,  H0 = sum h_x n_x + J*Dz sum n_x n_{x+1}  (diagonal),  T = sum_x (hop x<->x+1).
Recursion Y_0 = N_A, [H0, Y_k] = -J [T, Y_{k-1}]  (solved entrywise in the occupation basis),
Ntilde_n = sum_{k<=n} Y_k  satisfies  [H, Ntilde_n] = J [T, Y_n].
Checks: diag(J[T,Y_k]) == 0 (reality argument), eps_n = ||[H,Ntilde_n]|| decreasing,
delta_n = ||Ntilde_n - N_A||, and ||N_A(t)-N_A|| <= 2 delta_n + eps_n t.
"""
import numpy as np
from lbit import Sector

rng = np.random.default_rng(5)
L, J, Dz = 10, 1.0, 1.0
sec = Sector(L)
Dm = sec.D
T = np.zeros((Dm, Dm))
for i in range(L - 1):
    i00, i11, i01, i10 = sec.bond[i]
    T[i01, i10] = T[i10, i01] = 1.0
occ = np.array([[(s >> i) & 1 for i in range(L)] for s in sec.states], dtype=float)
NA = np.diag(sec.NA.astype(float))
for W in [6.0, 12.0, 24.0]:
    for rep in range(2):
        h = rng.uniform(-W, W, L)
        E0 = occ @ h + J*Dz*(occ[:, :-1]*occ[:, 1:]).sum(1)
        H = np.diag(E0) + J*T
        den = E0[:, None] - E0[None, :]
        Y = NA.copy(); Nt = NA.copy()
        out = []
        gam_min = np.inf
        for k in range(1, 9):
            R = J*(T @ Y - Y @ T)
            diagmax = np.abs(np.diag(R)).max()
            mask = np.abs(R) > 1e-14
            gam_min = min(gam_min, np.abs(den[mask]).min())
            Y = np.where(mask, -R/np.where(mask, den, 1.0), 0.0)
            Nt = Nt + Y
            C = H @ Nt - Nt @ H
            eps = np.abs(np.linalg.eigvalsh(1j*C)).max()
            dlt = np.abs(np.linalg.eigvalsh(Nt - NA)).max()
            out.append((k, eps, dlt, diagmax))
        best = min(out, key=lambda z: z[1])
        k, eps, dlt, _ = best
        w, v = np.linalg.eigh(H)
        worst = 0.0
        for t in [1, 10, 100, 1e3, 1e4, 1e5]:
            Ut = (v*np.exp(-1j*w*t)) @ v.conj().T
            dev = np.abs(np.linalg.eigvalsh(Ut.conj().T @ NA @ Ut - NA)).max()
            worst = max(worst, dev/(2*dlt + eps*t))
        print(f"W={W:5.1f} rep{rep}: min denominator met = {gam_min:.3f};  max|diag R_k| = {max(o[3] for o in out):.1e}")
        print("     n, eps_n=||[H,N~_n]||, delta_n=||N~_n-N_A||:",
              ", ".join(f"({o[0]}, {o[1]:.2e}, {o[2]:.3f})" for o in out))
        print(f"     best n={k}: max_t ||N_A(t)-N_A|| / (2 delta + eps t) over t in 1..1e5 = {worst:.3f} (<=1)")
