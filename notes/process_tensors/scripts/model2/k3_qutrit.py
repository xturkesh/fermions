import numpy as np, time
from kreplica import otoc_k
rng = np.random.default_rng(5)
Zm = (rng.normal(size=(9, 9)) + 1j*rng.normal(size=(9, 9)))/np.sqrt(2); Q, R = np.linalg.qr(Zm); U = Q*(np.diag(R)/abs(np.diag(R)))
M1 = np.zeros((9, 9), complex)                        # single-replica channel on A (qutrit)
for kk in range(9):
    E = np.zeros(9, complex); E[kk] = 1; Xa = E.reshape(3, 3)
    M1[:, kk] = (np.einsum('aibi->ab', (U.conj().T @ np.kron(Xa, np.eye(3)) @ U).reshape(3, 3, 3, 3))/3).reshape(-1)
lam = np.sort(np.abs(np.linalg.eigvals(M1)))[::-1][1]
A = np.diag([1., 1., -2.]).astype(complex); B = np.diag([1., 0., 0.]).astype(complex)
k3A, phiB = np.trace(A @ A @ A).real/3, 1/3
print("qutrit: lambda=%.4f ; free late value k3(A) phi(B)^3 = %.6f ; noisy -> 0" % (lam, k3A*phiB**3))
D = 10**4
for p in [0.0, 0.1]:
    t0 = time.time(); vals = {t: otoc_k(U, A, B, t, D, 3, p) for t in [12, 18, 24]}
    if p == 0: print("  p=0  C3(12,18,24) =", {t: round(v, 6) for t, v in vals.items()})
    else:
        r = abs(vals[24]/vals[18])**(1/6)
        print("  p=%.1f C3(12,18,24) = %s  decay rate=%.4f  (1-p)^3=%.4f (1-p)^2=%.4f  [%.0fs]" % (p, {t: round(v, 6) for t, v in vals.items()}, r, (1-p)**3, (1-p)**2, time.time()-t0))
