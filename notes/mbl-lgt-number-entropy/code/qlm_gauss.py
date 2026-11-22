"""Gauss-law kinematic bound on number entropy in U(1) quantum link models (staggered fermions, spin-S links).
In a fixed background-charge sector, N_A = E_R - E_L + #odd(A) - sum_A Q, so N_A takes at most 4S+1 values
(2S+1 if one end of A is an open boundary with fixed flux).  Check: random states in a sector."""
import numpy as np, itertools
rng = np.random.default_rng(0)
def check(N, S, Q, E0, A):
    Evals = np.arange(-S, S + 1)
    states = []
    for n in itertools.product([0, 1], repeat=N):
        # Gauss: E_j - E_{j-1} = rho_j + Q_j, rho_j = n_j - (1 - (-1)^j)/2 ; E_{-1} = E0 fixed
        E = []; prev = E0; ok = True
        for j in range(N):
            rho = n[j] - (1 - (-1) ** j) / 2
            cur = prev + rho + Q[j]
            if cur < -S - 1e-9 or cur > S + 1e-9: ok = False; break
            E.append(cur); prev = cur
        if ok: states.append((n, tuple(E)))
    dim = len(states)
    psi = rng.normal(size=dim) + 1j * rng.normal(size=dim); psi /= np.linalg.norm(psi)
    NA = np.array([sum(n[j] for j in A) for n, E in states])
    vals = sorted(set(NA[np.abs(psi) > 0]))
    p = np.array([np.sum(np.abs(psi[NA == v]) ** 2) for v in vals])
    SN = -(p * np.log(p)).sum()
    return dim, vals, SN
for S in [0.5, 1.0, 1.5]:
    N = 10
    for trial in range(3):
        Q = rng.choice([-1, 0, 1], size=N) if trial else np.zeros(N, int)
        E0 = -S if S != 0.5 else -0.5
        for A, name in [(range(0, 5), 'half chain (left end at fixed-flux boundary)'), (range(3, 7), 'interior interval')]:
            dim, vals, SN = check(N, S, Q, E0 if S == 0.5 else 0.0 if S == 1.0 else 0.5, A)
            if dim == 0: continue
            bound = np.log(2 * S + 1) if A[0] == 0 else np.log(4 * S + 1)
            print(f'S={S} Q={list(Q)} {name}: sector dim={dim}, #values of N_A={len(vals)}, S_N={SN:.3f} <= {bound:.3f}')
