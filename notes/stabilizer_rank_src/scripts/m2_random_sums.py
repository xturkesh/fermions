# M2 of normalized sums of k (pseudo)random stabilizer states, n = 12.
# Heuristic claim: M2 ~ 3 log2 k  (well below saturation ~ n - 2).
import numpy as np
rng = np.random.default_rng(7)
n = 12; N = 2**n

def apply1(psi, U, q):
    psi = psi.reshape([2]*n); psi = np.moveaxis(psi, q, 0)
    psi = np.tensordot(U, psi, axes=(1, 0)); return np.moveaxis(psi, 0, q).reshape(N)
def cnot(psi, a, b):
    psi = psi.reshape([2]*n).copy(); idx = [slice(None)]*n; idx[a] = 1
    sub = psi[tuple(idx)]; bb = b if b < a else b - 1
    psi[tuple(idx)] = np.flip(sub, axis=bb); return psi.reshape(N)
Hd = np.array([[1, 1], [1, -1]])/np.sqrt(2); Sg = np.diag([1, 1j])
def rand_stab(depth=400):
    psi = np.zeros(N, complex); psi[0] = 1
    for _ in range(depth):
        r = rng.integers(3)
        if r == 0: psi = apply1(psi, Hd, rng.integers(n))
        elif r == 1: psi = apply1(psi, Sg, rng.integers(n))
        else:
            a, b = rng.choice(n, 2, replace=False); psi = cnot(psi, a, b)
    return psi
def wht(M):                       # Walsh-Hadamard along last axis (unnormalized)
    B = M.shape[0]; M = M.reshape([B] + [2]*n)
    for ax in range(1, n+1):
        a0 = np.take(M, 0, axis=ax); a1 = np.take(M, 1, axis=ax)
        M = np.stack([a0 + a1, a0 - a1], axis=ax)
    return M.reshape(B, N)
def M2(psi):
    x = np.arange(N); tot = 0.0
    for a0 in range(0, N, 256):
        a = np.arange(a0, a0 + 256)[:, None]
        G = np.conj(psi[x[None, :] ^ a]) * psi[None, :]
        tot += np.sum(np.abs(wht(G))**4)
    return -np.log2(tot / N)
pool = [rand_stab() for _ in range(64)]
print("n = %d, saturation (Haar) ~ %.2f" % (n, np.log2((N + 3) / 4)))
for k in [1, 2, 3, 4, 6, 8]:
    vals = []
    for rep in range(4):
        idx = rng.choice(len(pool), k, replace=False)
        psi = sum(np.exp(2j*np.pi*rng.random()) * pool[i] for i in idx)
        vals.append(M2(psi / np.linalg.norm(psi)))
    print("k=%2d  M2 = %.3f +- %.3f   3*log2(k) = %.3f" % (k, np.mean(vals), np.std(vals), 3*np.log2(k)))
