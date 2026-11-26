import numpy as np, itertools as it
from functools import reduce
np.random.seed(1)
n = 3; N = 2**n
I2 = np.eye(2); X = np.array([[0,1],[1,0]]); Z = np.diag([1,-1]); Y = 1j*X@Z
Hd = np.array([[1,1],[1,-1]])/np.sqrt(2); Sg = np.diag([1,1j])
kron = lambda ms: reduce(np.kron, ms)
def on(g, q): return kron([g if i == q else I2 for i in range(n)])
def cnot(a, b):
    U = np.zeros((N, N))
    for x in range(N):
        bits = [(x >> (n-1-i)) & 1 for i in range(n)]
        if bits[a]: bits[b] ^= 1
        U[int("".join(map(str, bits)), 2), x] = 1
    return U
gens = [on(Hd, q) for q in range(n)] + [on(Sg, q) for q in range(n)] + [cnot(a, b) for a in range(n) for b in range(n) if a != b]
def key(v):
    k = np.argmax(np.abs(v) > 1e-9); v = v * np.conj(v[k]) / abs(v[k])
    return tuple(np.round(v, 6).view(float))
start = np.zeros(N, complex); start[0] = 1
seen = {key(start): start}; frontier = [start]
while frontier:
    new = []
    for v in frontier:
        for g in gens:
            w = g @ v; kk = key(w)
            if kk not in seen: seen[kk] = w; new.append(w)
    frontier = new
STAB = np.array(list(seen.values())); print("n=3 stabilizer states:", len(STAB))

paulis = [kron(p) for p in it.product([I2, X, Y, Z], repeat=n)]
def M2(psi):
    ev = np.array([np.real(np.conj(psi) @ P @ psi) for P in paulis])
    return -np.log2(np.sum(ev**4) / N)
def F(psi): return np.max(np.abs(STAB.conj() @ psi)**2)
def U3_8(psi):                       # ||2^{n/2} psi||_{U^3}^8
    f = np.sqrt(N) * psi; tot = 0
    for x, h1, h2, h3 in it.product(range(N), repeat=4):
        p = 1
        for w in it.product([0, 1], repeat=3):
            y = x ^ (h1*w[0]) ^ (h2*w[1]) ^ (h3*w[2])
            p *= f[y] if sum(w) % 2 == 0 else np.conj(f[y])
        tot += p
    return np.real(tot) / N**4
worst = 0
for trial in range(200):
    psi = np.random.randn(N) + 1j*np.random.randn(N); psi /= np.linalg.norm(psi)
    if trial < 3: print("  2^-M2 = %.6f   U3^8 = %.6f" % (2**-M2(psi), U3_8(psi)))
    worst = max(worst, M2(psi) / (-np.log2(F(psi))))
print("max over 200 random states of M2 / log2(1/F) = %.3f  (claim: <= 4)" % worst)
c, s = np.cos(np.pi/8), np.sin(np.pi/8); h = kron([np.array([c, s])]*n)
print("H^3: M2 = %.4f (3*log2(4/3) = %.4f), F = %.4f (cos^6 = %.4f)" % (M2(h), 3*np.log2(4/3), F(h), c**6))
