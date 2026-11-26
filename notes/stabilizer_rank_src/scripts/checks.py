import numpy as np, itertools as it
from functools import reduce
np.random.seed(0)
c, s = np.cos(np.pi/8), np.sin(np.pi/8)
H1 = np.array([c, s])                      # H-type magic state (Clifford-equiv. to T|+>)
D = [np.array(v, complex)/np.linalg.norm(v) for v in
     ([1,0],[0,1],[1,1],[1,-1],[1,1j],[1,-1j])]   # 6 single-qubit stabilizer states
kron = lambda vs: reduce(np.kron, vs)

def in_span(vecs, t, tol=1e-9):
    A = np.array(vecs).T
    x, *_ = np.linalg.lstsq(A, t, rcond=None)
    return np.linalg.norm(A @ x - t) < tol

# (a) explicit product decompositions
t2 = kron([H1, H1]); print("H^2 in span{|00>,|11>,|++>}:", in_span([kron([D[0],D[0]]), kron([D[1],D[1]]), kron([D[2],D[2]])], t2))
t5 = kron([H1]*5);  print("H^5 in span{d^5, d in D}:   ", in_span([kron([d]*5) for d in D], t5))

# (b) exhaustive product-stabilizer rank for n=2,3
def prod_rank(n, kmax):
    dic = np.array([kron(p) for p in it.product(D, repeat=n)])
    t = kron([H1]*n)
    for k in range(1, kmax+1):
        for S in it.combinations(range(len(dic)), k):
            if in_span(dic[list(S)], t): return k
    return None
print("prod-rank H^2 =", prod_rank(2, 4))

def prod_rank_batched(n, k):
    dic = np.array([kron(p) for p in it.product(D, repeat=n)]); t = kron([H1]*n)
    combos = np.array(list(it.combinations(range(len(dic)), k)))
    found = False
    for chunk in np.array_split(combos, 40):
        A = np.transpose(dic[chunk], (0, 2, 1))            # (B, 2^n, k)
        Q, _ = np.linalg.qr(A)
        r = t[None, :, None] - Q @ (np.conj(np.transpose(Q, (0, 2, 1))) @ t[None, :, None])
        rk = np.linalg.matrix_rank(A, tol=1e-9)
        ok = (np.linalg.norm(r[:, :, 0], axis=1) < 1e-9) & (rk == k)
        if ok.any(): found = True; break
    return found
print("prod-rank H^3 <= 3 ?", prod_rank_batched(3, 3))
