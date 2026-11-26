"""Exact (rational) version of the S_k transfer, to avoid cancellation errors at large D / k=3,4.
Noise: depolarizing with rational p on each env qudit (or none)."""
import sys, numpy as np, itertools as it
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from common import perms, compose, inverse, ncyc, cyc, transp, IDP, perm_op
def inv_exact(A):
    n = len(A); M = [list(row) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]; M[c] = [x/pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f*y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]
def depol_qudit_matrix(k, q, p):
    """g[gamma][beta] = Tr(P_gamma^{-1} N^{(x)k}(P_beta)), N = (1-p) id + p R, R(X)=Tr(X) I/q, exact."""
    P = perms(k); ops = [perm_op(s, q).reshape([q]*(2*k)) for s in P]
    g = [[Fr(0)]*len(P) for _ in P]
    for A in it.chain.from_iterable(it.combinations(range(k), r) for r in range(k+1)):
        w = (1-p)**(k-len(A)) * p**len(A)
        for bi, Tb in enumerate(ops):
            Z = Tb.copy()
            for i in A:                                   # replace factor i by Tr(.) I/q
                tr = np.trace(Z, axis1=i, axis2=k+i)      # removes axes i, k+i
                Z = np.expand_dims(np.expand_dims(tr, i), k+i) * (np.eye(q)/q).reshape([q if a in (i, k+i) else 1 for a in range(2*k)])
            for gi, Tg in enumerate(ops):
                val = np.sum(Tg*Z)                         # Tr(P_g^T Z) = sum_xy P_g[x,y] Z[x,y]
                g[gi][bi] += w*Fr(val).limit_denominator(q**(2*k))
    return g
class ExactModel:
    def __init__(self, k, d, q, L, p=None):
        self.k, self.d, self.q, self.L = k, d, q, L; self.D = q**L; self.P = perms(k)
        N = d*self.D
        Q = [[Fr(N)**ncyc(compose(inverse(a), b)) for b in self.P] for a in self.P]
        self.W = inv_exact(Q)
        if p is None: g = [[Fr(q)**ncyc(compose(inverse(a), b)) for b in self.P] for a in self.P]
        else: g = depol_qudit_matrix(k, q, Fr(p))
        self.G = [[x**L for x in row] for row in g]
        self.m0 = [Fr(q)**(L*(ncyc(s) - k)) for s in self.P]
        self.cache = {}
    def H(self, pin, pout):
        if (pin, pout) not in self.cache:
            k, d = self.k, self.d; n = len(self.P); pinv = inverse(pin)
            Hm = [[self.W[ai][bi]*Fr(d)**ncyc(compose(inverse(a), pinv))*Fr(d)**ncyc(compose(pout, b))
                   for ai, a in enumerate(self.P)] for bi, b in enumerate(self.P)]
            GH = [[sum(self.G[gi][bi]*Hm[bi][ai] for bi in range(n))/Fr(d)**k for ai in range(n)] for gi in range(n)]
            self.cache[(pin, pout)] = GH
        return self.cache[(pin, pout)]
    def value(self, pattern):
        v = list(self.m0); n = len(v)
        for (pin, pout) in pattern:
            M = self.H(pin, pout); v = [sum(M[i][j]*v[j] for j in range(n)) for i in range(n)]
        return v[self.P.index(IDP(self.k))]
if __name__ == "__main__":
    from common import Model, depolarizing
    m = ExactModel(2, 2, 2, 1); F = transp(2, 0, 1); print("exact P(n=1..3):", [str(m.value([(F, F)]*n)) for n in [1, 2, 3]], "(0.4, 0.28, 0.256)")
    for k in [3, 4]:
        me = ExactModel(k, 2, 2, 1, Fr(3, 10)); mf = Model(k, 2, 2, 1, depolarizing(2, 0.3)); c = cyc(k)
        print("k=%d noisy p=0.3 n=3: exact %.10f  float %.10f" % (k, float(me.value([(c, c)]*3)), mf.value([(c, c)]*3)))
