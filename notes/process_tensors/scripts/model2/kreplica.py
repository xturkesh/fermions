# Exact k-replica recursion (any D, Haar average with S_k Weingarten) for the bottleneck model,
# with depolarizing noise p on E applied copy-wise. C_k(t) = Tr[(A(t)B)^k]/(d^2 D).
import numpy as np, itertools as it
from functools import reduce
def cycles(p):
    seen, c = set(), 0
    for i in range(len(p)):
        if i not in seen:
            c += 1; j = i
            while j not in seen: seen.add(j); j = p[j]
    return c
comp = lambda a, b: tuple(a[b[i]] for i in range(len(a)))          # (a∘b)(i)=a(b(i))
inv = lambda a: tuple(np.argsort(a))
def perm_op(p, dim):                                                   # P_p |i_1..i_k> = |i_{p^-1(1)} ...>
    k = len(p); M = np.zeros((dim**k, dim**k))
    for idx in it.product(range(dim), repeat=k):
        out = tuple(idx[inv(p)[j]] for j in range(k))
        M[np.ravel_multi_index(out, [dim]*k), np.ravel_multi_index(idx, [dim]*k)] = 1
    return M
def ptrace_perm(tau, R):                                               # Tr_R(P_tau) = D^n * P_tau' (tau' on [k], identity on R)
    tau = list(tau); n = 0
    for j in R:
        if tau[j] == j: n += 1
        else:
            pre = tau.index(j); tau[pre] = tau[j]; tau[j] = j
    return tuple(tau), n
def TE(tau, sig, D, p):                                                # Tr[ N^{dag ⊗k}(P_tau) P_sig^dag ], N depolarizing(p)
    k = len(tau); tot = 0.0
    for r in range(k+1):
        for R in it.combinations(range(k), r):
            tp, n = ptrace_perm(tau, R)
            tot += (1-p)**(k-r) * p**r * D**(n - r) * D**cycles(comp(tp, inv(sig)))
    return tot
def otoc_k(U, A, B, t, D, k, p=0.0):
    d = A.shape[0]; q = d*D; perms = list(it.permutations(range(k)))
    G = np.array([[float(q)**cycles(comp(inv(s), tt)) for tt in perms] for s in perms]); Gi = np.linalg.inv(G)
    PC = {s: perm_op(s, d) for s in perms}
    # reorder U^{⊗k} from (a1 c1 a2 c2 ...) to (a1..ak, c1..ck)
    Uk = reduce(np.kron, [U]*k).reshape([d]*(4*k))
    order = [2*j for j in range(k)] + [2*j+1 for j in range(k)]
    Uk = Uk.transpose(order + [2*k + o for o in order]).reshape(d**(2*k), d**(2*k))
    Ak = reduce(np.kron, [A]*k); dk = d**k
    Xs = {s: (Ak.copy() if s == tuple(range(k)) else np.zeros((dk, dk), complex)) for s in perms}
    TEm = {(a, b): TE(a, b, D, p) for a in perms for b in perms}
    for step in range(t):
        Ys = {s: Uk.conj().T @ np.kron(Xs[s], PC[s]) @ Uk for s in perms}
        if step == t-1: break
        red = {(s, sg): np.einsum('aibj,ji->ab', Ys[s].reshape(dk, dk, dk, dk), PC[sg].T.conj()) for s in perms for sg in perms}
        new = {}
        for ti, tp in enumerate(perms):
            acc = np.zeros((dk, dk), complex)
            for si, sg in enumerate(perms):
                b = sum(TEm[(s, sg)] * red[(s, sg)] for s in perms)     # Tr_{CE}[Y P_sg^dag] (A-operator)
                acc += Gi[ti, si] * b
            new[tp] = acc
        Xs = new
    gam = tuple((j+1) % k for j in range(k))                              # cyclic: P|i1..ik> = |i2..ik i1>
    Pg_AC = None
    Bk = reduce(np.kron, [B]*k)
    PA, PCg = perm_op(inv(gam), d), perm_op(inv(gam), d)
    tot = 0.0
    for s in perms:
        T = np.kron(Bk @ PA, PCg)                                         # (B^{⊗k} ⊗ I) * P_gamma on (A..., C...)
        tot += np.trace(Ys[s] @ np.kron(Bk, np.eye(dk)) @ np.kron(PA, PCg)) * D**cycles(comp(s, inv(gam)))
    return np.real(tot)/(d*d*D)
if __name__ == "__main__":
    from bottleneck import gate, recursion, brute, X, Z
    from limit import otoc_inf
    U = gate(0.5, 0.3); I = np.eye(2)
    print("k=2 general code vs dedicated recursion:")
    for (t, D, p) in [(3, 2, 0.0), (4, 5, 0.3), (5, 1000, 0.2)]:
        print("  t=%d D=%d p=%.1f  general=%.6f  dedicated=%.6f" % (t, D, p, otoc_k(U, Z, Z, t, D, 2, p), recursion(U, Z, Z, t, D, pE=p)))
