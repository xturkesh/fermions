"""Noisy Haar collision model: S (dim d) + E (L qudits of dim q, D=q^L).
Each step: fresh S input (Choi leg), Haar U on S⊗E, then noise channel (Kraus list) on every env qudit.
Annealed k-replica quantities via S_k Weingarten transfer in 'moment space'
  m_alpha(Y) = Tr(P_alpha^{-1} Y)  (Y = k-replica operator on E).
Per step with S-leg contractions (pin, pout):
  coefficients y_beta = sum_alpha Wg(alpha^{-1} beta, dD) d^{#(alpha^{-1} pin^{-1})} d^{#(pout beta)} m_alpha
  moments after noise m'_gamma = sum_beta Tr(P_gamma^{-1} N^{(x)k}(P_beta)) y_beta   (per-qudit factor ^L)
normalized by d^k per step (Choi normalized to trace 1)."""
import numpy as np, itertools as it
from functools import reduce

def perms(k): return list(it.permutations(range(k)))
def compose(s, t): return tuple(s[t[i]] for i in range(len(s)))       # (s t)(i) = s(t(i))
def inverse(s):
    inv = [0]*len(s)
    for i, j in enumerate(s): inv[j] = i
    return tuple(inv)
def ncyc(s):
    seen, c = set(), 0
    for i in range(len(s)):
        if i not in seen:
            c += 1; j = i
            while j not in seen: seen.add(j); j = s[j]
    return c
def cyc(k): return tuple((i+1) % k for i in range(k))                  # i -> i+1
def transp(k, a, b):
    s = list(range(k)); s[a], s[b] = b, a; return tuple(s)
IDP = lambda k: tuple(range(k))

def perm_op(s, q):
    """P_s |x_1..x_k> = |x_{s^{-1}(1)} .. x_{s^{-1}(k)}>, so that P_s P_t = P_{st}."""
    k = len(s); dim = q**k; P = np.zeros((dim, dim))
    sinv = inverse(s)
    for x in it.product(range(q), repeat=k):
        y = tuple(x[sinv[i]] for i in range(k))
        P[np.ravel_multi_index(y, [q]*k), np.ravel_multi_index(x, [q]*k)] = 1
    return P

def wg_matrix(k, N):
    P = perms(k)
    Q = np.array([[float(N)**ncyc(compose(inverse(a), b)) for b in P] for a in P])
    return np.linalg.inv(Q)                                              # W[a,b] = Wg(a^{-1} b)

def qudit_noise_matrix(k, q, kraus):
    """g[gamma,beta] = Tr(P_gamma^{-1} N^{(x)k}(P_beta)) for one qudit."""
    P = perms(k); ops = [perm_op(s, q) for s in P]
    def Nk(X):
        out = np.zeros_like(X, dtype=complex)
        for ks in it.product(kraus, repeat=k):
            K = reduce(np.kron, ks); out += K @ X @ K.conj().T
        return out
    NP = [Nk(o) for o in ops]
    return np.real(np.array([[np.trace(ops[gi].T @ NP[bi]) for bi in range(len(P))] for gi in range(len(P))]))

def qudit_initial_moments(k, q, rho):
    P = perms(k); rk = reduce(np.kron, [rho]*k)
    return np.real(np.array([np.trace(perm_op(s, q).T @ rk) for s in P]))

class Model:
    def __init__(self, k, d, q, L, kraus=None, rho_q=None):
        self.k, self.d, self.q, self.L = k, d, q, L
        self.D = q**L; self.P = perms(k); self.idx = {s: i for i, s in enumerate(self.P)}
        self.W = wg_matrix(k, d*self.D)
        kraus = kraus if kraus is not None else [np.eye(q)]
        self.G = qudit_noise_matrix(k, q, kraus)**L                      # product over qudits
        rho_q = rho_q if rho_q is not None else np.eye(q)/q
        self.m0 = qudit_initial_moments(k, q, rho_q)**L
        self.cache = {}
    def H(self, pin, pout):
        key = (pin, pout)
        if key not in self.cache:
            k, d, P = self.k, self.d, self.P; n = len(P); Hm = np.zeros((n, n))
            pin_inv = inverse(pin)
            for bi, b in enumerate(P):
                for ai, a in enumerate(P):
                    Hm[bi, ai] = self.W[ai, bi] * d**ncyc(compose(inverse(a), pin_inv)) * d**ncyc(compose(pout, b))
            self.cache[key] = (self.G @ Hm) / d**k
        return self.cache[key]
    def run(self, pattern, m0=None):
        v = self.m0.copy() if m0 is None else m0
        for (pin, pout) in pattern: v = self.H(pin, pout) @ v
        return v
    def value(self, pattern): return self.run(pattern)[self.idx[IDP(self.k)]]

# ---------------- Monte Carlo: exact per-sample process tensor (density-matrix simulation) ----------------
def haar(N, rng):
    Z = (rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))/np.sqrt(2)
    Qm, R = np.linalg.qr(Z); return Qm * (np.diag(R)/np.abs(np.diag(R)))

def env_channel(rho_full_tensor, kraus, q, L, axes_ket, axes_bra):
    """apply single-qudit channel to every env qudit; tensor has separate axes for each env qudit (ket & bra)."""
    T = rho_full_tensor
    for l in range(L):
        ak, ab = axes_ket[l], axes_bra[l]
        out = 0
        for K in kraus:
            t = np.tensordot(K, T, axes=(1, ak)); t = np.moveaxis(t, 0, ak)
            t = np.tensordot(t, K.conj(), axes=(ab, 1)); t = np.moveaxis(t, -1, ab)
            out = out + t
        T = out
    return T

def sample_upsilon(d, q, L, n, kraus, rng, rho_q=None):
    """normalized Choi of the process tensor, registers ordered (R1 O1 ... Rn On)."""
    D = q**L; rho_q = np.eye(q)/q if rho_q is None else rho_q
    rho = reduce(np.kron, [rho_q]*L).astype(complex)                     # (E, E)
    phi = np.zeros(d*d, complex); phi[[i*d+i for i in range(d)]] = 1/np.sqrt(d)
    Phi = np.outer(phi, phi.conj())
    for j in range(n):
        dx = rho.shape[0]//D
        r = rho.reshape(dx, D, dx, D)
        full = np.einsum('aebf,rs->arebsf', r, Phi).reshape(dx, d, d, D, dx, d, d, D)
        U = haar(d*D, rng).reshape(d, D, d, D)
        full = np.einsum('tfse,xrseyRSE,TFSE->xrtfyRTF', U, full, U.conj(), optimize=True)
        if kraus is not None and len(kraus) > 1:
            T = full.reshape([dx, d, d] + [q]*L + [dx, d, d] + [q]*L)
            T = env_channel(T, kraus, q, L, list(range(3, 3+L)), list(range(6+L, 6+2*L)))
            full = T.reshape(dx, d, d, D, dx, d, d, D)
        rho = full.reshape(dx*d*d*D, dx*d*d*D)
    N = rho.shape[0]//D
    return np.einsum('aebe->ab', rho.reshape(N, D, N, D))

def marginals(Ups, d, n1, n2):
    P, F = d**(2*n1), d**(2*n2); U4 = Ups.reshape(P, F, P, F)
    return np.einsum('afbf->ab', U4), np.einsum('afag->fg', U4)

# ---------------- noise channels (qubit / qudit) ----------------
def depolarizing(q, p):
    # N(rho) = (1-p) rho + p Tr(rho) I/q via Weyl twirl
    w = np.exp(2j*np.pi/q); X = np.roll(np.eye(q), 1, axis=0); Z = np.diag(w**np.arange(q))
    weyl = [np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b) for a in range(q) for b in range(q)]
    ks = [np.sqrt(1 - p + p/q**2)*np.eye(q)] + [np.sqrt(p)/q*W for W in weyl[1:]]
    return ks
def dephasing(p):   # qubit: (1-p) rho + p Z rho Z
    return [np.sqrt(1-p)*np.eye(2), np.sqrt(p)*np.diag([1., -1.])]
def amp_damp(g):
    return [np.array([[1, 0], [0, np.sqrt(1-g)]]), np.array([[0, np.sqrt(g)], [0, 0]])]
