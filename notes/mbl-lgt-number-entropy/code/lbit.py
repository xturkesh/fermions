"""Common tools: fixed-N sector basis, U(1) two-site gates, brickwork circuits
aligned with the half-chain cut, l-bit energies D_c, number distribution of N_A.

Conventions: open chain, sites 0..L-1, bit i of the integer state = n_i.
A = {0,...,l-1}, l = L//2, cut bond b = (l-1, l).
Brickwork: half-layer 'a' contains the cut bond and all bonds (l-1+2k, l+2k);
half-layer 'b' the complementary bonds.  A 'full layer' = a then b, so each
bond is acted on d_u times by a circuit with d_u full layers.  Gates are
indexed by (layer, position relative to cut), so circuits for different L can
share exactly the same gates near the cut.
"""
import numpy as np
from itertools import combinations


def sector_basis(L, N):
    states = []
    for occ in combinations(range(L), N):
        s = 0
        for i in occ:
            s |= 1 << i
        states.append(s)
    return np.array(sorted(states), dtype=np.int64)


def popcount(a):
    a = a.copy()
    c = np.zeros_like(a)
    while np.any(a):
        c += a & 1
        a >>= 1
    return c


class Sector:
    def __init__(self, L, N=None):
        self.L = L
        self.N = L // 2 if N is None else N
        self.states = sector_basis(L, self.N)
        self.D = len(self.states)
        self.l = L // 2
        maskA = (1 << self.l) - 1
        self.NA = popcount(self.states & maskA)
        self.bond = {}
        for i in range(L - 1):
            mi, mj = 1 << i, 1 << (i + 1)
            bi = (self.states & mi) != 0
            bj = (self.states & mj) != 0
            i00 = np.where(~bi & ~bj)[0]
            i11 = np.where(bi & bj)[0]
            i01 = np.where(~bi & bj)[0]          # n_i=0, n_{i+1}=1
            partner = self.states[i01] ^ (mi | mj)
            i10 = np.searchsorted(self.states, partner)
            self.bond[i] = (i00, i11, i01, i10)

    def index(self, s):
        return int(np.searchsorted(self.states, s))

    def neel(self, parity):
        s = 0
        for i in range(self.L):
            if i % 2 == parity:
                s |= 1 << i
        return s


def random_u1_gate(rng, f):
    """g = exp(-i f M), M Hermitian, U(1) symmetric: blocks (00),(01,10),(11).
    Returns (g00, g11, G2x2)."""
    a, e = rng.normal(size=2)
    b, d = rng.normal(size=2)
    c = (rng.normal() + 1j*rng.normal())/np.sqrt(2)
    M2 = np.array([[b, c], [np.conj(c), d]])
    w, v = np.linalg.eigh(M2)
    G = (v*np.exp(-1j*f*w)) @ v.conj().T
    return np.exp(-1j*f*a), np.exp(-1j*f*e), G


def gate_dist_from_identity(g):
    g00, g11, G = g
    return max(abs(g00 - 1), abs(g11 - 1), np.linalg.norm(G - np.eye(2), 2))


def make_circuit_gates(rng, d_u, f, rel_range=64):
    """Gates for all relative bond positions r in [-rel_range, rel_range]
    (bond (l-1+r, l+r)).  Layer 2k   ('a'): r even.  Layer 2k+1 ('b'): r odd."""
    gates = {}
    for layer in range(2*d_u):
        par = layer % 2
        for r in range(-rel_range, rel_range + 1):
            if (r - par) % 2 == 0:
                gates[(layer, r)] = random_u1_gate(rng, f)
    return gates


def circuit_list(sec, gates, d_u):
    """Ordered list of (bond_site_i, gate) for U = G_T ... G_1 (applied first to last)."""
    out = []
    l = sec.l
    for layer in range(2*d_u):
        par = layer % 2
        for i in range(sec.L - 1):
            r = i - (l - 1)
            if (r - par) % 2 == 0:
                out.append((i, gates[(layer, r)]))
    return out


def apply_gate(sec, psi, i, g, dagger=False):
    g00, g11, G = g
    if dagger:
        g00, g11, G = np.conj(g00), np.conj(g11), G.conj().T
    i00, i11, i01, i10 = sec.bond[i]
    psi[i00] *= g00
    psi[i11] *= g11
    a = psi[i01].copy()
    b = psi[i10].copy()
    # basis order in G: (01, 10)
    psi[i01] = G[0, 0]*a + G[0, 1]*b
    psi[i10] = G[1, 0]*a + G[1, 1]*b
    return psi


def apply_U(sec, psi, clist):
    for i, g in clist:
        apply_gate(sec, psi, i, g)
    return psi


def apply_Udag(sec, psi, clist):
    for i, g in reversed(clist):
        apply_gate(sec, psi, i, g, dagger=True)
    return psi


def lbit_energies_full(rng, L, W=1.0, J0=1.0, kappa=0.5, halfpow=True):
    """D_c = sum_X J_X prod_{i in X} c_i over ALL subsets X (|X|>=1):
    J_{i} = h_i ~ U[-W,W];  |X|>=2: J_X = J0 * u_X * exp(-(diam X - 1)/kappa) * 2^{-(|X|-2)},
    u_X ~ U[-1,1].  Returns the full 2^L array (zeta transform)."""
    n = 1 << L
    masks = np.arange(n, dtype=np.int64)
    pc = popcount(masks)
    low = np.zeros(n, dtype=np.int64)
    high = np.zeros(n, dtype=np.int64)
    for i in range(L):
        bit = (masks >> i) & 1
        high = np.where(bit == 1, i, high)
    for i in reversed(range(L)):
        bit = (masks >> i) & 1
        low = np.where(bit == 1, i, low)
    diam = high - low
    Jc = np.zeros(n)
    single = pc == 1
    Jc[single] = rng.uniform(-W, W, size=single.sum())
    multi = pc >= 2
    fac = np.exp(-(diam[multi] - 1)/kappa)
    if halfpow:
        fac = fac * 0.5**(pc[multi] - 2)
    Jc[multi] = J0*rng.uniform(-1, 1, size=multi.sum())*fac
    # zeta transform: D[c] = sum_{X subset c} Jc[X]
    Dfull = Jc.copy()
    for i in range(L):
        Dfull = Dfull.reshape(-1, 2, 1 << i)
        Dfull[:, 1, :] += Dfull[:, 0, :]
        Dfull = Dfull.reshape(-1)
    return Dfull


def number_distribution(sec, Psi):
    """Psi: (D,) or (D, m) -> p: (l+1,) or (l+1, m)."""
    P = np.abs(Psi)**2
    if P.ndim == 1:
        return np.bincount(sec.NA, weights=P, minlength=sec.l + 1)
    out = np.zeros((sec.l + 1, P.shape[1]))
    np.add.at(out, sec.NA, P)
    return out


def shannon(p, axis=0):
    p = np.clip(p, 1e-300, None)
    return -(p*np.log(p)).sum(axis=axis)


def variance(p, axis=0):
    n = np.arange(p.shape[0]).reshape((-1,) + (1,)*(p.ndim - 1))
    m1 = (p*n).sum(axis=0)
    m2 = (p*n**2).sum(axis=0)
    return m2 - m1**2
