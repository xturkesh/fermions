"""Exact small-L checks of T1 (symmetric circuit U), T1' (non-symmetric U with
U(1)-invariant LIOMs) and Proposition T2.1 (dephased number distribution is
independent of D and, for L >= 8 d_u, of L).
"""
import numpy as np
from lbit import *

rng = np.random.default_rng(7)

# ---------------------------------------------------------------- T1 (a)
def dense_U(sec, clist):
    M = np.eye(sec.D, dtype=complex)
    for i, g in clist:
        apply_gate(sec, M, i, g)
    return M

print("=== T1(a): symmetric brickwork, L=8, all N sectors ===")
L, d_u, f = 8, 3, 0.6
gates = make_circuit_gates(rng, d_u, f, rel_range=12)
cut_gates = [gates[(lay, 0)] for lay in range(0, 2*d_u, 2)]
delta_circ = sum(min(1.0, gate_dist_from_identity(g)) for g in cut_gates)
deltaA = 0.0
for N in range(L + 1):
    sec = Sector(L, N)
    U = dense_U(sec, circuit_list(sec, gates, d_u))
    NA = np.diag(sec.NA.astype(float))
    Delta = U.conj().T @ NA @ U - NA
    deltaA = max(deltaA, np.abs(np.linalg.eigvalsh(Delta)).max())
print(f"delta_A = ||U^+ N_A U - N_A|| = {deltaA:.4f};  circuit bound sum_cut min(1,||g-1||) = {delta_circ:.4f};  #cut gates d_u = {d_u}")

sec = Sector(L)
clist = circuit_list(sec, gates, d_u)
U = dense_U(sec, clist)
Dfull = lbit_energies_full(rng, L)
Dc = Dfull[sec.states]
NA = np.diag(sec.NA.astype(float))
worst = 0.0
for t in 10.0**rng.uniform(-1, 8, size=40):
    Ut = U @ np.diag(np.exp(-1j*Dc*t)) @ U.conj().T
    NAt = Ut.conj().T @ NA @ Ut
    worst = max(worst, np.abs(np.linalg.eigvalsh(NAt - NA)).max())
print(f"max_t ||N_A(t)-N_A|| over 40 random t (half filling) = {worst:.4f}  <=  2 delta_A = {2*deltaA:.4f}")

# support bound for product states
psi0 = np.zeros(sec.D, complex); psi0[sec.index(sec.neel(0))] = 1
n0 = sec.NA[sec.index(sec.neel(0))]
maxdev = 0
for t in 10.0**rng.uniform(-1, 8, size=40):
    psit = U @ (np.exp(-1j*Dc*t) * (U.conj().T @ psi0))
    p = number_distribution(sec, psit)
    supp = np.where(p > 1e-14)[0]
    maxdev = max(maxdev, np.abs(supp - n0).max())
print(f"support: max |n - n0| with p_t(n) > 1e-14 = {maxdev}  <= 2 d_u = {2*d_u}  (trivial bound {min(n0, L//2 - n0) if False else L//2})")

# ---------------------------------------------------------------- T1' non-symmetric U
print("\n=== T1': U = V C with C a classical reversible (non-U(1)) circuit, L=8 full space ===")
L = 8
dim = 1 << L
nvec = np.array([[(s >> i) & 1 for i in range(L)] for s in range(dim)], dtype=float)

def full_gate_matrix(i, g):
    g00, g11, G = g
    M = np.zeros((dim, dim), complex)
    for s in range(dim):
        bi, bj = (s >> i) & 1, (s >> (i+1)) & 1
        if bi == 0 and bj == 0: M[s, s] = g00
        elif bi == 1 and bj == 1: M[s, s] = g11
        elif bi == 0 and bj == 1:
            s2 = s ^ ((1 << i) | (1 << (i+1)))
            M[s, s] += G[0, 0]; M[s2, s] += G[1, 0]
        else:
            s2 = s ^ ((1 << i) | (1 << (i+1)))
            M[s, s] += G[1, 1]; M[s2, s] += G[0, 1]
    return M

V = np.eye(dim, dtype=complex)
d_u = 2
gatesV = make_circuit_gates(rng, d_u, 0.5, rel_range=12)
for layer in range(2*d_u):
    for i in range(L - 1):
        r = i - (L//2 - 1)
        if (r - layer % 2) % 2 == 0:
            V = full_gate_matrix(i, gatesV[(layer, r)]) @ V
# classical reversible circuit: random Toffoli / CNOT / X gates on nearby sites
perm = np.arange(dim)
for k in range(12):
    typ = rng.integers(3); i = rng.integers(0, L - 2)
    new = perm.copy()
    for s in range(dim):
        if typ == 0: t = s ^ (1 << i)                                        # X
        elif typ == 1: t = s ^ (((s >> i) & 1) << (i+1))                      # CNOT i->i+1
        else: t = s ^ ((((s >> i) & 1) & ((s >> (i+1)) & 1)) << (i+2))        # Toffoli
        new[s] = t
    perm = new[perm]
C = np.zeros((dim, dim)); C[perm, np.arange(dim)] = 1
U = V @ C
Dfull = lbit_energies_full(rng, L)
H = U @ np.diag(Dfull) @ U.conj().T
Ntot = np.diag(nvec.sum(1))
print("||[H,N]|| =", np.abs(H @ Ntot - Ntot @ H).max(), "  ||[U,N]|| =", np.abs(U @ Ntot - Ntot @ U).max(), "(U not symmetric)")
l = L//2
NAop = np.diag(nvec[:, :l].sum(1))
m = [U.conj().T @ np.diag(nvec[:, x]) @ U for x in range(L)]
g = [mx - np.diag(np.diag(mx)) for mx in m]
print("||sum_x g_x|| =", np.abs(sum(g)).max(), "(=0 since U^+ N U is diagonal)")
BA = sum(g[:l])
nBA = np.abs(np.linalg.eigvalsh(BA)).max()
worst = 0.0
for t in 10.0**rng.uniform(-1, 8, size=30):
    Ut = U @ np.diag(np.exp(-1j*Dfull*t)) @ U.conj().T
    worst = max(worst, np.abs(np.linalg.eigvalsh(Ut.conj().T @ NAop @ Ut - NAop)).max())

# eta(j) = sup_x || m_x - Pi_{B_j(x)} m_x ||  (normalized partial trace)
def ptrace_keep(O, keep):
    """Pi_Y(O): normalized partial trace over complement of Y, tensored with identity."""
    O = O.reshape([2]*L*2)   # axes: out bits (L-1..0 order?) -> use explicit
    return None

def Pi(O, keep):
    # bit i of index s <-> site i ; reshape with axis order (bit L-1,...,bit 0)
    T = O.reshape([2]*(2*L))
    out_ax = {i: L-1-i for i in range(L)}
    in_ax = {i: 2*L-1-i for i in range(L)}
    Tm = T.copy()
    for i in range(L):
        if i not in keep:
            Tm = (np.take(Tm, 0, axis=out_ax[i]).take(0, axis=in_ax[i]-1) +
                  np.take(Tm, 1, axis=out_ax[i]).take(1, axis=in_ax[i]-1))/2
            Tm = np.expand_dims(np.expand_dims(Tm, out_ax[i]), in_ax[i])
            eye = np.eye(2).reshape([2 if a in (out_ax[i], in_ax[i]) else 1 for a in range(2*L)])
            Tm = Tm*eye
    return Tm.reshape(dim, dim)

eta = []
for j in range(L):
    e = 0.0
    for x in range(L):
        keep = set(range(max(0, x-j), min(L, x+j+1)))
        e = max(e, np.abs(np.linalg.eigvalsh(m[x] - Pi(m[x], keep))).max())
    eta.append(e)
bound = 24*1*sum((j+1)*eta[j] for j in range(L))
print(f"||B_A|| = {nBA:.4f};  max_t ||N_A(t)-N_A|| = {worst:.4f} <= 2||B_A|| = {2*nBA:.4f};  "
      f"crude bound 24|dA| sum (j+1) eta(j) = {bound:.2f}")
print("eta(j) =", np.round(eta, 4))

# ---------------------------------------------------------------- Prop T2.1
print("\n=== Prop T2.1: dephased distribution p_bar is independent of D and (L >= 8 d_u) of L ===")
def pbar_exact(sec, clist, psi):
    U = dense_U(sec, clist)
    phi = U.conj().T @ psi
    q = number_distribution(sec, U)           # (l+1, D): distribution of N_A in U|c>
    return q @ (np.abs(phi)**2)

for d_u in [1, 2]:
    gates = make_circuit_gates(rng, d_u, 0.8, rel_range=20)
    res = {}
    for L in [8, 10, 12, 14, 16] if d_u == 1 else [12, 14, 16]:
        if L == 16 and d_u == 1:
            pass
        sec = Sector(L)
        clist = circuit_list(sec, gates, d_u)
        # Neel with site l (first site right of cut) occupied, relative to cut
        par = sec.l % 2
        psi = np.zeros(sec.D, complex); psi[sec.index(sec.neel(par))] = 1
        n0 = sec.NA[sec.index(sec.neel(par))]
        if sec.D > 5000:
            continue
        pb = pbar_exact(sec, clist, psi)
        # relative to n0
        rel = {n - n0: pb[n] for n in range(len(pb)) if pb[n] > 1e-15}
        res[L] = rel
        print(f"d_u={d_u} L={L}: p_bar(n-n0) =", {k: round(v, 10) for k, v in sorted(rel.items())})
    # D-independence is automatic in formula; check against long-time average with two different D
sec = Sector(10)
gates = make_circuit_gates(rng, 1, 0.8, rel_range=20)
clist = circuit_list(sec, gates, 1)
par = sec.l % 2
psi = np.zeros(sec.D, complex); psi[sec.index(sec.neel(par))] = 1
U = dense_U(sec, clist)
phi = U.conj().T @ psi
pb = pbar_exact(sec, clist, psi)
for trial in range(2):
    Dc = lbit_energies_full(rng, 10, W=1.0 + trial, kappa=0.5 + trial)[sec.states]
    ts = rng.uniform(0, 1e7, size=4000)
    acc = np.zeros(sec.l + 1)
    for t0 in range(0, len(ts), 500):
        tt = ts[t0:t0+500]
        Psi = U @ (np.exp(-1j*np.outer(Dc, tt)) * phi[:, None])
        acc += number_distribution(sec, Psi).sum(1)
    acc /= len(ts)
    print(f"D draw {trial}: time-averaged p (4000 random t in [0,1e7]) =", np.round(acc, 4), " exact p_bar =", np.round(pb, 4))
