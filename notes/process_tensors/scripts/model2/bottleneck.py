# Bottleneck model (Fritzsch-Claeys type): sites A (d), C (d), E (D).
# Step (forward): rho -> V_s (U rho U^dag) V_s^dag, then noise. U fixed on A⊗C, V_s Haar on C⊗E.
# OTOC C2(t) = Tr(A(t) B A(t) B)/(d^2 D), A(t) = Heisenberg (noisy) evolution, noise per copy.
import numpy as np
rng = np.random.default_rng(11)
def haar(N):
    Z = (rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))/np.sqrt(2)
    Q, R = np.linalg.qr(Z); return Q*(np.diag(R)/np.abs(np.diag(R)))
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
def expmh(H):                                               # exp(-i H), H Hermitian
    w, v = np.linalg.eigh(H); return (v * np.exp(-1j*w)) @ v.conj().T
def gate(tau, J, h=(0.3, 0.7)):
    H = np.kron(X, X) + np.kron(Y, Y) + J*np.kron(Z, Z)
    r1 = expmh(h[0]*X + 0.4*Z); r2 = expmh(h[1]*Y - 0.2*Z)
    return np.kron(r1, r2) @ expmh(tau*H) @ np.kron(r2.conj().T, r1)

# ---------- exact 2-replica recursion (any D), labels tau in {0=e,1=F} on C⊗E ----------
def perm_C(d, t):
    P = np.eye(d*d).reshape(d, d, d, d)                    # [c1',c2',c1,c2]
    return P if t == 0 else P.transpose(1, 0, 2, 3)          # swap: <c1' c2'|F|c1 c2> = d(c1',c2)d(c2',c1)
def dep_on_pair(Op, p, which, d):
    """apply N^dag ⊗ N^dag (depolarizing p, unital) on legs 'which' (A or C) of a 2-copy op on (a1 c1 a2 c2)."""
    if p == 0: return Op
    T = Op.reshape([d]*8)                                    # out a1 c1 a2 c2, in a1 c1 a2 c2
    i1, i2 = (0, 2) if which == 'A' else (1, 3)
    def one(T, i):
        tr = np.trace(T, axis1=i, axis2=i+4)                 # Tr over that leg
        tr = np.expand_dims(np.expand_dims(tr, i), i+4) * np.eye(d).reshape([d if k in (i, i+4) else 1 for k in range(8)]) / d
        return (1-p)*T + p*tr
    return one(one(T, i1), i2).reshape(Op.shape)
def recursion(U, A, B, t, D, pE=0.0, pC=0.0, pA=0.0):
    d = A.shape[0]; q = d*D
    Wg = {0: 1/(q*q-1), 1: -1/(q*(q*q-1))}
    uE = (1-pE)**2
    def TE(tau, sig):                                        # Tr(N^dag⊗N^dag(tau_E) sig_E)
        if tau == 0: return D*D if sig == 0 else D
        return uE*D*D + (1-uE) if sig == 1 else D        # F -> uE F + (1-uE) I/D
    Uc = U.conj().T
    # X_tau: operators on (a1,a2) as [a1',a2',a1,a2]
    Xs = {0: np.einsum('ij,kl->ikjl', A, A), 1: np.zeros((d, d, d, d), complex)}
    def U_layer(Xt, tau):                                    # (U^dag⊗U^dag)(X_tau ⊗ tau_C)(U⊗U) on (a1 c1 a2 c2)
        Pc = perm_C(d, tau)
        full = np.einsum('ikjl,mnop->imknjolp', Xt, Pc).reshape(d, d, d, d, d, d, d, d)   # [a1 c1 a2 c2 | a1 c1 a2 c2]
        U4 = U.reshape(d, d, d, d); Ud4 = Uc.reshape(d, d, d, d)
        full = np.einsum('abxy,xyzwpqrs->abzwpqrs', Ud4, full)                   # copy1 left
        full = np.einsum('cdzw,abzwpqrs->abcdpqrs', Ud4, full)                   # copy2 left
        full = np.einsum('abcdxyrs,xypq->abcdpqrs', full, U4)                    # copy1 right
        full = np.einsum('abcdpqxy,xyrs->abcdpqrs', full, U4)                    # copy2 right
        return full.reshape(d**4, d**4)
    for step in range(t):
        Ys = {tau: U_layer(Xs[tau], tau) for tau in (0, 1)}
        if step == t-1: break
        # next (lower) step begins with noise on C/A (Heisenberg order: noise, then V, then U)
        Ys = {tau: dep_on_pair(dep_on_pair(Ys[tau], pC, 'C', d), pA, 'A', d) for tau in (0, 1)}
        new = {0: np.zeros((d, d, d, d), complex), 1: np.zeros((d, d, d, d), complex)}
        for tau in (0, 1):
            T = Ys[tau].reshape([d]*8)
            for sig in (0, 1):
                Pc = perm_C(d, sig)                          # Tr_{C1C2}[Y (I_A ⊗ sig_C)]
                red = np.einsum('aibjckdl,klij->abcd', T, Pc)                     # [a1' a2' a1 a2]
                for taup in (0, 1):
                    new[taup] += Wg[sig ^ taup] * TE(tau, sig) * red
        Xs = new
    # bottom: Tr[Y_tau (B⊗I ⊗ B⊗I)(F_AC)] * Tr(tau_E F_E) / (d^2 D)
    BB = np.kron(np.kron(B, I2 if d == 2 else np.eye(d)), np.kron(B, np.eye(d)))
    Fac = np.eye(d**4).reshape(d*d, d*d, d*d, d*d).transpose(1, 0, 2, 3).reshape(d**4, d**4)
    tot = sum(np.trace(Ys[tau] @ BB @ Fac) * (D if tau == 0 else D*D) for tau in (0, 1))   # Tr(tau_E F_E), no noise after last step
    return np.real(tot)/(d*d*D)

# ---------- brute force (finite D, Monte Carlo over V) ----------
def brute(U, A, B, t, D, S, pE=0.0, pC=0.0, pA=0.0):
    d = A.shape[0]; N = d*d*D; vals = []
    UU = np.kron(U, np.eye(D))
    def Ndag(O, p, site):
        if p == 0: return O
        T = O.reshape(d, d, D, d, d, D); ax = {'A': 0, 'C': 1, 'E': 2}[site]; dim = T.shape[ax]
        tr = np.trace(T, axis1=ax, axis2=ax+3)
        tr = np.expand_dims(np.expand_dims(tr, ax), ax+3) * np.eye(dim).reshape([dim if k in (ax, ax+3) else 1 for k in range(6)]) / dim
        return ((1-p)*T + p*tr).reshape(N, N)
    for _ in range(S):
        Vs = [np.kron(np.eye(d), haar(d*D)) for _ in range(t)]
        O = np.kron(A, np.eye(d*D))
        for s in reversed(range(t)):
            O = Ndag(Ndag(Ndag(O, pE, 'E'), pC, 'C'), pA, 'A')
            O = Vs[s].conj().T @ O @ Vs[s]; O = UU.conj().T @ O @ UU
        BB = np.kron(B, np.eye(d*D))
        vals.append(np.real(np.trace(O @ BB @ O @ BB))/N)
    return np.mean(vals), np.std(vals)/np.sqrt(S)
if __name__ == "__main__":
    U = gate(0.5, 0.3); A, B = Z, Z
    for (D, pE, pC) in [(2, 0, 0), (3, 0, 0), (2, 0.3, 0), (2, 0, 0.3)]:
        for t in [1, 2, 3]:
            m, e = brute(U, A, B, t, D, 3000 if t < 3 else 1500, pE=pE, pC=pC)
            print("D=%d pE=%.1f pC=%.1f t=%d  MC=%.4f+-%.4f  exact=%.4f" % (D, pE, pC, t, m, e, recursion(U, A, B, t, D, pE, pC)))
