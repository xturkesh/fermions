# Task 2: quenched vs annealed Renyi-2 temporal MI (per-sample exact 2-replica propagation), + von Neumann MI for tiny sizes.
import sys, time, numpy as np; sys.path.insert(0, '.')
from common import *
rng = np.random.default_rng(5)
def apply_env_kraus(T, kraus, L, ket_axes, bra_axes):
    for l in range(L):
        out = 0
        for K in kraus:
            t = np.moveaxis(np.tensordot(K, T, axes=(1, ket_axes[l])), 0, ket_axes[l])
            t = np.moveaxis(np.tensordot(t, K.conj(), axes=(bra_axes[l], 1)), -1, bra_axes[l])
            out = out + t
        T = out
    return T
def noise2(Y, kraus, q, L):        # Y[e1,e2,f1,f2] on E⊗E, apply N⊗N
    D = q**L; T = Y.reshape([q]*(4*L))
    T = apply_env_kraus(T, kraus, L, list(range(0, L)), list(range(2*L, 3*L)))          # replica 1
    T = apply_env_kraus(T, kraus, L, list(range(L, 2*L)), list(range(3*L, 4*L)))        # replica 2
    return T.reshape(D, D, D, D)
def noise1(rho, kraus, q, L):
    D = q**L; T = rho.reshape([q]*(2*L))
    return apply_env_kraus(T, kraus, L, list(range(L)), list(range(L, 2*L))).reshape(D, D)
def stepF(Y, U, d):                # U[a,b,s,e] (out S,E ; in S,E); F contraction of S legs in and out
    Uc = U.conj()
    return np.einsum('abse,ABSE,eEfF,AgSf,aGsF->bBgG', U, U, Y, Uc, Uc, optimize=True)/d**2
def stepI1(rho, U, d):             # single replica, maximally mixed S input, S output traced
    return np.einsum('abse,ef,agsf->bg', U, rho, U.conj(), optimize=True)/d
def sample_I2(d, q, L, n1, n2, kraus):
    D = q**L; Us = [haar(d*D, rng).reshape(d, D, d, D) for _ in range(n1+n2)]
    rho0 = np.eye(D, dtype=complex)/D; Y = np.einsum('ef,EF->eEfF', rho0, rho0)
    noisy = kraus is not None and len(kraus) > 1
    for j in range(n1+n2):
        Y = stepF(Y, Us[j], d)
        if noisy: Y = noise2(Y, kraus, q, L)
        if j == n1-1: Pp = np.real(np.einsum('eEeE->', Y))
    Pa = np.real(np.einsum('eEeE->', Y))
    rho = rho0
    for j in range(n1):
        rho = stepI1(rho, Us[j], d)
        if noisy: rho = noise1(rho, kraus, q, L)
    Y = np.einsum('ef,EF->eEfF', rho, rho)
    for j in range(n1, n1+n2):
        Y = stepF(Y, Us[j], d)
        if noisy: Y = noise2(Y, kraus, q, L)
    Pf = np.real(np.einsum('eEeE->', Y))
    return Pa, Pp, Pf
def annealed(d, q, L, n1, n2, kraus):
    m = Model(2, d, q, L, kraus); F = transp(2, 0, 1); I = IDP(2)
    return m.value([(F, F)]*(n1+n2)), m.value([(F, F)]*n1 + [(I, I)]*n2), m.value([(I, I)]*n1 + [(F, F)]*n2)
d, q, n1, n2 = 2, 2, 8, 8
print("Renyi-2 temporal MI, d=2, q=2, n1=n2=%d (long-time), local depolarizing p; S samples per point" % n1)
print("   L  p     annealed I2   quenched <I2>   median I2   std    | check <P_all>/annealed")
t0 = time.time()
for p in [0.0, 0.2, 0.3, 0.42, 0.55]:
    for L in [1, 2, 3]:
        kr = depolarizing(q, p) if p > 0 else None
        S = 300 if L < 3 else 120
        vals = np.array([sample_I2(d, q, L, n1, n2, kr) for _ in range(S)])
        Ia = np.log2(np.mean(vals[:, 0])/(np.mean(vals[:, 1])*np.mean(vals[:, 2])))
        an = annealed(d, q, L, n1, n2, kr); Iex = np.log2(an[0]/(an[1]*an[2]))
        Iq = np.log2(vals[:, 0]/(vals[:, 1]*vals[:, 2]))
        print("   %d  %.2f   %8.3f      %8.3f        %8.3f   %.3f  | %.3f" % (L, p, Iex, Iq.mean(), np.median(Iq), Iq.std(), np.mean(vals[:, 0])/an[0]))
print("   (time %.0fs)" % (time.time()-t0))
# --- von Neumann vs Renyi-2 for tiny sizes, n1=n2=2, full Upsilon
def S_vn(A):
    ev = np.linalg.eigvalsh((A + A.conj().T)/2); ev = ev[ev > 1e-14]; return -np.sum(ev*np.log2(ev))
print("\nTiny sizes, n1=n2=2 (full Upsilon, 4^4=256 dim): quenched von Neumann MI vs Renyi-2 (quenched, annealed)")
print("   L  p     I_vN(quenched)   I2(quenched)   I2(annealed exact)")
for L in [1, 2]:
    for p in [0.0, 0.2, 0.42, 0.7]:
        kr = depolarizing(q, p) if p > 0 else None
        S = 150 if L == 1 else 60; Iv, I2s, P = [], [], []
        for _ in range(S):
            U = sample_upsilon(d, q, L, 4, kr, rng); Up, Uf = marginals(U, d, 2, 2)
            Iv.append(S_vn(Up) + S_vn(Uf) - S_vn(U))
            pa, pp, pf = [np.real(np.trace(X@X)) for X in (U, Up, Uf)]; I2s.append(np.log2(pa/(pp*pf)))
        an = annealed(d, q, L, 2, 2, kr)
        print("   %d  %.2f   %8.3f +- %.3f   %8.3f       %8.3f" % (L, p, np.mean(Iv), np.std(Iv)/np.sqrt(S), np.mean(I2s), np.log2(an[0]/(an[1]*an[2]))))
