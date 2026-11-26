# Annealed Renyi-2 temporal mutual information I2(past:future) of a Haar-random process tensor with noise on E.
import numpy as np
exec(open('pt_purity.py').read().split("for (d, D, p) in")[0])      # reuse haar(), process_tensor-like code
def step_matrix(d, D, pin, pout):
    q = d*D; Wg = {0: 1/(q*q-1), 1: -1/(q*(q*q-1))}
    trS = lambda a, b: d*d if a == b else d             # Tr(a_S b_S), a,b in {0=I,1=F}
    trE = {(0, 0): D*D, (0, 1): D, (1, 0): D, (1, 1): D*D}
    M = np.zeros((2, 2))                                 # (m_I, m_F) -> (m_I', m_F'), normalized by d^2
    for col, Yb in enumerate([None, None]):
        pass
    for j in range(2):                                   # input moment j: m_j = Tr(tau_E^j ... ) basis
        for sig in range(2):
            for tau in range(2):
                w = Wg[sig ^ tau] * trS(tau, pin) * trS(sig, pout)
                # Tr(tau_E Y) = m_tau ; output sigma_E contributes Tr(I sigma)=trE[0,sig], Tr(F sigma)=trE[1,sig]
                if tau == j:
                    M[0, j] += w * trE[(0, sig)]
                    M[1, j] += w * trE[(1, sig)]
    return M / d**2
def noise(D, p):
    u = (1-p)**2; return np.array([[1, 0], [(1-u)/D, u]])
def Z(d, D, p, pattern, PE):
    v = np.array([1.0, PE])
    for (pin, pout) in pattern:
        v = noise(D, p) @ (step_matrix(d, D, pin, pout) @ v)
    return v[0]
def I2(d, D, p, n1, n2, PE=None):
    PE = 1/D if PE is None else PE
    F, I = (1, 1), (0, 0)
    za = Z(d, D, p, [F]*(n1+n2), PE); zp = Z(d, D, p, [F]*n1 + [I]*n2, PE); zf = Z(d, D, p, [I]*n1 + [F]*n2, PE)
    return np.log2(za) - np.log2(zp) - np.log2(zf), za, zp, zf
# --- Monte Carlo check of the marginal purities (d=2, D=2, n1=1, n2=1 and n1=2,n2=1)
def mc(d, D, n1, n2, p, S):
    out = []
    for _ in range(S):
        n = n1 + n2
        rho = np.eye(D, dtype=complex)/D
        phi = np.zeros(d*d, complex); phi[[i*d+i for i in range(d)]] = 1/np.sqrt(d); Phi = np.outer(phi, phi.conj())
        for j in range(n):
            dx = rho.shape[0]//D; r = rho.reshape(dx, D, dx, D)
            full = np.einsum('aebf,rs->arebsf', r, Phi).reshape(dx, d, d, D, dx, d, d, D)
            U = haar(d*D).reshape(d, D, d, D)
            full = np.einsum('tfse,xrseyRSE,TFSE->xrtfyRTF', U, full, U.conj())
            dep = np.einsum('xrtfyRTf->xrtyRT', full)
            full = (1-p)*full + p*np.einsum('xrtyRT,fF->xrtfyRTF', dep, np.eye(D)/D)
            rho = full.reshape(dx*d*d*D, dx*d*d*D)
        N = rho.shape[0]//D; Ups = np.einsum('aebe->ab', rho.reshape(N, D, N, D))
        P, Fd = d**(2*n1), d**(2*n2)
        U4 = Ups.reshape(P, Fd, P, Fd)
        past = np.einsum('afbf->ab', U4); fut = np.einsum('afag->fg', U4)
        out.append([np.real(np.trace(Ups@Ups)), np.real(np.trace(past@past)), np.real(np.trace(fut@fut))])
    return np.mean(out, axis=0), np.std(out, axis=0)/np.sqrt(S)
for (n1, n2, p) in [(1, 1, 0.0), (2, 1, 0.0), (1, 1, 0.3)]:
    m, e = mc(2, 2, n1, n2, p, 2000)
    _, za, zp, zf = I2(2, 2, p, n1, n2)
    print("n1=%d n2=%d p=%.1f  MC (all,past,fut) = %s +- %s   exact = [%.4f %.4f %.4f]" % (n1, n2, p, np.round(m, 4), np.round(e, 4), za, zp, zf))
print("\nI2(past:future) in bits, d=2, n1=n2=20 (saturated):")
for D in [2, 4, 8, 16, 64]:
    print("D=%3d  " % D + "  ".join("p=%.2f: %.4f" % (p, I2(2, D, p, 20, 20)[0]) for p in [0.0, 0.01, 0.05, 0.2]))
