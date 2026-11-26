# Exact annealed purity of a Haar-random process tensor (collision model with memory), with depolarizing noise on E.
# Steps: fresh S input (Choi), U_j Haar on S⊗E, then depolarizing(p) on E. rho_E(0) = I/D.
# Claim: x_{j+1} = a x_j + b y_j,  y_{j+1} = u y_j' ... (see derivation); noiseless: P_n = a^n + (P_E/D)(1-a^n)
import numpy as np
rng = np.random.default_rng(3)
def haar(N):
    Z = (rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))/np.sqrt(2)
    Q, R = np.linalg.qr(Z); return Q * (np.diag(R)/np.abs(np.diag(R)))
def process_tensor(d, D, n, p):
    rho = np.eye(D, dtype=complex)/D                        # registers: (X, E), X empty
    phi = np.zeros(d*d, complex); phi[[i*d+i for i in range(d)]] = 1/np.sqrt(d)
    Phi = np.outer(phi, phi.conj())                         # on (R_j, S)
    paulis_like = [np.eye(D)]                               # depolarizing via twirl over a unitary 1-design (generalized Paulis)
    w = np.exp(2j*np.pi/D); X = np.roll(np.eye(D), 1, axis=0); Zm = np.diag(w**np.arange(D))
    weyl = [np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Zm, b) for a in range(D) for b in range(D)]
    for j in range(n):
        dx = rho.shape[0]//D
        r = rho.reshape(dx, D, dx, D)
        full = np.einsum('aebf,rs->arebsf', r, Phi).reshape(dx, d*d, D, dx, d*d, D)
        # split (R_j,S): index rs -> (R,S)
        full = full.reshape(dx, d, d, D, dx, d, d, D)
        U = haar(d*D).reshape(d, D, d, D)                     # U[s', e', s, e]
        full = np.einsum('tfse,xrseyRSE,TFSE->xrtfyRTF', U, full, U.conj())
        # depolarizing on E: (1-p) rho + p Tr_E(rho) ⊗ I/D  ==  (1-p) rho + p * avg over Weyl
        dep = np.einsum('xrtfyRTf->xrtyRT', full)
        full = (1-p)*full + p*np.einsum('xrtyRT,fF->xrtfyRTF', dep, np.eye(D)/D)
        rho = full.reshape(dx*d*d*D, dx*d*d*D)
    Ups = np.einsum('aebe->ab', rho.reshape(rho.shape[0]//D, D, rho.shape[0]//D, D))
    return np.real(np.trace(Ups @ Ups))                      # normalized purity (Tr Ups = 1 here)
def theory(d, D, n, p, PE):
    q = d*D; a = (D*D-1)/(q*q-1); b = D*(d*d-1)/(q*q-1); u = (1-p)**2
    x, y = 1.0, PE
    for _ in range(n):
        xh = a*x + b*y; yh = y
        x, y = xh, u*yh + (1-u)*xh/D
    return x
for (d, D, p) in [(2, 2, 0.0), (2, 3, 0.0), (2, 2, 0.3)]:
    for n in [1, 2, 3]:
        S = 2000 if n < 3 else 800
        vals = [process_tensor(d, D, n, p) for _ in range(S)]
        print("d=%d D=%d p=%.1f n=%d  MC = %.4f +- %.4f   exact = %.4f" % (d, D, p, n, np.mean(vals), np.std(vals)/np.sqrt(S), theory(d, D, n, p, 1/D)))
