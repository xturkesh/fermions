# Minimal self-contained copy of the model-1 k=2 transfer (from pt_mi.py / pt_local.py)
import numpy as np
def step_matrix(d, D, pin, pout):
    D = float(D); q = d*D; Wg = {0: 1/(q*q-1), 1: -1/(q*(q*q-1))}
    trS = lambda a, b: d*d if a == b else d
    trE = {(0, 0): D*D, (0, 1): D, (1, 0): D, (1, 1): D*D}
    M = np.zeros((2, 2))
    for j in range(2):
        for sig in range(2):
            for tau in range(2):
                w = Wg[sig ^ tau]*trS(tau, pin)*trS(sig, pout)
                if tau == j:
                    M[0, j] += w*trE[(0, sig)]; M[1, j] += w*trE[(1, sig)]
    return M/d**2
def noise_local(q, L, p):
    D = float(q)**int(L); u = (1-p)**2; G = (u*q*q + 1 - u)**int(L)
    return np.array([[1, 0], [(D*D - G)/(D*(D*D-1)), (G-1)/(D*D-1)]])
def Zl(d, q, L, p, pattern):
    D = float(q)**int(L); v = np.array([1.0, 1/D]); Nz = noise_local(q, L, p)
    Ms = {pp: step_matrix(d, D, *pp) for pp in [(0, 0), (1, 1)]}
    for pp in pattern: v = Nz @ (Ms[pp] @ v)
    return v[0]
def I2l(d, q, L, p, n1, n2):
    F, I = (1, 1), (0, 0)
    return (np.log2(Zl(d, q, L, p, [F]*(n1+n2))) - np.log2(Zl(d, q, L, p, [F]*n1+[I]*n2)) - np.log2(Zl(d, q, L, p, [I]*n1+[F]*n2)))
