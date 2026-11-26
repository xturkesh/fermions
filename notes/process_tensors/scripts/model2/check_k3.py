import numpy as np
from bottleneck import gate, haar, X, Z
from kreplica import otoc_k
U = gate(0.5, 0.3); I = np.eye(2)
def brute_k(U, A, B, t, D, S, k, p=0.0):
    d = A.shape[0]; N = d*d*D; UU = np.kron(U, np.eye(D)); vals = []
    for _ in range(S):
        Vs = [np.kron(np.eye(d), haar(d*D)) for _ in range(t)]
        O = np.kron(A, np.eye(d*D))
        for s in reversed(range(t)):
            if p > 0:
                T = O.reshape(d, d, D, d, d, D); tr = np.trace(T, axis1=2, axis2=5)
                O = ((1-p)*T + p*np.einsum('abcd,ef->abecdf', tr, np.eye(D))/D).reshape(N, N)
            O = Vs[s].conj().T @ O @ Vs[s]; O = UU.conj().T @ O @ UU
        M = O @ np.kron(B, np.eye(d*D))
        vals.append(np.real(np.trace(np.linalg.matrix_power(M, k)))/N)
    return np.mean(vals), np.std(vals)/np.sqrt(S)
A, B = Z + 0.5*I, X + 0.5*I
for (t, D, p) in [(2, 2, 0.0), (3, 2, 0.0), (3, 2, 0.4), (3, 3, 0.4)]:
    m, e = brute_k(U, A, B, t, D, 2500, 3, p)
    print("k=3 t=%d D=%d p=%.1f  MC=%.4f+-%.4f  exact=%.4f" % (t, D, p, m, e, otoc_k(U, A, B, t, D, 3, p)))
