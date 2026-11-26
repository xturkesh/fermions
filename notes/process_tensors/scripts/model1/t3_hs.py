# Task 3: annealed Hilbert-Schmidt distance from the Markovian (product) process via k=2,3,4 transfers,
# plus Monte Carlo trace distance for tiny sizes.
import sys, time, numpy as np; sys.path.insert(0, '.')
from common import *
rng = np.random.default_rng(7)
def hs_terms(d, q, L, n1, n2, kraus=None):
    F = transp(2, 0, 1); I = IDP(2)
    m2 = Model(2, d, q, L, kraus); m3 = Model(3, d, q, L, kraus); m4 = Model(4, d, q, L, kraus)
    t12, t13 = transp(3, 0, 1), transp(3, 0, 2); s12, s34 = transp(4, 0, 1), transp(4, 2, 3)
    A = m2.value([(F, F)]*(n1+n2))
    B = m3.value([(t12, t12)]*n1 + [(t13, t13)]*n2)
    C = m4.value([(s12, s12)]*n1 + [(s34, s34)]*n2)
    return A, B, C
if __name__ == "__main__":
    print("Noiseless: delta_HS = E||U - Up x Uf||_2^2 / E Tr U^2  (annealed ratio), d=2;  compare d^{2n}/(D^2 + d^{2n})")
    print("   D    n1=n2=1             n1=n2=2             n1=n2=3             n1=n2=4")
    for L in [1, 2, 3, 4, 5, 6]:
        D = 2**L; row = []
        for n1 in [1, 2, 3, 4]:
            A, B, C = hs_terms(2, 2, L, n1, n1); dl = (A - 2*B + C)/A; n = 2*n1
            row.append("%.4f (%.4f)" % (dl, 4.0**n/(D*D + 4.0**n)))
        print("   %-4d " % D + "  ".join(row))
    print("\nWith local depolarizing noise (long time n1=n2=6), d=2: delta_HS")
    for p in [0.0, 0.1, 0.3]:
        row = []
        for L in [1, 2, 3, 4, 5]:
            kr = depolarizing(2, p) if p > 0 else None
            A, B, C = hs_terms(2, 2, L, 6, 6, kr); row.append("L=%d: %.4f" % (L, (A - 2*B + C)/A))
        print("   p=%.1f  " % p + "  ".join(row))
    print("\nMonte Carlo, tiny sizes: trace distance T = 0.5||U - Up x Uf||_1 and per-sample HS ratio, noiseless, d=2")
    def trnorm(X): return np.sum(np.abs(np.linalg.eigvalsh((X + X.conj().T)/2)))
    t0 = time.time()
    for (n1, L, S) in [(1, 1, 400), (1, 2, 400), (1, 3, 300), (1, 4, 150), (2, 1, 150), (2, 2, 100), (2, 3, 40)]:
        T, R = [], []
        for _ in range(S):
            U = sample_upsilon(2, 2, L, 2*n1, None, rng); Up, Uf = marginals(U, 2, n1, n1); X = U - np.kron(Up, Uf)
            T.append(0.5*trnorm(X)); R.append(np.real(np.trace(X@X))/np.real(np.trace(U@U)))
        A, B, C = hs_terms(2, 2, L, n1, n1)
        print("   n1=n2=%d D=%-3d  <T> = %.4f +- %.4f   <HS ratio> = %.4f   annealed delta_HS = %.4f" % (n1, 2**L, np.mean(T), np.std(T)/np.sqrt(S), np.mean(R), (A-2*B+C)/A))
    print("   (time %.0fs)" % (time.time()-t0))
