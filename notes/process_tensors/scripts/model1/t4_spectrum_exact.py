# Task 4 (exact arithmetic): E Tr U^k vs induced measure (free Poisson) ; HS distance at large D re-checked exactly.
import sys, time, numpy as np; sys.path.insert(0, '.')
from fractions import Fraction as Fr
from exact import ExactModel
from common import perms, compose, inverse, ncyc, cyc, transp, IDP
def induced(k, N, M):
    g = cyc(k); num = sum(Fr(N)**ncyc(s) * Fr(M)**ncyc(compose(inverse(s), g)) for s in perms(k))
    den = Fr(1)
    for j in range(k): den *= (Fr(N)*Fr(M) + j)
    return num/den
def moments(d, q, L, n, p=None):
    out = []
    for k in [2, 3, 4]:
        m = ExactModel(k, d, q, L, p); c = cyc(k); out.append(m.value([(c, c)]*n))
    return out
d = 2; t0 = time.time()
print("Noiseless (exact): ratio E Tr U^k / induced(N=d^{2n}, M=D^2), k=2,3,4")
for L in [1, 2, 3, 4, 6]:
    D = 2**L; row = []
    for n in [1, 2, 4, 8]:
        mo = moments(d, 2, L, n); N = d**(2*n); M = D*D
        row.append(" ".join("%.3f" % float(mo[i]/induced(k, N, M)) for i, k in enumerate([2, 3, 4])))
    print("   D=%-3d  " % D + "   ".join("n=%d: %s" % (n, r) for n, r in zip([1, 2, 4, 8], row)))
print("\nNoisy (exact, local depolarizing p): M_eff from purity, then k=3,4 exact/predicted (free-Poisson test), d=2, n=8")
def Meff(P, N): return (Fr(N) - P)/(P*Fr(N) - 1)
for p in [Fr(1, 20), Fr(1, 5), Fr(2, 5)]:
    for L in [1, 2, 3]:
        mo = moments(d, 2, L, 8, p); N = d**16; Me = Meff(mo[0], N)
        print("   p=%.2f L=%d  M_eff=%10.2f (D^2=%3d)  k=3: %.4f   k=4: %.4f" % (float(p), L, float(Me), 4**L, float(mo[1]/induced(3, N, Me)), float(mo[2]/induced(4, N, Me))))
print("\nExact delta_HS = (A - 2B + C)/A  (A=E TrU^2, B=E Tr[U(Up x Uf)], C=E[TrUp^2 TrUf^2]), noiseless, d=2")
F = transp(2, 0, 1); t12, t13 = transp(3, 0, 1), transp(3, 0, 2); s12, s34 = transp(4, 0, 1), transp(4, 2, 3)
def hs(L, n1, n2, p=None):
    A = ExactModel(2, d, 2, L, p).value([(F, F)]*(n1+n2))
    B = ExactModel(3, d, 2, L, p).value([(t12, t12)]*n1 + [(t13, t13)]*n2)
    C = ExactModel(4, d, 2, L, p).value([(s12, s12)]*n1 + [(s34, s34)]*n2)
    return A, B, C
for (n1, L) in [(1, 6), (1, 8), (2, 6), (2, 8), (3, 6), (3, 9)]:
    A, B, C = hs(L, n1, n1); D = 2**L; dl = (A - 2*B + C)/A
    print("   n1=n2=%d D=%-4d delta=%.6f  delta*D^2/d^{2n}=%.4f  (1-d^-2n1)^2=%.4f   B==C? %s" % (n1, D, float(dl), float(dl)*D*D/4.0**(2*n1), (1-4.0**-n1)**2, B == C))
print("\nExact delta_HS with local depolarizing noise, n1=n2=6")
for p in [Fr(1, 10), Fr(3, 10)]:
    print("   p=%.1f  " % float(p) + "  ".join("L=%d: %.4f" % (L, float((lambda A, B, C: (A-2*B+C)/A)(*hs(L, 6, 6, p)))) for L in [1, 2, 3, 4, 5]))
print("(time %.0fs)" % (time.time()-t0))
