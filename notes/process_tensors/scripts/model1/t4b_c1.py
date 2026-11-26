# Free-Poisson test at the 'Page point' N = d^{2n} = M = D^2 (c = 1), exact arithmetic; and delta_HS = 1 - 2^{-I2} identity (unital).
import sys, numpy as np; sys.path.insert(0, '.')
from fractions import Fraction as Fr
from exact import ExactModel
from t4_spectrum_exact import induced, moments, hs
from common import transp, IDP
d = 2
print("c = N/M = 1 (N=d^{2n}=D^2): exact/induced for k=2,3,4")
for (L, n) in [(2, 2), (4, 4), (6, 6), (8, 8)]:
    D = 2**L; mo = moments(d, 2, L, n); N = d**(2*n); M = D*D
    print("   D=%-4d n=%d  " % (D, n) + "  ".join("k=%d: %.5f" % (k, float(mo[i]/induced(k, N, M))) for i, k in enumerate([2, 3, 4])))
print("identity check (unital): delta_HS == 1 - 2^{-I2_annealed}")
F = transp(2, 0, 1); I = IDP(2)
for (L, n1, n2, p) in [(2, 2, 3, None), (3, 3, 2, Fr(1, 10)), (2, 4, 4, Fr(3, 10))]:
    A, B, C = hs(L, n1, n2, p); m = ExactModel(2, d, 2, L, p)
    Pp = m.value([(F, F)]*n1 + [(I, I)]*n2); Pf = m.value([(I, I)]*n1 + [(F, F)]*n2)
    print("   L=%d n1=%d n2=%d p=%s : B==C %s ; delta=%.6f ; 1-Pp*Pf/Pall=%.6f" % (L, n1, n2, p, B == C, float((A-2*B+C)/A), float(1 - Pp*Pf/A)))
