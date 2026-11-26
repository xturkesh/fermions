import numpy as np
from bottleneck import gate, recursion, X, Z
U = gate(0.5, 0.3); I = np.eye(2); a, b = 0.5, 0.5; A, B = Z + a*I, X + b*I
free, mono = a*a + b*b + a*a*b*b, a*a*(1 + b*b)
print("late-time C2 (t=60, D=1e4): free=%.4f monotone=%.4f" % (free, mono))
for (lab, kw) in [("noise on E p=0.3", dict(pE=0.3)), ("noise on C p=0.3", dict(pC=0.3)), ("noise on C p=0.8", dict(pC=0.8)), ("noise on A p=0.3", dict(pA=0.3))]:
    print("  %-18s C2(30)=%.4f  C2(60)=%.4f" % (lab, recursion(U, A, B, 30, 10**4, **kw), recursion(U, A, B, 60, 10**4, **kw)))
