import numpy as np
from bottleneck import gate, X, Z
from kreplica import otoc_k
U = gate(0.5, 0.3); I = np.eye(2); Dbig = 10**4; lam = 0.46219
a, b = 0.5, 0.5; A, B = Z + a*I, X + b*I
free = a**3*(3*b + b**3) + 3*a*1*b*(1 + b*b); mono = a**3*(3*b + b**3)
print("k=3, d=2: free=%.6f  monotone(A>B)=%.6f" % (free, mono))
for p in [0.0, 0.1, 0.3]:
    vals = {t: otoc_k(U, A, B, t, Dbig, 3, p) for t in [12, 18, 40]}
    r = abs((vals[18]-vals[40])/(vals[12]-vals[40]))**(1/6)
    print("  p=%.1f  C3(12,18,40) = %.6f %.6f %.6f   approach rate=%.4f  (1-p)^2=%.4f lam=%.4f" % (p, vals[12], vals[18], vals[40], r, (1-p)**2, lam))
