import numpy as np
from bottleneck import gate, X, Y, Z
from limit import otoc_inf
U = gate(0.5, 0.3); I = np.eye(2)
print("late-time C2 = phi(A B A B), A = Z + a I, B = X + b I  (phi(A)=a, k2(A)=1, phi(B)=b, k2(B)=1)")
print("free: a^2+b^2+a^2 b^2 ; boolean: a^2 b^2 ; 'A-monotone' a^2(1+b^2) ; 'B-monotone' b^2(1+a^2)")
for (a, b) in [(0.5, 0.0), (0.0, 0.5), (0.5, 0.7), (1.0, 1.0)]:
    A = Z + a*I; B = X + b*I
    row = []
    for p in [0.0, 0.02, 0.1, 0.3, 0.6, 1.0]:
        u = (1-p)**2; row.append(otoc_inf(U, A, B, 60, u))
    print("a=%.1f b=%.1f  C2(t=60) for p=0,0.02,0.1,0.3,0.6,1: %s   free=%.4f bool=%.4f" % (
        a, b, np.round(row, 4), a*a+b*b+a*a*b*b, a*a*b*b))
