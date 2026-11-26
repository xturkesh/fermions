import numpy as np, time
from bottleneck import gate, X, Z
from kreplica import otoc_k
U = gate(0.5, 0.3); I = np.eye(2); b = 0.5; A, B = Z, X + b*I
print("k=4, A=Z (kappa1=0, kappa2=1, kappa4=-1), B=X+%.1f : free late value 2b^2+b^4 = %.6f" % (b, 2*b*b + b**4))
t0 = time.time()
print("  p=0  : C4(14) = %.6f   C4(20) = %.6f" % (otoc_k(U, A, B, 14, 10**4, 4, 0.0), otoc_k(U, A, B, 20, 10**4, 4, 0.0)), "[%.0fs]" % (time.time()-t0))
v1, v2 = otoc_k(U, A, B, 14, 10**4, 4, 0.1), otoc_k(U, A, B, 20, 10**4, 4, 0.1)
print("  p=0.1: C4(14) = %.6f   C4(20) = %.6f   rate = %.4f   (1-p)^4 = %.4f  (1-p)^2 = %.4f" % (v1, v2, abs(v2/v1)**(1/6), 0.9**4, 0.81))
