import numpy as np
from bottleneck import gate, X, Y, Z
from limit import otoc_inf
U = gate(0.5, 0.3); I = np.eye(2); lam = 0.46219
def rate(A, B, u, t1=14, t2=24):
    inf = otoc_inf(U, A, B, 200, u)
    g1, g2 = otoc_inf(U, A, B, t1, u) - inf, otoc_inf(U, A, B, t2, u) - inf
    return abs(g2/g1)**(1/(t2-t1)), inf
print("lambda=%.4f lambda^2=%.4f ; predicted crossing (1-p)^2 = lambda^2 -> p_c = %.3f" % (lam, lam**2, 1-lam))
for (lab, A, B) in [("A traceless, B=X+0.5", Z, X+0.5*I), ("A=Z+0.5, B=X+0.5", Z+0.5*I, X+0.5*I)]:
    print(lab)
    for p in [0.0, 0.1, 0.3, 0.45, 0.5, 0.55, 0.6, 0.7, 0.9]:
        u = (1-p)**2; r, inf = rate(A, B, u)
        print("   p=%.2f u=%.4f  relaxation rate=%.4f  max(u,lam^2)=%.4f  max(u,lam)=%.4f  plateau=%.4f" % (p, u, r, max(u, lam**2), max(u, lam), inf))
