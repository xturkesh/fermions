# Scaling limit p->0, L->inf at fixed x = pL: I2_inf/(2 log2 D) becomes a step at x* = ln(d^2)/kappa.
import numpy as np
def I2_inf_scaled(d, q, L, rF, rI_u):
    """rF = g_F/q^2, rI_u = g_I/q (=1 for unital). Uses normalized combinations to avoid overflow."""
    D = float(q)**L; iD2 = D**-2
    a = (1 - iD2)/(d*d - iD2); b = (d*d - 1)/(D*(d*d - iD2))
    GF_D2 = rF**L; GI_D = rI_u**L                     # G_F/D^2 and G_I/D
    cF = (GF_D2 - GI_D*iD2)/(1 - iD2)                 # (D G_F - G_I)/(D(D^2-1))
    cI = (GI_D - GF_D2/D)/(D*(1 - iD2))               # (D G_I - G_F)/(D(D^2-1))
    T = np.array([[a, b], [cI*a, cI*b + cF]])
    ev, R = np.linalg.eig(T); i = np.argmax(ev.real); r = R[:, i].real
    evl, Lm = np.linalg.eig(T.T); j = np.argmax(evl.real); l = Lm[:, j].real; l = l/(l @ r)
    ystar = (cI + cF*b)/(1 - cF*a); w = np.array([1.0, ystar])
    return -np.log2(r[0]*(l @ w))/(2*np.log2(D))
d = q = 2
print("depolarizing: x* = ln(d^2)/(2(1-q^-2)) = %.4f ;  I2_inf/(2 log2 D) at x = pL:" % (np.log(d*d)/(2*(1-q**-2))))
xs = [0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.3]
print("   p      " + "  ".join("x=%.2f" % x for x in xs))
for p in [0.1, 0.03, 0.01, 0.005]:
    rF = (1-p)**2 + (1-(1-p)**2)/q**2
    row = [I2_inf_scaled(d, q, x/p, rF, 1.0) for x in xs]
    print("   %.3f  " % p + "  ".join("%6.3f" % v for v in row))
print("dephasing: rF=(1-p)^2+p^2, x*=ln(d^2)/2 = %.4f" % (np.log(4)/2))
for p in [0.03, 0.01, 0.005]:
    rF = (1-p)**2 + p**2; print("   %.3f  " % p + "  ".join("%6.3f" % I2_inf_scaled(d, q, x/p, rF, 1.0) for x in [0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 1.0]) + "   (x=0.5,0.6,0.65,0.7,0.75,0.8,1.0)")
print("amplitude damping (non-unital): rF=1-g+g^2/2, g_I/q = 1+g^2, x*=ln(d^2)/1 = %.4f" % np.log(4))
for g in [0.03, 0.01, 0.005]:
    rF = 1 - g + g*g/2; print("   %.3f  " % g + "  ".join("%6.3f" % I2_inf_scaled(d, q, x/g, rF, 1 + g*g) for x in [1.0, 1.2, 1.3, 1.35, 1.4, 1.45, 1.6]) + "   (x=1.0,1.2,1.3,1.35,1.4,1.45,1.6)")
