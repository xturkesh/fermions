# Task 1: general noise, closed 2x2 forms, long-time I2, crossing and avoided-crossing width.
import sys, numpy as np; sys.path.insert(0, '.')
from common import *
def gIF(kraus, q):
    NI = sum(K @ K.conj().T for K in kraus)
    gI = np.real(np.trace(NI @ NI))
    gF = sum(abs(np.trace(Ki.conj().T @ Kj))**2 for Ki in kraus for Kj in kraus)
    return gI, gF
def closed(d, D, GI, GF):
    q2 = (d*D)**2; a = (D*D-1)/(q2-1); b = D*(d*d-1)/(q2-1)
    cI = (D*GI - GF)/(D*(D*D-1)); cF = (D*GF - GI)/(D*(D*D-1))
    TFF = np.array([[a, b], [cI*a, cI*b + cF]]); TII = np.array([[1, 0], [cI + cF*b, cF*a]])
    return a, b, cI, cF, TFF, TII
def I2_finite(d, D, GI, GF, PE, n1, n2):
    a, b, cI, cF, TFF, TII = closed(d, D, GI, GF); v0 = np.array([1.0, PE]); e = np.array([1.0, 0])
    P = lambda M, v: e @ M @ v
    Pall = e @ np.linalg.matrix_power(TFF, n1+n2) @ v0; Pp = e @ np.linalg.matrix_power(TFF, n1) @ v0
    Pf = e @ np.linalg.matrix_power(TFF, n2) @ np.linalg.matrix_power(TII, n1) @ v0
    return np.log2(Pall/(Pp*Pf))
def I2_inf(d, D, GI, GF, PE=None):
    a, b, cI, cF, TFF, TII = closed(d, D, GI, GF)
    ev, R = np.linalg.eig(TFF); i = np.argmax(ev.real); lam = ev[i].real
    evl, Lm = np.linalg.eig(TFF.T); j = np.argmax(evl.real)
    r = R[:, i].real; l = Lm[:, j].real; l = l/(l @ r)
    # w = fixed point of TII acting on (1, y): y* = (cI + cF b)/(1 - cF a)
    ystar = (cI + cF*b)/(1 - cF*a); w = np.array([1.0, ystar])
    return -np.log2((r[0])*(l @ w)), lam, ev
# --- (i) g_I, g_F from the formula vs from the k=2 permutation computation
print("(i) per-qudit g_I = Tr N(I)^2, g_F = sum_ij |Tr K_i^dag K_j|^2  vs  k=2 noise matrix:")
for name, kr, q in [("depol p=0.2 (q=2)", depolarizing(2, 0.2), 2), ("depol p=0.2 (q=3)", depolarizing(3, 0.2), 3),
                    ("dephasing p=0.3", dephasing(0.3), 2), ("amp.damp g=0.4", amp_damp(0.4), 2)]:
    gI, gF = gIF(kr, q); g = qudit_noise_matrix(2, q, kr)
    print("   %-20s formula (gI,gF)=(%.6f, %.6f)   perm-matrix (g[F,e], g[F,F])=(%.6f, %.6f)" % (name, gI, gF, g[1, 0], g[1, 1]))
print("   analytic: depol gF = u q^2 + 1 - u (u=(1-p)^2), gI = q ; dephasing gF = 4[(1-p)^2+p^2], gI = 2 ;"
      " amp.damp gF = (2-g)^2 + g^2, gI = 2(1+g^2)")
# --- (ii) closed 2x2 forms vs general S_2 transfer
print("\n(ii) closed 2x2 transfer vs general Model(k=2): P_all, P_past, P_fut")
F = transp(2, 0, 1); I = IDP(2)
for name, kr, L, n1, n2 in [("depol 0.1", depolarizing(2, 0.1), 3, 4, 3), ("amp.damp 0.3", amp_damp(0.3), 2, 3, 5), ("dephasing 0.2", dephasing(0.2), 3, 2, 2)]:
    m = Model(2, 2, 2, L, kr); D = 2**L; gI, gF = gIF(kr, 2)
    a, b, cI, cF, TFF, TII = closed(2, D, gI**L, gF**L); v0 = np.array([1.0, 1.0/D]); e = np.array([1.0, 0])
    c = [e @ np.linalg.matrix_power(TFF, n1+n2) @ v0, e @ np.linalg.matrix_power(TFF, n1) @ v0,
         e @ np.linalg.matrix_power(TFF, n2) @ np.linalg.matrix_power(TII, n1) @ v0]
    g = [m.value([(F, F)]*(n1+n2)), m.value([(F, F)]*n1 + [(I, I)]*n2), m.value([(I, I)]*n1 + [(F, F)]*n2)]
    print("   %-13s L=%d closed %s  general %s" % (name, L, np.round(c, 8), np.round(g, 8)))
# --- (iii) noiseless temporal Page curve
print("\n(iii) noiseless: I2(n1,n2) = log2[P(n1+n2)/(P(n1)P(n2))], P(n) = a^n + (1-a^n)/D^2 ; examples d=2:")
for D in [4, 16, 64]:
    vals = [I2_finite(2, D, D, D*D, 1/D, n1, n1) for n1 in [1, 2, 3, 4, 6, 10]]
    print("   D=%3d  n1=n2=1,2,3,4,6,10:  " % D + "  ".join("%.3f" % v for v in vals) + "   (2log2D=%.0f)" % (2*np.log2(D)))
# --- (iv) long-time I2 vs L for the three channels; bare crossing L* from cF + cI b = a ; width
print("\n(iv) long-time I2_inf [bits] vs L (q=2, d=2); L* = bare crossing (cF + cI*b = a), Lc = steepest drop, width = L(25%)-L(75%) of I2/(2L)")
def scan(kr, Ls):
    gI, gF = gIF(kr, 2); out = []
    for L in Ls:
        D = 2.0**L; v, lam, ev = I2_inf(2, D, gI**L, gF**L); out.append(v)
    return np.array(out)
for name, mk, ps in [("depolarizing", lambda p: depolarizing(2, p), [0.05, 0.1, 0.2, 0.3]),
                     ("dephasing", dephasing, [0.1, 0.25, 0.5]),
                     ("amp.damping", amp_damp, [0.1, 0.3, 0.6])]:
    for p in ps:
        kr = mk(p); gI, gF = gIF(kr, 2)
        Ls = np.linspace(0.5, 40, 4000); I2 = scan(kr, Ls); frac = I2/(2*Ls)
        # bare crossing
        def bare(L):
            D = 2.0**L; a, b, cI, cF, *_ = closed(2, D, gI**L, gF**L); return cF + cI*b - a
        Lb = [L for L, L2 in zip(Ls[:-1], Ls[1:]) if bare(L) > 0 >= bare(L2)]
        L75 = Ls[np.argmax(frac < 0.75)] if np.any(frac < 0.75) else np.nan
        L25 = Ls[np.argmax(frac < 0.25)] if np.any(frac < 0.25) else np.nan
        Lpred = np.log(4)/np.log(4/gF)
        imax = np.argmax(I2)
        print("   %-12s p=%.2f  gF/q^2=%.4f  L*_bare=%s  L*_pred=ln d^2/ln(q^2/gF)=%.2f  I2max=%.2f bits at L=%.2f  width=%.3f" %
              (name, p, gF/4, ("%.2f" % Lb[0]) if Lb else "none", Lpred, I2[imax], Ls[imax], L25 - L75))
# --- (v) minimum eigenvalue gap at the avoided crossing vs prediction 2 sqrt(a b cI)
print("\n(v) avoided crossing: min |lam+ - lam-| over L vs 2 sqrt(a b c_I) at L*  (depolarizing, q=2)")
for p in [0.3, 0.2, 0.1, 0.05]:
    kr = depolarizing(2, p); gI, gF = gIF(kr, 2); best = (9, None)
    for L in np.linspace(0.5, 40, 20000):
        D = 2.0**L; a, b, cI, cF, TFF, _ = closed(2, D, gI**L, gF**L); ev = np.linalg.eigvals(TFF)
        gap = abs(ev[0] - ev[1])
        if gap < best[0]: best = (gap, L, 2*np.sqrt(a*b*cI), D)
    print("   p=%.2f  L*=%.2f  D=%.1f  min gap=%.3e  2sqrt(a b cI)=%.3e" % (p, best[1], best[3], best[0], best[2]))
