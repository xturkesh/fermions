import sys, numpy as np
sys.path.insert(0, '../'); sys.path.insert(0, '../temporal_ipr')
from m1core import I2l
from records_m1 import exact
d, q = 2, 2
out = {}
# --- P1: long-time I2 vs L for several p (exact finite n1=n2=30)
Ls = np.arange(1, 17); ps = [0.0, 0.05, 0.1, 0.3]
out['P1_L'] = Ls
for p in ps: out['P1_p%.2f' % p] = np.array([I2l(d, q, L, p, 30, 30) for L in Ls])
# --- P2: scaling collapse, closed-form long-time I2 (fork formula)
def I2inf(L, p):
    D = float(q)**L; u = (1-p)**2; GI, GF = D, (u*q*q + 1 - u)**L
    a = (D*D-1)/(d*d*D*D-1); b = D*(d*d-1)/(d*d*D*D-1)
    cI = (D*GI - GF)/(D*(D*D-1)); cF = (D*GF - GI)/(D*(D*D-1))
    T = np.array([[a, b], [cI*a, cI*b + cF]])
    w, R = np.linalg.eig(T); i = np.argmax(w.real); r = R[:, i].real
    wl, Lv = np.linalg.eig(T.T); j = np.argmax(wl.real); l = Lv[:, j].real
    l = l/(l @ r); ys = (cI + cF*b)/(1 - cF*a); wv = np.array([1.0, ys])
    return -np.log2((r[0])*(l @ wv))
# consistency check against finite-n exact at moderate L
for L, p in [(4, 0.1), (8, 0.1), (3, 0.3)]:
    print("check I2inf vs exact(n=60):", L, p, round(I2inf(L, p), 4), round(I2l(d, q, L, p, 60, 60), 4))
for p in [0.1, 0.03, 0.01, 0.005]:
    Lmax = int(1.6/p); Lg = np.arange(1, Lmax+1)
    out['P2_x_p%.3f' % p] = Lg*p
    out['P2_y_p%.3f' % p] = np.array([I2inf(L, p)/(2*L) for L in Lg])
# --- P6: normalized temporal IPR vs n (closed form), pure and mixed bath
ns = np.unique(np.round(np.logspace(0, np.log10(400), 60)).astype(int)); out['P6_n'] = ns
for D in [1, 2, 8, 32, 128]:
    lp = (D+1)/(d*D+1); lm = (D-1)/(d*D-1) if D > 1 else 0.0
    for tag, PE in [('pure', 1.0), ('mixed', 1.0/D)]:
        out['P6_%s_D%d' % (tag, D)] = np.array([float(d)**int(n)*((1+PE)/2*lp**int(n) + (1-PE)/2*lm**int(n)) for n in ns])
# --- P8: record memory witness vs n1 (pure bath)
n1s = np.unique(np.round(np.logspace(0, np.log10(600), 30)).astype(int)); out['P8_n1'] = n1s
for D in [2, 8, 32, 128]:
    vals = []
    for n1 in n1s:
        R = lambda rec: exact(d, D, 0.0, 2*n1, 1.0, rec)
        vals.append(np.log2(R([True]*(2*n1))/(R([True]*n1+[False]*n1)*R([False]*n1+[True]*n1))))
    out['P8_D%d' % D] = np.array(vals)
np.savez('data_m1_ipr.npz', **out); print("saved", len(out), "arrays")
