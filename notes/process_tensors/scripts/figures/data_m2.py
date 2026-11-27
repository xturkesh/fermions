import sys, time, numpy as np
sys.path.insert(0, '../model2')
from bottleneck import gate, X, Z
from limit import otoc_inf
from kreplica import otoc_k
U = gate(0.5, 0.3); I = np.eye(2); lam = 0.46219
out = {'lam': lam}; t0 = time.time()
# (a) |C2(t)| for A=B=Z, several p
ts = np.arange(1, 31); out['A_t'] = ts
for p in [0.0, 0.2, 0.7]:
    out['A_p%.1f' % p] = np.array([otoc_inf(U, Z, Z, int(t), (1-p)**2) for t in ts])
# (b) relaxation rates vs p
def rate(A, B, u, t1=14, t2=24):
    inf = otoc_inf(U, A, B, 200, u)
    g1, g2 = otoc_inf(U, A, B, t1, u) - inf, otoc_inf(U, A, B, t2, u) - inf
    return abs(g2/g1)**(1/(t2-t1))
pg = np.round(np.arange(0.04, 0.96, 0.04), 2); out['B_p'] = pg
out['B_traceless'] = np.array([rate(Z, X + 0.5*I, (1-p)**2) for p in pg])
out['B_generic'] = np.array([rate(Z + 0.5*I, X + 0.5*I, (1-p)**2) for p in pg])
print("rates done %.0fs" % (time.time()-t0))
# (c) C2(t) approach to plateaus, A=Z+0.5, B=X+0.5
tc = np.arange(1, 121); out['C_t'] = tc
for p in [0.0, 0.01, 0.03, 0.1]:
    out['C_p%.2f' % p] = np.array([otoc_inf(U, Z + 0.5*I, X + 0.5*I, int(t), (1-p)**2) for t in tc])
print("plateaus done %.0fs" % (time.time()-t0))
# (d) cumulant spectroscopy at p = 0.1
p = 0.1
t2 = np.arange(2, 41, 2); out['D2_t'] = t2
out['D2'] = np.array([otoc_inf(U, Z, X + 0.5*I, int(t), (1-p)**2) for t in t2])
rng = np.random.default_rng(5)
Zm = (rng.normal(size=(9, 9)) + 1j*rng.normal(size=(9, 9)))/np.sqrt(2); Q, R = np.linalg.qr(Zm); U3 = Q*(np.diag(R)/abs(np.diag(R)))
A3 = np.diag([1., 1., -2.]).astype(complex); B3 = np.diag([1., 0., 0.]).astype(complex)
t3 = np.arange(4, 25, 4); out['D3_t'] = t3
out['D3'] = np.array([otoc_k(U3, A3, B3, int(t), 10**4, 3, p) for t in t3])
print("k=3 done %.0fs" % (time.time()-t0))
t4 = np.arange(4, 21, 4); out['D4_t'] = t4
out['D4'] = np.array([otoc_k(U, Z, X + 0.5*I, int(t), 10**4, 4, p) for t in t4])
print("k=4 done %.0fs" % (time.time()-t0))
np.savez('data_m2.npz', **out)
