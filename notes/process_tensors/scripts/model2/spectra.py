import numpy as np
from bottleneck import gate, X, Y, Z
from limit import maps, otoc_inf, single
U = gate(0.5, 0.3); d = 2; Ulayer, R = maps(U, d)
def as_matrix(f):
    M = np.zeros((d**4, d**4), complex)
    for k in range(d**4):
        E = np.zeros(d**4, complex); E[k] = 1
        M[:, k] = f(E.reshape(d, d, d, d)).reshape(-1)
    return M
Me = as_matrix(lambda Xt: R(Ulayer(Xt, 0), 0)/d**2)
MF = as_matrix(lambda Xt: R(Ulayer(Xt, 1), 1)/d**2)
M1 = np.zeros((4, 4), complex)
for k in range(4):
    E = np.zeros(4, complex); E[k] = 1; Xa = E.reshape(2, 2)
    M1[:, k] = (np.einsum('aibi->ab', (U.conj().T @ np.kron(Xa, np.eye(2)) @ U).reshape(2, 2, 2, 2))/2).reshape(-1)
ev = lambda M: np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
print("single-replica channel |eig|:", np.round(ev(M1), 5))
print("M_e (box) |eig| top 6:", np.round(ev(Me)[:6], 5))
print("M_F (echo) |eig| top 6:", np.round(ev(MF)[:6], 5))
lam = ev(M1)[1]; print("lambda = %.5f, lambda^2 = %.5f" % (lam, lam**2))
print("\nC2(t) for A=B=Z and u=(1-p)^2 (noise on E):")
for p in [0.0, 0.05, 0.2, 0.5]:
    u = (1-p)**2; vals = [otoc_inf(U, Z, Z, t, u) for t in [10, 20, 30, 40]]
    rates = [abs(vals[i+1]/vals[i])**(1/10) for i in range(3)]
    print("p=%.2f  C2(10,20,30,40) = %s   rate est = %s" % (p, np.array2string(np.array(vals), precision=3), np.round(rates, 4)))
print("C1(t) rate:", np.round([abs(single(U, Z, Z, t+10)/single(U, Z, Z, t))**0.1 for t in [10, 20]], 4))
