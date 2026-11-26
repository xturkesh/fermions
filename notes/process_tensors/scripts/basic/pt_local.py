import numpy as np
exec(open('pt_mi.py').read().split("# --- Monte Carlo check")[0])
def noise_local(q, L, p):                 # depolarizing(p) on each of L env qudits, in (m_I, m_F) moment space
    D = q**L; u = (1-p)**2; G = (u*q*q + 1 - u)**L
    return np.array([[1, 0], [(D*D - G)/(D*(D*D-1)), (G-1)/(D*D-1)]])
def Zl(d, q, L, p, pattern):
    D = q**L; v = np.array([1.0, 1/D])
    for (pin, pout) in pattern: v = noise_local(q, L, p) @ (step_matrix(d, D, pin, pout) @ v)
    return v[0]
def I2l(d, q, L, p, n1, n2):
    F, I = (1, 1), (0, 0)
    return (np.log2(Zl(d, q, L, p, [F]*(n1+n2))) - np.log2(Zl(d, q, L, p, [F]*n1+[I]*n2)) - np.log2(Zl(d, q, L, p, [I]*n1+[F]*n2)))
print("I2(past:future) [bits], d=2, q=2, n1=n2=30, local depolarizing p on each of L env qubits")
for p in [0.0, 0.01, 0.05, 0.1, 0.3]:
    print("p=%.2f  " % p + "  ".join("L=%d: %6.3f" % (L, I2l(2, 2, L, p, 30, 30)) for L in [1, 2, 4, 6, 8, 10, 12]))
