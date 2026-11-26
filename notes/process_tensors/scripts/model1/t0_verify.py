# Monte Carlo verification of the S_k transfer (k=2,3,4), unital and non-unital noise.
import sys, numpy as np; sys.path.insert(0, '.')
from common import *
rng = np.random.default_rng(11)
def mc(d, q, L, n1, n2, kraus, S, rho_q=None):
    acc = []
    for _ in range(S):
        U = sample_upsilon(d, q, L, n1+n2, kraus, rng, rho_q)
        Up, Uf = marginals(U, d, n1, n2)
        tr = lambda A, k: np.real(np.trace(np.linalg.matrix_power(A, k)))
        cross = np.real(np.trace(U @ np.kron(Up, Uf)))
        acc.append([tr(U, 2), tr(Up, 2), tr(Uf, 2), tr(U, 3), cross, tr(Up, 2)*tr(Uf, 2), tr(U, 4)])
    a = np.array(acc); return a.mean(0), a.std(0)/np.sqrt(S)
def exact(d, q, L, n1, n2, kraus, rho_q=None):
    n = n1 + n2
    m2 = Model(2, d, q, L, kraus, rho_q); F = transp(2, 0, 1); I2 = IDP(2)
    m3 = Model(3, d, q, L, kraus, rho_q); c3 = cyc(3); t12 = transp(3, 0, 1); t13 = transp(3, 0, 2)
    m4 = Model(4, d, q, L, kraus, rho_q); c4 = cyc(4); s12 = transp(4, 0, 1); s34 = transp(4, 2, 3)
    return np.array([m2.value([(F, F)]*n), m2.value([(F, F)]*n1 + [(I2, I2)]*n2), m2.value([(I2, I2)]*n1 + [(F, F)]*n2),
                     m3.value([(c3, c3)]*n), m3.value([(t12, t12)]*n1 + [(t13, t13)]*n2),
                     m4.value([(s12, s12)]*n1 + [(s34, s34)]*n2), m4.value([(c4, c4)]*n)])
labels = ['TrU^2', 'TrUp^2', 'TrUf^2', 'TrU^3', 'Tr[U(Up x Uf)]', 'TrUp^2*TrUf^2', 'TrU^4']
cases = [("noiseless d=2,D=2, n1=n2=1", 2, 2, 1, 1, 1, None, None),
         ("depolarizing p=0.3, D=2, n1=2,n2=1", 2, 2, 1, 2, 1, depolarizing(2, 0.3), None),
         ("dephasing p=0.25, D=4 (L=2), n1=n2=1", 2, 2, 2, 1, 1, dephasing(0.25), None),
         ("amp.damp. g=0.5 (non-unital), D=2, n1=2,n2=1", 2, 2, 1, 2, 1, amp_damp(0.5), None),
         ("amp.damp. g=0.7, D=4, n1=1,n2=2, rho_E(0)=|0><0|", 2, 2, 2, 1, 2, amp_damp(0.7), np.diag([1., 0.]))]
if __name__ == "__main__":
  for (name, d, q, L, n1, n2, kr, rq) in cases:
      S = 1500 if q**L*d**(2*(n1+n2)) <= 256 else 600
      m, e = mc(d, q, L, n1, n2, kr, S, rq); ex = exact(d, q, L, n1, n2, kr, rq)
      print("==", name, " (samples=%d)" % S)
      for lab, a, b, c in zip(labels, m, e, ex):
          print("   %-16s MC %.5f +- %.5f   exact %.5f   (%.1f sigma)" % (lab, a, b, c, abs(a-c)/max(b, 1e-12)))
