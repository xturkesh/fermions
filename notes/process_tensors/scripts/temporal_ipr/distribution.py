# Distribution of record probabilities N*p(x) for model 1 with a pure bath (state-vector branching).
import numpy as np
rng = np.random.default_rng(8)
def haar(N):
    Z = (rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))/np.sqrt(2)
    Q, R = np.linalg.qr(Z); return Q*(np.diag(R)/np.abs(np.diag(R)))
def record_probs(d, D, n):
    psi = np.zeros((1, d*D), complex); psi[0, 0] = 1          # S=|0>, E=|0>
    for t in range(n):
        U = haar(d*D); psi = psi @ U.T                         # apply U to every branch
        psi = psi.reshape(-1, d, D)
        new = np.zeros((psi.shape[0]*d, d, D), complex)
        for o in range(d):                                     # project S on |o>, S -> |o>
            new[o::d, o, :] = psi[:, o, :]
        psi = new.reshape(-1, d*D)
    return np.sum(np.abs(psi)**2, axis=1)                      # p(x), x in d^n records
if __name__ == "__main__":
  d = 2
  print("statistics of ln(N p): PT reference mean=-0.577 var=1.645 ; uniform: var=0 ; log-normal: var grows ~ n")
  for D in [1, 4, 16, 64]:
      for n in [4, 8, 12]:
          S = 30 if n == 12 else 60
          lnNp, niprs = [], []
          for _ in range(S):
              p = record_probs(d, D, n); N = d**n
              lnNp.append(np.log(N*p[p > 0])); niprs.append(N*np.sum(p*p))
          lnNp = np.concatenate(lnNp); niprs = np.array(niprs)
          lp = (D+1)/(d*D+1); ann = (d*lp)**n
          print("D=%2d n=%2d  <ln Np>=%6.3f var=%6.3f | N*IPR annealed(exact)=%7.3f  MC mean=%7.3f  quenched exp<ln>=%7.3f" % (
              D, n, lnNp.mean(), lnNp.var(), ann, niprs.mean(), np.exp(np.mean(np.log(niprs)))))
