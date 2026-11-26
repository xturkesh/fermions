# per-sample check: Tr[U (Up x Uf)] vs TrUp^2 TrUf^2 ; and re-run of the non-unital case with a new seed
import sys, numpy as np; sys.path.insert(0, '.')
from common import *
rng = np.random.default_rng(2024)
for (name, q, L, n1, n2, kr) in [("noiseless D=2 n1=1 n2=2", 2, 1, 1, 2, None), ("noiseless D=4 n1=2 n2=2", 2, 2, 2, 2, None),
                                  ("depol 0.2 D=2 n1=2 n2=2", 2, 1, 2, 2, depolarizing(2, 0.2)), ("ampdamp 0.5 D=2 n1=1 n2=2", 2, 1, 1, 2, amp_damp(0.5))]:
    diffs = []
    for _ in range(20):
        U = sample_upsilon(2, q, L, n1+n2, kr, rng); Up, Uf = marginals(U, 2, n1, n2)
        c = np.real(np.trace(U @ np.kron(Up, Uf))); p = np.real(np.trace(Up@Up)*np.trace(Uf@Uf))
        diffs.append((c - p)/p)
    print("%-28s per-sample (cross-prod)/prod: max |.| = %.2e" % (name, np.max(np.abs(diffs))))
# re-run amp damp g=0.5 with more samples, new seed
from t0_verify import mc, exact, labels
import t0_verify; t0_verify.rng = np.random.default_rng(99)
m, e = mc(2, 2, 1, 2, 1, amp_damp(0.5), 3000); ex = exact(2, 2, 1, 2, 1, amp_damp(0.5))
for lab, a, b, c in zip(labels, m, e, ex): print("   amp.damp rerun %-16s MC %.5f +- %.5f exact %.5f (%.1f sigma)" % (lab, a, b, c, abs(a-c)/b))
