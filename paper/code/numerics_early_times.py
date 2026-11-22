"""Early-time behaviour of S_N - S_G (Section sec:numerics, paragraph on S_G versus S_N).

Two checks, both quoted in sections/numerics.tex:
  (1) the q-dependence of Var_q(N_A) at small t:
      Var_q(N_A) = J^2 t^2 - 2 J^4 t^4 - (1/3) J^2 h^2 (q_l - q_{l+1})^2 t^4 + O(t^6),
      so that equal and opposite charges at the cut differ by -(4/3) J^2 h^2 t^4 (here J = 1);
  (2) the comparison of the stored data for S_N - S_G with the leading-order Jensen gap lambda^4 (Jt)^6 / 72,
      for N = 10 (exact, all 1024 sectors) on 0.01 <= Jt <= 0.1 and for N = 256 (400 sectors) on 0.02 <= Jt <= 0.1,
      both at lambda = 2.  The remainder of the expansion is not controlled; this is a numerical comparison only.
Adapted at assembly from the referee scripts for the numerics section.  Uses data/numerics_entropies_small.npz and
data/numerics_entropies_early_lam2.npz produced by numerics_entropies_small.py and numerics_entropies_large.py.
Run:  python3 numerics_early_times.py > numerics_early_times.out
"""
import os
import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def check_variance_expansion():
    rng = np.random.default_rng(1)
    print("=== (1) Var_q(N_A) at Jt = 0.01: equal versus opposite charges at the cut (J = 1, CDW on even sites) ===")
    for N, ell in ((256, 128), (10, 5)):
        S = np.arange(1, N, 2)          # occupied sites x = 2, 4, ... (0-based indices 1, 3, ...)
        A = np.arange(ell)
        B = np.arange(ell, N)
        for h in (1.0, 0.5):
            out = {}
            for pair in ((1, 1), (1, -1)):
                vals = []
                for _ in range(20):
                    q = rng.choice([-1.0, 1.0], N)
                    q[ell - 1], q[ell] = pair
                    hq = -(np.eye(N, k=1) + np.eye(N, k=-1)) + 2 * h * np.diag(q)
                    t = 0.01
                    U = expm(-1j * t * hq)
                    K = U[:, S] @ U[:, S].conj().T
                    vals.append((np.abs(K[np.ix_(B, A)]) ** 2).sum())    # Var_q(N_A) = ||P_B K P_A||_2^2
                out[pair] = (np.mean(vals), np.std(vals))
            d = out[(1, -1)][0] - out[(1, 1)][0]
            print(f"  N={N} h={h}: equal charges {out[(1, 1)][0]:.10e} (spread over samples {out[(1, 1)][1]:.1e}), "
                  f"opposite {out[(1, -1)][0]:.10e}; difference {d:.4e}, predicted -(4/3) h^2 t^4 = {-(4/3) * h**2 * 1e-8:.4e}")


def compare_with_leading_order():
    lam = 2.0
    print("=== (2) S_N - S_G versus lambda^4 (Jt)^6 / 72 at lambda = 2 ===")
    sm = np.load(os.path.join(DATA, "numerics_entropies_small.npz"))
    t, d = sm["times"], sm["SN"] - sm["SG"]
    w = (t >= 0.0099) & (t <= 0.1001)
    r = d[w] / (lam**4 * t[w] ** 6 / 72)
    print(f"  N=10 (exact): grid times {np.round(t[w], 4).tolist()}")
    print(f"     ratio data/formula: {np.round(r, 3).tolist()}; max |ratio - 1| = {np.abs(r - 1).max():.3f}")
    ea = np.load(os.path.join(DATA, "numerics_entropies_early_lam2.npz"))
    t, d = ea["times"], ea["SN_N256"] - ea["SG_N256"]
    w = (t >= 0.0199) & (t <= 0.1001)
    r = d[w] / (lam**4 * t[w] ** 6 / 72)
    print(f"  N=256 (400 sectors): grid times {np.round(t[w], 4).tolist()}")
    print(f"     ratio data/formula: {np.round(r, 3).tolist()}; max |ratio - 1| = {np.abs(r - 1).max():.3f}")


if __name__ == "__main__":
    check_variance_expansion()
    compare_with_leading_order()
