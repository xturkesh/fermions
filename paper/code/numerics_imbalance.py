"""
Numerics for Section sec:numerics, Figure fig:imbalance (CDW imbalance of the Z2 LGT).

Conventions (paper/BRIEF.md, Sec. 2): sites x = 1..N (array index x-1), J = 1,
    h_q = -J Delta_N + 2h diag(q)   (Dirichlet ends),   q i.i.d. uniform on {+1,-1}^N,   lambda = 2h/J.
Initial state Psi_0 = |up...up>_z (x) CDW, CDW: n_x = 1 iff x is even;  Gamma = diag((-1)^x).
    I(t) = (2/N) sum_x (-1)^x <Psi(t), n_x Psi(t)>,   I(0) = 1.
By the exact sector averaging (cor:sector) and the proof of thm:imbalance(i), for even N
    I(t) = E_q N^{-1} Tr(Gamma U_q(t)^* Gamma U_q(t)) = E_q N^{-1} sum_{a,b} |G_ab|^2 cos((E_a-E_b)t),
with G = V^T Gamma V (V: real eigenvectors of h_q), and the infinite-time average is
    Ibar = E_q N^{-1} sum_a a_a^2,   a_a = <phi_a, Gamma phi_a> = G_aa.
Everything is exact up to round-off for each sampled sector; the only error is the Monte Carlo error of
the sector average (standard error over i.i.d. sectors).

kappa_N of thm:imbalance (eq:main-kappa) is max over ALL intervals B of E_q sum_{x in B, y notin B} R_q(x,y),
with R_q(x,y) = sum_a phi_a(x)^2 phi_a(y)^2; since it is linear in R we average R over sectors first and
maximize over intervals with 2D prefix sums.  Its statistical error is estimated by batch means at the
maximizing interval.  The lower bounds of thm:imbalance(iii) (global and block bounds) and (iv) are then
evaluated with the measured kappa (see the section text for what this does and does not prove).

Usage:
    python3 numerics_imbalance.py traces          -> data/numerics_imbalance_traces.npz
    python3 numerics_imbalance.py ibar N [ns]     -> data/numerics_imbalance_ibar_N{N}.npz
Recommended: OPENBLAS_NUM_THREADS=2 (or 1).  Runtime: traces ~3 min, ibar N=400 ~6 min.
Nothing computed here is used in a proof.
"""
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import j0

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
J = 1.0
BASESEED = 20261122
LAMS_TRACES = (0.5, 1.0, 2.0, 4.0)
LAMS_IBAR = (0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.6, 8.0)


def P(*a):
    print(*a, flush=True)


def rng_for(tag, lam, N):
    return np.random.default_rng(np.random.SeedSequence([BASESEED, tag, int(round(1000 * lam)), N]))


def gamma_diag(N):
    return np.array([(-1.0) ** x for x in range(1, N + 1)])          # (-1)^x, x = 1..N


def eig_q(q, h):
    N = len(q)
    return eigh_tridiagonal(2.0 * h * q, -J * np.ones(N - 1))


# ----------------------------------------------------------------------------------------- traces
def imbalance_trace(E, V, gam, times, iu):
    """I_q(t) for one sector (exact), and Ibar_q."""
    N = len(E)
    G = V.T @ (gam[:, None] * V)
    a = np.diag(G)
    ibq = float((a ** 2).sum() / N)
    w = 2.0 * G[iu] ** 2 / N
    om = E[iu[1]] - E[iu[0]]
    return ibq + np.cos(np.outer(times, om)) @ w, ibq


def run_traces(N=200, ns=1000):
    times = np.logspace(-1, 6, 211)
    gam = gamma_diag(N)
    iu = np.triu_indices(N, 1)
    P(f"=== I(t): N={N}, J={J}, CDW on even sites, {ns} i.i.d. sectors per lambda, "
      f"{len(times)} times in [{times[0]:.1e}, {times[-1]:.1e}] ===")
    res = dict(times=times, N=N, ns=ns, lams=np.array(LAMS_TRACES))
    mean = np.zeros((len(LAMS_TRACES), len(times)))
    se = np.zeros_like(mean)
    ibm = np.zeros(len(LAMS_TRACES))
    ibse = np.zeros(len(LAMS_TRACES))
    for il, lam in enumerate(LAMS_TRACES):
        h = lam * J / 2
        rng = rng_for(1, lam, N)
        s1 = np.zeros(len(times))
        s2 = np.zeros(len(times))
        ib = np.zeros(ns)
        err0 = 0.0
        t0 = time.time()
        for k in range(ns):
            q = rng.choice([-1.0, 1.0], size=N)
            E, V = eig_q(q, h)
            Iq, ib[k] = imbalance_trace(E, V, gam, times, iu)
            if k < 5:   # I_q(0) = 1 and direct evaluation at one time
                I0, _ = imbalance_trace(E, V, gam, np.array([0.0]), iu)
                U = (V * np.exp(-1j * E * times[50])) @ V.T
                nx = (np.abs(U[:, 1::2]) ** 2).sum(1)                  # occupied: x even <-> index odd
                Idir = (2.0 / N) * (gam * nx).sum()
                err0 = max(err0, abs(I0[0] - 1.0), abs(Idir - Iq[50]))
            s1 += Iq
            s2 += Iq ** 2
        mean[il] = s1 / ns
        se[il] = np.sqrt(np.maximum(s2 / ns - mean[il] ** 2, 0.0) / (ns - 1))
        ibm[il] = ib.mean()
        ibse[il] = ib.std(ddof=1) / np.sqrt(ns)
        late = times >= 1e4
        dev = (mean[il][late] - ibm[il]) / np.sqrt(se[il][late] ** 2 + ibse[il] ** 2)
        P(f"  lambda={lam:4.1f}: Ibar = {ibm[il]:.5f} +- {ibse[il]:.5f};  I(t) at t=1,10,1e2,1e3,1e6: "
          + ", ".join(f"{mean[il][np.argmin(abs(times - tt))]:.4f}+-{se[il][np.argmin(abs(times - tt))]:.4f}"
                      for tt in (1.0, 10.0, 1e2, 1e3, 1e6))
          + f";  t>=1e4: max |I(t)-Ibar|/err = {np.max(np.abs(dev)):.2f};  self-check err = {err0:.1e}"
          + f"  [{time.time() - t0:.0f}s]")
    res.update(mean=mean, se=se, ibar=ibm, ibar_se=ibse)
    # clean sector q = (1,...,1): the potential is constant, so this is the h = 0 chain
    E, V = eig_q(np.ones(N), 0.0)
    Ic, ibc = imbalance_trace(E, V, gam, times, iu)
    res.update(clean=Ic, clean_ibar=ibc, j0=j0(4 * J * times))
    k = times <= 40
    P(f"  clean sector q=1: Ibar = {ibc:.1e} (exactly 0, rem:main-imbalance(ii)); "
      f"max_(t<=40) |I_clean(t) - J_0(4Jt)| = {np.max(np.abs(Ic[k] - j0(4 * J * times[k]))):.1e}; "
      f"max_(t>=1e3) |I_clean(t)| = {np.max(np.abs(Ic[times >= 1e3])):.3f}")
    np.savez_compressed(os.path.join(DATA, "numerics_imbalance_traces.npz"), **res)


# ----------------------------------------------------------------------------------------- Ibar, kappa
def interval_cross(R, a, b):
    """sum_{x in [a,b], y notin [a,b]} R(x,y) (0-based, inclusive)."""
    blk = R[a:b + 1]
    return float(blk.sum() - blk[:, a:b + 1].sum())


def kappa_all_intervals(R):
    """max over intervals B=[a,b] of |B| - sum_{x,y in B} R(x,y) (= crossing sum, since rows of R sum to 1)."""
    N = R.shape[0]
    C = np.zeros((N + 1, N + 1))
    C[1:, 1:] = np.cumsum(np.cumsum(R, 0), 1)
    A = np.arange(N)[:, None]
    B = np.arange(N)[None, :]
    inner = C[B + 1, B + 1] - C[A, B + 1] - C[B + 1, A] + C[A, A]
    val = np.where(A <= B, (B - A + 1) - inner, -np.inf)
    i = int(np.argmax(val))
    a, b = divmod(i, N)
    return float(val[a, b]), int(a), int(b)


def imbalance_bounds(N, h, kap):
    """thm:imbalance(iii): global bound and best block bound (max over 1 <= L <= N)."""
    c = h * h / (J + abs(h)) ** 2
    glob = c / N
    L = np.arange(1, N + 1)
    blk = (N // L) / N * c * np.clip(np.sqrt(L / 2.0) - 2.0 * kap, 0.0, None) ** 2 / L
    i = int(np.argmax(blk))
    return glob, float(blk[i]), int(L[i])


def uniform_bound(h, kap):
    """thm:imbalance(iv) with kappa := kap."""
    Ls = max(2, int(np.ceil(32.0 * kap ** 2)))
    return h * h / (16.0 * Ls * (J + abs(h)) ** 2), Ls


def run_ibar(N, ns=2000, K=10):
    assert N % 2 == 0 and ns % K == 0
    gam = gamma_diag(N)
    P(f"=== Ibar, kappa_N: N={N}, {ns} i.i.d. sectors per lambda ({K} batches for the error of kappa_N) ===")
    nl = len(LAMS_IBAR)
    out = dict(N=N, ns=ns, lams=np.array(LAMS_IBAR), ibar=np.zeros(nl), ibar_se=np.zeros(nl),
               kappa=np.zeros(nl), kappa_se=np.zeros(nl), kap_a=np.zeros(nl, int), kap_b=np.zeros(nl, int),
               Femp=np.zeros((nl, N)), S1emp=np.zeros(nl))
    for il, lam in enumerate(LAMS_IBAR):
        h = lam * J / 2
        rng = rng_for(2, lam, N)
        t0 = time.time()
        ib = np.zeros(ns)
        Rb = np.zeros((K, N, N))
        Qs = np.zeros((N, N))
        for k in range(ns):
            q = rng.choice([-1.0, 1.0], size=N)
            E, V = eig_q(q, h)
            Pm = V ** 2
            a = gam @ Pm
            ib[k] = (a ** 2).sum() / N
            Rb[k % K] += Pm @ Pm.T
            Av = np.abs(V)
            Qs += Av @ Av.T
        Rb /= ns // K
        R = Rb.mean(0)
        kap, a0, b0 = kappa_all_intervals(R)
        kb = np.array([interval_cross(Rb[j], a0, b0) for j in range(K)])
        Qm = Qs / ns
        Fd = np.array([np.max(np.diagonal(Qm, d)) for d in range(N)])      # max_x E Q(x, x+d)
        F = np.minimum(1.0, np.maximum.accumulate(Fd[::-1])[::-1])          # non-increasing envelope
        d = np.arange(N)
        out["ibar"][il] = ib.mean()
        out["ibar_se"][il] = ib.std(ddof=1) / np.sqrt(ns)
        out["kappa"][il] = kap
        out["kappa_se"][il] = kb.std(ddof=1) / np.sqrt(K)
        out["kap_a"][il], out["kap_b"][il] = a0 + 1, b0 + 1                   # 1-based sites
        out["Femp"][il] = F
        out["S1emp"][il] = float((d[1:] * F[1:]).sum())
        P(f"  lambda={lam:4.2f}: Ibar = {ib.mean():.5f} +- {out['ibar_se'][il]:.5f}, "
          f"kappa_N = {kap:.3f} +- {out['kappa_se'][il]:.3f} (B=[{a0+1},{b0+1}]), "
          f"Sigma_1(F_emp) = {out['S1emp'][il]:.3g} (F_emp(N-1) = {F[-1]:.1e})  [{time.time() - t0:.0f}s]")
    np.savez_compressed(os.path.join(DATA, f"numerics_imbalance_ibar_N{N}.npz"), **out)


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    what = sys.argv[1] if len(sys.argv) > 1 else "traces"
    if what == "traces":
        run_traces(ns=int(sys.argv[2]) if len(sys.argv) > 2 else 1000)
    elif what == "ibar":
        run_ibar(int(sys.argv[2]), ns=int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
    else:
        raise SystemExit("usage: numerics_imbalance.py traces [ns] | ibar N [ns]")
