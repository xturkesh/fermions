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
    python3 numerics_imbalance.py spectral        -> printed only (heuristic check of the slow oscillations)
Recommended: OPENBLAS_NUM_THREADS=1 or 2.  Runtime: traces ~8 min, ibar N=400 ~5 min.
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


def run_traces(N=200, ns=1000, K=10, dt=0.05, tlin=200.0, M=2 ** 23):
    """I(t) on a uniform grid t = 0, dt, ..., tlin (all oscillations resolved) and on a logarithmic grid
    tlin < t <= 1e6.  Uniform grid: the frequencies omega_ab = E_b - E_a >= 0 and weights 2 G_ab^2 / N of all
    sampled sectors are deposited on a frequency grid of spacing dw = 2 pi / (M dt) by linear interpolation,
    and sum_b W_b cos(omega_b t_n) is evaluated by one real FFT per batch.  For |t| <= tlin the interpolation
    error of cos(omega t) is at most (dw t)^2 / 8 per unit weight (< 2e-6 here); it is also measured below.
    Logarithmic grid: exact evaluation.  Errors: batch means (K batches) on the uniform grid, standard
    error over sectors on the logarithmic grid."""
    assert ns % K == 0
    nd = int(round(tlin / dt)) + 1
    tlin_grid = dt * np.arange(nd)
    tlog = np.logspace(np.log10(tlin), 6, 187)[1:]
    gam = gamma_diag(N)
    iu = np.triu_indices(N, 1)
    dw = 2 * np.pi / (M * dt)
    P(f"=== I(t): N={N}, J={J}, CDW on even sites, {ns} i.i.d. sectors per lambda; uniform grid dt={dt} on "
      f"[0,{tlin:g}] ({nd} times, FFT with dw={dw:.2e}), {len(tlog)} log-spaced times in ({tlin:g}, 1e6] ===")
    nl = len(LAMS_TRACES)
    res = dict(tlin=tlin_grid, tlog=tlog, N=N, ns=ns, lams=np.array(LAMS_TRACES),
               mean_lin=np.zeros((nl, nd)), se_lin=np.zeros((nl, nd)),
               mean_log=np.zeros((nl, len(tlog))), se_log=np.zeros((nl, len(tlog))),
               ibar=np.zeros(nl), ibar_se=np.zeros(nl))
    for il, lam in enumerate(LAMS_TRACES):
        h = lam * J / 2
        nb = int(np.ceil((4 * J + 4 * abs(h)) / dw)) + 3
        rng = rng_for(1, lam, N)
        Wb = np.zeros((K, nb))
        ibb = np.zeros(K)
        s1 = np.zeros(len(tlog))
        s2 = np.zeros(len(tlog))
        ib = np.zeros(ns)
        err0 = 0.0
        errfft = 0.0
        t0 = time.time()
        for k in range(ns):
            q = rng.choice([-1.0, 1.0], size=N)
            E, V = eig_q(q, h)
            G = V.T @ (gam[:, None] * V)
            a = np.diag(G)
            ib[k] = (a ** 2).sum() / N
            w = 2.0 * G[iu] ** 2 / N
            om = E[iu[1]] - E[iu[0]]
            Iq = ib[k] + np.cos(np.outer(tlog, om)) @ w
            s1 += Iq
            s2 += Iq ** 2
            u = om / dw
            i0 = np.floor(u).astype(np.int64)
            f = u - i0
            Wk = np.bincount(i0, weights=w * (1.0 - f), minlength=nb) + np.bincount(i0 + 1, weights=w * f, minlength=nb)
            Wb[k % K] += Wk[:nb]
            ibb[k % K] += ib[k]
            if k < 3:   # self-checks: I_q(0) = 1, direct evaluation, FFT reconstruction of one sector
                tchk = np.array([0.0, 1.3, 37.0, tlin])
                Iex = ib[k] + np.cos(np.outer(tchk, om)) @ w
                U = (V * np.exp(-1j * E * 37.0)) @ V.T
                nx = (np.abs(U[:, 1::2]) ** 2).sum(1)                      # occupied: x even <-> index odd
                err0 = max(err0, abs(Iex[0] - 1.0), abs((2.0 / N) * (gam * nx).sum() - Iex[2]))
                Wpad = np.zeros(M)
                Wpad[:nb] = Wk[:nb]
                Ifft = ib[k] + np.fft.rfft(Wpad).real[:nd]
                idx = np.round(tchk / dt).astype(int)
                errfft = max(errfft, np.max(np.abs(Ifft[idx] - Iex)))
        nk = ns // K
        Ib_batches = np.zeros((K, nd))
        for j in range(K):
            Wpad = np.zeros(M)
            Wpad[:nb] = Wb[j]
            Ib_batches[j] = (ibb[j] + np.fft.rfft(Wpad).real[:nd]) / nk
        res["mean_lin"][il] = Ib_batches.mean(0)
        res["se_lin"][il] = Ib_batches.std(0, ddof=1) / np.sqrt(K)
        res["mean_log"][il] = s1 / ns
        res["se_log"][il] = np.sqrt(np.maximum(s2 / ns - (s1 / ns) ** 2, 0.0) / (ns - 1))
        res["ibar"][il] = ib.mean()
        res["ibar_se"][il] = ib.std(ddof=1) / np.sqrt(ns)
        ml, sl = res["mean_lin"][il], res["se_lin"][il]
        mg, sg = res["mean_log"][il], res["se_log"][il]

        def at(tt):
            if tt <= tlin:
                j = int(round(tt / dt))
                return f"{ml[j]:.4f}({sl[j]:.4f})"
            j = int(np.argmin(np.abs(tlog - tt)))
            return f"{mg[j]:.4f}({sg[j]:.4f}) [t={tlog[j]:.3g}]"
        late = tlog >= 1e4
        dev = (mg[late] - res["ibar"][il]) / np.sqrt(sg[late] ** 2 + res["ibar_se"][il] ** 2)
        P(f"  lambda={lam:4.1f}: Ibar = {res['ibar'][il]:.5f}({res['ibar_se'][il]:.5f}); I(t) at t = 1, 10, 29, 100, "
          f"1e3, 1e4, 2.9e4, 1e5, 1e6: " + ", ".join(at(tt) for tt in (1.0, 10.0, 29.0, 100.0, 1e3, 1e4, 2.9e4, 1e5, 1e6)))
        P(f"      t>=1e4: max |I(t)-Ibar|/err = {np.max(np.abs(dev)):.2f};  max over t in [100,200] of |I(t)-Ibar| = "
          f"{np.max(np.abs(ml[tlin_grid >= 100] - res['ibar'][il])):.4f};  self-checks: |I(0)-1|, direct = {err0:.1e}, "
          f"FFT reconstruction error = {errfft:.1e}  [{time.time() - t0:.0f}s]")
    # clean sector q = (1,...,1): constant potential, i.e. the h = 0 chain (exact, t <= 40)
    E, V = eig_q(np.ones(N), 0.0)
    G = V.T @ (gam[:, None] * V)
    w = 2.0 * G[iu] ** 2 / N
    om = E[iu[1]] - E[iu[0]]
    tc = tlin_grid[tlin_grid <= 40.0]
    Ic = (np.diag(G) ** 2).sum() / N + np.cos(np.outer(tc, om)) @ w
    res.update(tclean=tc, clean=Ic, j0=j0(4 * J * tc))
    P(f"  clean sector q=1: Ibar = {(np.diag(G) ** 2).sum() / N:.1e} (exactly 0, rem:gaussproj(ii)); "
      f"max_(t<=40) |I_clean(t) - J_0(4Jt)| = {np.max(np.abs(Ic - j0(4 * J * tc))):.1e}")
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


def run_spectral(N=200, ns=400, lam=4.0):
    """spectral weight of the oscillating part of I_q(t): w(omega) = E_q N^{-1} sum_{a != b} G_ab^2 delta(omega - |E_a-E_b|),
    histogrammed in bins of 0.02 decades.  Heuristic check of the slow oscillations at strong coupling (see the text):
    effective k-th order hoppings J (2 lambda)^{-(k-1)} between equal charges at odd distance k, separated by k-1
    opposite charges, give oscillations at omega_k = 2 J (2 lambda)^{-(k-1)}."""
    gam = gamma_diag(N)
    iu = np.triu_indices(N, 1)
    h = lam * J / 2
    rng = rng_for(4, lam, N)
    bins = np.logspace(-6, 0, 301)
    W = np.zeros(len(bins) - 1)
    for _ in range(ns):
        q = rng.choice([-1.0, 1.0], size=N)
        E, V = eig_q(q, h)
        G = V.T @ (gam[:, None] * V)
        W += np.histogram(np.abs(E[iu[1]] - E[iu[0]]), bins=bins, weights=2.0 * G[iu] ** 2 / N)[0] / ns
    dens = W / 0.02
    P(f"=== spectral weight of I_q(t) - Ibar_q: N={N}, lambda={lam}, {ns} sectors, bins of 0.02 decades ===")
    for k, (lo, hi) in zip((7, None, 5, None, 3, None), ((1e-6, 1e-5), (1e-5, 1e-4), (1e-4, 1e-3), (1e-3, 1e-2),
                                                        (1e-2, 1e-1), (1e-1, 1.0))):
        m = (bins[:-1] >= lo * 0.999) & (bins[1:] <= hi * 1.001)
        i = int(np.argmax(dens[m]))
        c = np.sqrt(bins[:-1][m][i] * bins[1:][m][i])
        line = (f"  omega in [{lo:.0e},{hi:.0e}]: weight {W[m].sum():.4f}; peak density {dens[m].max():.4f}/decade at "
                f"omega = {c:.2e}; median density {np.median(dens[m]):.4f}/decade")
        if k is not None:
            om = 2.0 * J * (2.0 * lam) ** (-(k - 1))
            mk = (bins[:-1] >= om / 1.2) & (bins[1:] <= om * 1.2)
            line += f";  omega_{k} = {om:.2e} (period {2*np.pi/om:.2g}/J): weight within a factor 1.2 = {W[mk].sum():.4f}"
        P(line)


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    what = sys.argv[1] if len(sys.argv) > 1 else "traces"
    if what == "traces":
        run_traces(ns=int(sys.argv[2]) if len(sys.argv) > 2 else 1000)
    elif what == "spectral":
        run_spectral()
    elif what == "ibar":
        run_ibar(int(sys.argv[2]), ns=int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
    else:
        raise SystemExit("usage: numerics_imbalance.py traces [ns] | ibar N [ns] | spectral")
