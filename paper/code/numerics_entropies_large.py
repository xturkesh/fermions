"""
Numerics for Section sec:numerics, Figure fig:entropies(b,c): sector-averaged conditional entanglement S_G(t)
and number entropy S_N(t) of the Z2 LGT for N = 64, 128, 256 (free-fermion numerics via the exact duality).

Conventions (paper/BRIEF.md, Sec. 2): sites x = 1..N, J = 1, h_q = -J Delta_N + 2h diag(q), lambda = 2h/J,
Psi_0 = |up...up>_z (x) CDW (n_x = 1 iff x even), A = [1, l] with l = N/2.
For each sampled sector q (i.i.d. uniform; exact sector averaging, cor:sector), psi_q(t) is the Slater
determinant with correlation matrix K_q(t) = U_q(t) P_S U_q(t)^*.  With C_A = P_A K_q(t) P_A and its
eigenvalues nu_i (lem:slater):
    S(rho^q_A(t)) = sum_i h(nu_i),    p^q_t = Poisson-binomial law of N_A with parameters nu_i,
    S_G(t) = E_q S(rho^q_A(t)),       S_N(t) = H(p_t),  p_t = E_q p^q_t   (prop:entropies(c)).
U_q(t) = V exp(-itE) V^T is exact for every t (no time stepping).  Errors: standard error over sectors for
S_G and E_q H(p^q_t); delta method for S_N = H(E_q p^q_t) (standard error of -sum_k p^q_t(k) ln p_t(k)),
cross-checked against a bootstrap for one (lambda, N).

Usage:  python3 numerics_entropies_large.py LAMBDA [ns]   -> data/numerics_entropies_large_lam{LAMBDA}.npz
Recommended OPENBLAS_NUM_THREADS=2.  Runtime ~4-6 min per lambda with ns = 1000.  Not used in any proof.
"""
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
J = 1.0
BASESEED = 20261122
NS_LIST = (64, 128, 256)


def P(*a):
    print(*a, flush=True)


def hbin(nu):
    nu = np.clip(nu, 1e-300, 1.0 - 1e-16)
    return -nu * np.log(nu) - (1.0 - nu) * np.log1p(-nu)


def poibin_batch(nu):
    """nu: (nt, l) -> Poisson-binomial pmf (nt, l+1)."""
    nt, l = nu.shape
    p = np.zeros((nt, l + 1))
    p[:, 0] = 1.0
    for i in range(l):
        v = nu[:, i:i + 1]
        p[:, 1:i + 2] = p[:, 1:i + 2] * (1.0 - v) + p[:, 0:i + 1] * v
        p[:, 0] *= (1.0 - v[:, 0])
    return p


def xlogx_entropy(p, axis=-1):
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(p > 1e-300, -p * np.log(np.where(p > 1e-300, p, 1.0)), 0.0)
    return t.sum(axis)


def run(lam, N, ns, times):
    h = lam * J / 2
    l = N // 2
    S = np.arange(1, N, 2)                       # occupied: x even <-> 0-based index odd
    rng = np.random.default_rng(np.random.SeedSequence([BASESEED, 3, int(round(1000 * lam)), N]))
    nt = len(times)
    SGq = np.zeros((ns, nt))
    EHq = np.zeros((ns, nt))
    Vq = np.zeros((ns, nt))
    pq = np.zeros((ns, nt, l + 1))
    t0 = time.time()
    for k in range(ns):
        q = rng.choice([-1.0, 1.0], size=N)
        E, V = eigh_tridiagonal(2.0 * h * q, -J * np.ones(N - 1))
        VA = V[:l]
        VS = V[S]
        ph = np.exp(-1j * np.outer(times, E))                         # (nt, N)
        W = (VA[None, :, :] * ph[:, None, :]) @ VS.T                  # U_q(t)[A, S]
        C = W @ np.conj(np.transpose(W, (0, 2, 1)))                    # C_A = P_A K P_A
        nu = np.clip(np.linalg.eigvalsh(C), 0.0, 1.0)
        SGq[k] = hbin(nu).sum(1)
        Vq[k] = (nu * (1.0 - nu)).sum(1)
        p = poibin_batch(nu)
        pq[k] = p
        EHq[k] = xlogx_entropy(p)
    pbar = pq.mean(0)                                                  # (nt, l+1)
    SN = xlogx_entropy(pbar)
    with np.errstate(divide="ignore"):
        logp = np.where(pbar > 1e-300, np.log(np.where(pbar > 1e-300, pbar, 1.0)), 0.0)
    g = -(pq * logp[None]).sum(2)                                      # (ns, nt): linearization of H at pbar
    SN_se = g.std(0, ddof=1) / np.sqrt(ns)
    kk = np.arange(l + 1)
    mean = (pbar * kk).sum(1)
    var_mix = (pbar * kk ** 2).sum(1) - mean ** 2
    late = times >= 1e3                                                # plateau window
    sgp = SGq[:, late].mean(1)
    gp = g[:, late].mean(1)
    out = dict(SG=SGq.mean(0), SG_se=SGq.std(0, ddof=1) / np.sqrt(ns),
               SG_plateau=sgp.mean(), SG_plateau_se=sgp.std(ddof=1) / np.sqrt(ns),
               SN_plateau=SN[late].mean(), SN_plateau_se=gp.std(ddof=1) / np.sqrt(ns),
               EHq=EHq.mean(0), EHq_se=EHq.std(0, ddof=1) / np.sqrt(ns),
               SN=SN, SN_se=SN_se, var_mix=var_mix, EVq=Vq.mean(0), mean_NA=mean,
               DNG_se=(g - SGq).std(0, ddof=1) / np.sqrt(ns),                 # error of S_N - S_G (delta method)
               SGsup=SGq.max(1).mean(), SGsup_se=SGq.max(1).std(ddof=1) / np.sqrt(ns))
    boot = None
    if N == 64:   # bootstrap check of the delta-method error of S_N
        br = np.random.default_rng(7)
        bs = np.array([xlogx_entropy(pq[br.integers(0, ns, ns)].mean(0)) for _ in range(200)])
        boot = bs.std(0, ddof=1)
        out["SN_boot_se"] = boot
    P(f"  N={N:3d}: {ns} sectors, {time.time() - t0:.0f}s")
    return out


def main():
    lam = float(sys.argv[1])
    ns = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    times = np.concatenate([[0.0], np.logspace(-2, 8, 61)])
    P(f"=== S_G(t), S_N(t): lambda={lam}, J={J}, CDW on even sites, A=[1,N/2], N in {NS_LIST}, "
      f"{len(times)} times in [0] u [1e-2, 1e8] ===")
    res = dict(times=times, lam=lam, ns=ns, Ns=np.array(NS_LIST))
    for N in NS_LIST:
        o = run(lam, N, ns, times)
        for key, val in o.items():
            res[f"{key}_N{N}"] = val
        late = times >= 1e3
        P(f"    S_G : t=1: {o['SG'][np.argmin(abs(times-1))]:.3f}+-{o['SG_se'][np.argmin(abs(times-1))]:.3f}, "
          f"t=10: {o['SG'][np.argmin(abs(times-10))]:.3f}, plateau mean(t>=1e3) = {o['SG'][late].mean():.3f} "
          f"(typical err {o['SG_se'][late].mean():.3f}); range over t>=1e3: [{o['SG'][late].min():.3f}, {o['SG'][late].max():.3f}]")
        P(f"    S_N : t=1: {o['SN'][np.argmin(abs(times-1))]:.3f}+-{o['SN_se'][np.argmin(abs(times-1))]:.3f}, "
          f"t=10: {o['SN'][np.argmin(abs(times-10))]:.3f}, plateau mean(t>=1e3) = {o['SN'][late].mean():.3f} "
          f"(typical err {o['SN_se'][late].mean():.3f}); range over t>=1e3: [{o['SN'][late].min():.3f}, {o['SN'][late].max():.3f}]")
        P(f"    E_qH(p^q) plateau = {o['EHq'][late].mean():.3f};  Var_mix(N_A) plateau = {o['var_mix'][late].mean():.3f}, "
          f"max_t = {o['var_mix'].max():.3f};  E_q Var_q plateau = {o['EVq'][late].mean():.3f};  "
          f"E_q max_(t in grid) S(rho^q_A) = {o['SGsup']:.3f}+-{o['SGsup_se']:.3f}")
        d = o["SN"] - o["SG"]
        pos = times[(d > 0) & (times > 0)]
        P(f"    S_N - S_G > 0 on grid points: " + (f"t in [{pos.min():.3g}, {pos.max():.3g}] ({len(pos)} points)" if len(pos) else "none")
          + f"; max(S_N - S_G) = {d.max():.3f} at t = {times[np.argmax(d)]:.3g}")
        if "SN_boot_se" in o:
            P(f"    delta-method vs bootstrap error of S_N (t>0): max ratio = {np.max(o['SN_se'][1:] / o['SN_boot_se'][1:]):.2f}, "
              f"min ratio = {np.min(o['SN_se'][1:] / o['SN_boot_se'][1:]):.2f}")
    # N-independence of the plateau
    late = times >= 1e3
    for key in ("SG", "SN"):
        vals = [(float(res[f"{key}_plateau_N{N}"]), float(res[f"{key}_plateau_se_N{N}"])) for N in NS_LIST]
        P(f"  {key} plateau (time average over the grid points t>=1e3, error over sectors) for N=64,128,256: "
          + ", ".join(f"{v:.4f}+-{e:.4f}" for v, e in vals))
    np.savez_compressed(os.path.join(DATA, f"numerics_entropies_large_lam{lam:g}.npz"), **res)


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    main()
