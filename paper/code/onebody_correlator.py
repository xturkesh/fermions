"""
Numerical illustration for Section sec:onebody (Theorem thm:EC, Figure fig:correlator).

Estimates the sector-averaged eigenfunction correlator
    E_q Q^{[1,N]}(x, x+d),   Q(x,y) = sum_alpha |phi_alpha(x) phi_alpha(y)|,
for the one-body operator h_q / J = -Delta_N + lambda diag(q), q i.i.d. uniform on {-1,+1},
at the boundary (x = 1) and in the bulk (x = N/2), for lambda in {1, 2, 4} and N in {100, 200, 400}.
(Q is invariant under h -> J h and under the sign change Delta -> -Delta, so we diagonalise
-Delta + lambda q directly.)  Also estimates the Lyapunov exponent gamma_lambda(E) on a grid
of energies (transfer-matrix products), for comparison with the fitted annealed decay rates.

Output: ../figs/correlator.pdf and onebody_correlator.out (next to this script).
Usage:  python3 onebody_correlator.py [samples]   (default 4000 samples per (lambda, N))
"""
import os
import sys
import numpy as np
from multiprocessing import Pool
from scipy.linalg import eigh_tridiagonal

HERE = os.path.dirname(os.path.abspath(__file__))
LAMBDAS = [1.0, 2.0, 4.0]
NS = [100, 200, 400]
DMAX = 80
FIT = (10, 40)          # fit window for the decay rate (distances d)


def worker(args):
    lam, N, nsamp, seed = args
    rng = np.random.default_rng(seed)
    dmax_b = min(DMAX, N - 1)
    xb = N // 2 - 1                      # 0-based index of the bulk site x = N/2
    dmax_c = min(DMAX, N - 1 - xb)
    sb = np.zeros(dmax_b + 1); sb2 = np.zeros(dmax_b + 1)
    sc = np.zeros(dmax_c + 1); sc2 = np.zeros(dmax_c + 1)
    for _ in range(nsamp):
        q = rng.choice([-1.0, 1.0], size=N)
        _, V = eigh_tridiagonal(lam * q, -np.ones(N - 1))
        A = np.abs(V)                     # A[x, alpha] = |phi_alpha(x)|
        Qb = A @ A[0]                     # Q(1, .)
        Qc = A @ A[xb]                    # Q(N/2, .)
        b = Qb[: dmax_b + 1]; c = Qc[xb: xb + dmax_c + 1]
        sb += b; sb2 += b * b; sc += c; sc2 += c * c
    return lam, N, nsamp, sb, sb2, sc, sc2


def lyapunov(lam, energies, n=200000, seed=7):
    """gamma(E) = lim (1/n) log ||T_n(E)||, T(E,v) = [[E - v, -1], [1, 0]], v = lambda q."""
    rng = np.random.default_rng(seed)
    v = lam * rng.choice([-1.0, 1.0], size=n)
    u = np.ones_like(energies); w = np.zeros_like(energies)   # (u_k, u_{k-1})
    acc = np.zeros_like(energies)
    for k in range(n):
        u, w = (energies - v[k]) * u - w, u
        r = np.sqrt(u * u + w * w)
        acc += np.log(r); u /= r; w /= r
    return acc / n


def main():
    nsamp = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    chunks = 8
    jobs = []
    seed = 12345
    for lam in LAMBDAS:
        for N in NS:
            for c in range(chunks):
                seed += 1
                jobs.append((lam, N, nsamp // chunks, seed))
    with Pool(min(4, os.cpu_count() or 1)) as pool:
        res = pool.map(worker, jobs)
    data = {}
    for lam, N, ns, sb, sb2, sc, sc2 in res:
        key = (lam, N)
        if key not in data:
            data[key] = [0, 0, 0, 0, 0]
        d = data[key]
        d[0] += ns; d[1] = d[1] + sb; d[2] = d[2] + sb2; d[3] = d[3] + sc; d[4] = d[4] + sc2
    stats = {}
    for key, (ns, sb, sb2, sc, sc2) in data.items():
        mb = sb / ns; mc = sc / ns
        eb = np.sqrt(np.maximum(sb2 / ns - mb ** 2, 0) / ns)
        ec = np.sqrt(np.maximum(sc2 / ns - mc ** 2, 0) / ns)
        stats[key] = (ns, mb, eb, mc, ec)

    out = open(os.path.join(HERE, 'onebody_correlator.out'), 'w')
    def pr(s):
        print(s); out.write(s + '\n')
    pr(f'# E_q Q^[1,N](x, x+d), boundary x=1 and bulk x=N/2; samples per (lambda,N) = {nsamp}')
    pr(f'# rate = -slope of log E Q over d in [{FIT[0]},{FIT[1]}] (least squares); +- = 1 standard error of the mean')
    rates = {}
    for lam in LAMBDAS:
        for N in NS:
            ns, mb, eb, mc, ec = stats[(lam, N)]
            dd = np.arange(FIT[0], FIT[1] + 1)
            rb = -np.polyfit(dd, np.log(mb[dd]), 1)[0]
            rc = -np.polyfit(dd, np.log(mc[dd]), 1)[0]
            rates[(lam, N)] = (rb, rc)
            sel = [1, 10, 20, 40]
            pr(f'lambda={lam:3.1f} N={N:3d} | boundary E Q(d=1,10,20,40) = ' +
               ' '.join(f'{mb[k]:.3e}+-{eb[k]:.1e}' for k in sel) + f'  rate={rb:.3f}')
            pr(f'                 |  bulk    E Q(d=1,10,20,40) = ' +
               ' '.join(f'{mc[k]:.3e}+-{ec[k]:.1e}' for k in sel) + f'  rate={rc:.3f}')
    pr('# Lyapunov exponent gamma_lambda(E) on a grid of E in [-2-|lambda|, 2+|lambda|] (n = 2e5 steps)')
    gam = {}
    for lam in LAMBDAS:
        Es = np.linspace(-2 - lam, 2 + lam, 161)
        g = lyapunov(lam, Es)
        gam[lam] = (Es, g)
        pr(f'lambda={lam:3.1f}: min_E gamma = {g.min():.3f} (at E={Es[g.argmin()]:+.2f}), gamma(0) = {g[80]:.3f}, '
           f'max_E gamma = {g.max():.3f}')
    out.close()

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                         'ytick.major.width': 0.6, 'legend.frameon': False})
    colors = {1.0: '#2a78d6', 2.0: '#eb6834', 4.0: '#1baf7a'}   # validated categorical slots 1-3
    grid_c = '#d9d8d4'
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), constrained_layout=True)

    ax = axs[0]
    N = 400
    for lam in LAMBDAS:
        ns, mb, eb, mc, ec = stats[(lam, N)]
        d = np.arange(len(mb)); dc = np.arange(len(mc))
        ax.semilogy(d[1:], mb[1:], '-', color=colors[lam], lw=1.6, label=fr'$\lambda={lam:g}$, $x=1$')
        ax.semilogy(dc[1:], mc[1:], '--', color=colors[lam], lw=1.6, label=fr'$\lambda={lam:g}$, $x=N/2$')
    ax.set_xlabel(r'distance $d=|x-y|$')
    ax.set_ylabel(r'$\mathbb{E}_q\, Q^{[1,N]}(x,x+d)$')
    ax.set_title(r'(a) boundary vs bulk, $N=400$', fontsize=9)
    ax.set_xlim(0, DMAX); ax.set_ylim(1e-7, 1.2)
    ax.grid(True, which='major', color=grid_c, lw=0.5)
    ax.legend(fontsize=7, ncol=1, loc='lower left')

    ax = axs[1]
    lam = 2.0
    styles = {100: ':', 200: '--', 400: '-'}
    for N in NS:
        ns, mb, eb, mc, ec = stats[(lam, N)]
        d = np.arange(len(mb))
        ax.semilogy(d[1:], mb[1:], styles[N], color=colors[lam], lw=1.6, label=fr'$N={N}$, $x=1$')
        ax.fill_between(d[1:], np.maximum(mb[1:] - 2 * eb[1:], 1e-12), mb[1:] + 2 * eb[1:],
                        color=colors[lam], alpha=0.15, lw=0)
    ns, mb, eb, mc, ec = stats[(lam, 400)]
    dc = np.arange(len(mc))
    ax.semilogy(dc[1:], mc[1:], '-', color='#52514e', lw=1.0, label=r'$N=400$, $x=N/2$')
    ax.set_xlabel(r'distance $d=|x-y|$')
    ax.set_title(r'(b) $N$-independence at the boundary, $\lambda=2$', fontsize=9)
    ax.set_xlim(0, DMAX); ax.set_ylim(1e-7, 1.2)
    ax.grid(True, which='major', color=grid_c, lw=0.5)
    ax.legend(fontsize=7, loc='lower left')
    figdir = os.path.join(os.path.dirname(HERE), 'figs')
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, 'correlator.pdf'))


if __name__ == '__main__':
    main()
