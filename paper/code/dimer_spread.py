"""
Numerics for Section sec:dimer, Figure fig:dimer (figs/dimer.pdf).

Single dressed fermion in the Z2 LGT, started at site x0 (odd, i.e. the left site of a dimer).  By
Prop. prop:dimerlaw / Cor. cor:sector its density is the sector average of |<delta_y, e^{-i t h_q} delta_x0>|^2
with h_q = -J Delta + 2h diag(q), for
   law 'dimer':  q_{2k-1} = q_{2k} = eta_k, eta_k iid uniform (dimer initial state of Def. def:dimerstate),
   law 'iid'  :  q_j iid uniform (SKKM initial state |up...up>_z x psi, Cor. cor:sector).
J = 1, h/J in {0.25, 0.75}, i.e. lambda = 2h/J in {0.5, 1.5}.

The one-body dynamics is computed on the infinite chain Z: Chebyshev propagation (exact to ~1e-14 per step)
on an adaptive window that always contains every site where |psi| > 1e-18 for some sample, plus the light cone
of the next step, so no boundary is ever reached (the largest window has about 4 J t_max + O(t_max^{1/3}) sites).

Outputs
  - <(X-x0)^2>_t averaged over S sectors (+- standard error), for t in [0, 1e4];
  - the exponential time average M_2(T) = int_0^infty (dt/T) e^{-t/T} <(X-x0)^2>_t (the quantity of Thm thm:dimer),
    per sector and averaged, for T in [1, 1e3] (trapezoidal rule on the time grid; validated below against the
    exact eigenfunction formula on a finite chain);
  - figs/dimer.pdf, code/dimer_spread.out, code/dimer_spread_data.npz (mean/s.e. of <X^2>_t and
    M_2(T) per sector: everything the figure needs).
Usage: python3 dimer_spread.py          (compute + plot)
       python3 dimer_spread.py --plot   (plot from dimer_spread_data.npz)
"""
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import jv

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'dimer_spread_data.npz')
J = 1.0
HS = [0.25, 0.75]
LAWS = ['dimer', 'iid']
S = 64                       # sectors (samples) per case
TMAX = 1.0e4
TGRID = np.concatenate([[0.0], np.logspace(-2, 4, 601)])
TAVG = np.logspace(0, 3, 31)
CUT = 1e-18


def potential(law, lam, L, rng):
    """(L, S) array of 2h q_x / J = lam q_x; array index i <-> site x = i+1; dimers {2k-1,2k} <-> (2k-2,2k-1)."""
    if law == 'dimer':
        eta = rng.choice([-1.0, 1.0], size=((L + 1) // 2, S))
        return lam * np.repeat(eta, 2, axis=0)[:L]
    return lam * rng.choice([-1.0, 1.0], size=(L, S))


def cheb_step(V, psi, lo, hi, dt, r):
    """psi <- e^{-i h dt} psi with h = -J(shift+ + shift-) + V, on the window [lo-K, hi+K] (exact: the k-th
    Chebyshev iterate is supported in [lo-k, hi+k]).  Spectrum of h in [-r, r]."""
    a = r * dt
    K = int(a + 10.0 * a ** (1.0 / 3.0) + 30)
    while abs(jv(K, a)) > 1e-17:
        K += 10
    L = psi.shape[0]
    a0, b0 = max(lo - K, 0), min(hi + K, L - 1)
    if a0 == 0 or b0 == L - 1:
        raise RuntimeError('window reached the end of the array; enlarge L')
    Vw = V[a0:b0 + 1]
    w0 = psi[a0:b0 + 1].copy()

    def Ht(v):
        out = Vw * v
        out[:-1] -= J * v[1:]
        out[1:] -= J * v[:-1]
        return out * (1.0 / r)
    w1 = Ht(w0)
    coef = jv(np.arange(K + 1), a)
    res = coef[0] * w0 + 2.0 * (-1j) * coef[1] * w1
    ph = -1j
    for k in range(2, K + 1):
        w2 = 2.0 * Ht(w1) - w0
        ph *= -1j
        res += (2.0 * ph * coef[k]) * w2
        w0, w1 = w1, w2
    psi[a0:b0 + 1] = res
    big = np.nonzero(np.max(np.abs(res), axis=1) > CUT)[0]
    new_lo, new_hi = a0 + big[0], a0 + big[-1]
    psi[a0:new_lo] = 0.0
    psi[new_hi + 1:b0 + 1] = 0.0
    return new_lo, new_hi


def run_case(law, h, seed):
    lam = 2 * h / J
    rng = np.random.default_rng(seed)
    half = int(2 * J * TMAX + 1500)
    L = 2 * half + 2
    i0 = half if half % 2 == 0 else half + 1            # even index <-> odd site x0 = i0+1 (left site of a dimer)
    V = potential(law, lam, L, rng)
    psi = np.zeros((L, S), complex)
    psi[i0] = 1.0
    r = 2 * J + abs(lam) * J + 1e-9
    lo = hi = i0
    m2 = np.zeros((len(TGRID), S))
    norm_err = 0.0
    width = 0
    for it in range(1, len(TGRID)):
        lo, hi = cheb_step(V, psi, lo, hi, TGRID[it] - TGRID[it - 1], r)
        y = np.arange(lo, hi + 1) - i0
        p = np.abs(psi[lo:hi + 1]) ** 2
        m2[it] = (y[:, None] ** 2 * p).sum(0)
        norm_err = max(norm_err, float(np.max(np.abs(p.sum(0) - 1.0))))
        width = max(width, hi - lo + 1)
    return m2, norm_err, width


def time_average(m2, Ts):
    """M(T) = int_0^inf (dt/T) e^{-t/T} m(t), trapezoidal rule on TGRID (tail beyond TMAX neglected)."""
    out = np.zeros((len(Ts),) + m2.shape[1:])
    for i, T in enumerate(Ts):
        w = np.exp(-TGRID / T) / T
        out[i] = np.trapezoid(w[:, None] * m2, TGRID, axis=0)
    return out


def validate_time_average(pr):
    """Compare the quadrature of M_2(T) with the exact eigenfunction formula
    M(T) = sum_{a,b} c_a c_b (X^2)_{ab} / (1 + (E_a-E_b)^2 T^2) on a finite chain (dimer, lambda = 0.5)."""
    global S
    S_saved = S
    S = 1
    lam = 0.5
    rng = np.random.default_rng(99)
    N = 2001
    V = potential('dimer', lam, N, rng)[:, 0]
    i0 = 1000
    E, W = eigh_tridiagonal(V, -J * np.ones(N - 1))
    c = W[i0]
    wv = (np.arange(N) - i0) ** 2.0
    A_ = W * c[None, :]
    B = (A_ * wv[:, None]).T @ A_
    # Chebyshev on a window of the same chain (no boundary reached for t <= 400)
    psi = np.zeros((N, 1), complex)
    psi[i0] = 1.0
    r = 2 * J + lam + 1e-9
    lo = hi = i0
    grid = TGRID[TGRID <= 400.0]
    m2 = np.zeros((len(grid), 1))
    for it in range(1, len(grid)):
        lo, hi = cheb_step(V[:, None], psi, lo, hi, grid[it] - grid[it - 1], r)
        y = np.arange(lo, hi + 1) - i0
        m2[it] = (y[:, None] ** 2 * np.abs(psi[lo:hi + 1]) ** 2).sum(0)
    # exact m2(t) at a few times from the eigen-decomposition
    errm = 0.0
    for t in [1.0, 10.0, 100.0, 400.0]:
        k = np.argmin(np.abs(grid - t))
        u = W @ (np.exp(-1j * E * grid[k]) * c)
        errm = max(errm, abs(float(wv @ np.abs(u) ** 2) - m2[k, 0]) / m2[k, 0])
    rel = []
    for T in [1.0, 3.0, 10.0, 30.0]:
        K = 1.0 / (1.0 + (E[:, None] - E[None, :]) ** 2 * T ** 2)
        ex = float(np.sum(B * K))
        w = np.exp(-grid / T) / T
        qd = float(np.trapezoid(w * m2[:, 0], grid))
        rel.append(abs(qd - ex) / ex)
    S = S_saved
    pr(f'# validation (dimer, lambda=0.5, N=2001 chain): Chebyshev vs eigendecomposition, max rel. error of m2(t) '
       f'at t=1,10,100,400: {errm:.1e};  quadrature vs exact M_2(T), T=1,3,10,30: ' +
       ' '.join(f'{x:.1e}' for x in rel))


def compute():
    out = open(os.path.join(HERE, 'dimer_spread.out'), 'w')

    def pr(s):
        print(s)
        out.write(s + '\n')
        out.flush()
    pr(f'# single dressed fermion, sector-averaged <(X-x0)^2>_t; J=1; S={S} sectors per case; x0 odd (left site of a dimer)')
    pr('# one-body dynamics on Z by Chebyshev propagation on an adaptive window (no boundary)')
    validate_time_average(pr)
    res = {}
    for law in LAWS:
        for h in HS:
            t0 = time.time()
            seed = {('dimer', 0.25): 1, ('dimer', 0.75): 2, ('iid', 0.25): 3, ('iid', 0.75): 4}[(law, h)]
            m2, nerr, width = run_case(law, h, seed)
            res[(law, h)] = m2
            MT = time_average(m2, TAVG)
            mean = m2.mean(1)
            se = m2.std(1, ddof=1) / np.sqrt(S)
            pick = [np.argmin(np.abs(TGRID - t)) for t in [10, 100, 1000, 10000]]
            pr(f'{law:5s} h/J={h:.2f} (lambda={2 * h:.1f}): <X^2>_t at t=10,1e2,1e3,1e4: ' +
               '  '.join(f'{mean[k]:.4g}+-{se[k]:.2g}' for k in pick) +
               f'   [max |norm-1| = {nerr:.1e}, max window = {width} sites, {time.time() - t0:.0f}s]')
            # local slope of log <X^2> between decades
            sl = [np.log(mean[pick[i + 1]] / mean[pick[i]]) / np.log(10.0) for i in range(3)]
            pr(f'      slope of log<X^2> vs log t on [10,1e2],[1e2,1e3],[1e3,1e4]: ' + ' '.join(f'{s:.3f}' for s in sl))
            kT = [np.argmin(np.abs(TAVG - T)) for T in [1, 10, 100, 1000]]
            pr(f'      M_2(T) (sector average) at T=1,10,1e2,1e3: ' + '  '.join(f'{MT.mean(1)[k]:.4g}' for k in kT) +
               ';  min over sectors of M_2(T)/T: ' + '  '.join(f'{(MT[k] / TAVG[k]).min():.4g}' for k in kT))
    save_summary(res)
    out.close()


def save_summary(res):
    """Store what the figure needs: mean and standard error of <X^2>_t, and M_2(T) per sector."""
    d = dict(TGRID=TGRID, TAVG=TAVG)
    for (law, h), m2 in res.items():
        d[f'mean_{law}_{h}'] = m2.mean(1)
        d[f'se_{law}_{h}'] = m2.std(1, ddof=1) / np.sqrt(m2.shape[1])
        d[f'MT_{law}_{h}'] = time_average(m2, TAVG)
    np.savez_compressed(DATA, **d)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = np.load(DATA)
    tg, Ts = d['TGRID'], d['TAVG']
    plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                         'ytick.major.width': 0.6, 'legend.frameon': False})
    col = {'iid': '#2a78d6', 'dimer': '#eb6834'}      # validated categorical slots 1-2
    ls = {0.25: '-', 0.75: '--'}
    ink2, grid_c, ref_c = '#52514e', '#d9d8d4', '#8a8984'
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.9), constrained_layout=True)

    ax = axs[0]
    sel = tg > 0
    t = tg[sel]
    ax.loglog(t, 2 * J ** 2 * t ** 2, '-', color=ref_c, lw=0.8, label=r'clean sector $q\equiv1$: $2J^2t^2$')
    for law in ['dimer', 'iid']:
        for h in HS:
            mean = d[f'mean_{law}_{h}'][sel]
            se = d[f'se_{law}_{h}'][sel]
            ax.loglog(t, mean, ls[h], color=col[law], lw=1.6,
                      label=('dimer law' if law == 'dimer' else 'i.i.d. law') + fr', $h={h:g}J$')
            ax.fill_between(t, np.maximum(mean - 2 * se, 1e-6), mean + 2 * se, color=col[law], alpha=0.15, lw=0)
    tt = np.logspace(2, 4, 10)
    ax.loglog(tt, 0.6 * tt ** 1.5, ':', color='#0b0b0b', lw=1.0)
    ax.text(1.5e3, 0.6 * 1.5e3 ** 1.5 * 2.2, r'$\propto t^{3/2}$', fontsize=7.5, color='#0b0b0b')
    ax.set_xlim(1e-1, 1e4)
    ax.set_ylim(1e-2, 1e8)
    ax.set_xlabel(r'time $Jt$')
    ax.set_ylabel(r'$\langle (X-x_0)^2\rangle_t$ (sector average)')
    ax.set_title(r'(a) mean-square displacement of one dressed fermion', fontsize=9)
    ax.grid(True, which='major', color=grid_c, lw=0.5)
    ax.legend(fontsize=7, loc='upper left')

    ax = axs[1]
    for law, h in [('dimer', 0.25), ('iid', 0.25), ('dimer', 0.75)]:
        MT = d[f'MT_{law}_{h}']
        if law == 'dimer' and h == 0.25:
            for s in range(MT.shape[1]):
                ax.loglog(Ts, MT[:, s] / Ts, '-', color=col[law], lw=0.4, alpha=0.35)
        ax.loglog(Ts, MT.mean(1) / Ts, ls[h], color=col[law], lw=1.8,
                  label=('dimer law' if law == 'dimer' else 'i.i.d. law') + fr', $h={h:g}J$' +
                  (' (thin: single sectors)' if (law == 'dimer' and h == 0.25) else ''))
    TT = np.logspace(1, 3, 10)
    ax.loglog(TT, 3.2 * TT ** 0.5, ':', color='#0b0b0b', lw=1.0)
    ax.text(1.3e2, 3.2 * 1.3e2 ** 0.5 * 1.5, r'$\propto T^{1/2}$', fontsize=7.5, color='#0b0b0b')
    ax.set_xlabel(r'averaging time $JT$')
    ax.set_ylabel(r'$M_2(T)/(JT)$')
    ax.set_title(r'(b) time-averaged second moment', fontsize=9)
    ax.set_xlim(1, 1e3)
    ax.grid(True, which='major', color=grid_c, lw=0.5)
    ax.legend(fontsize=7, loc='lower left')
    figdir = os.path.join(os.path.dirname(HERE), 'figs')
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, 'dimer.pdf'))


if __name__ == '__main__':
    if '--plot' not in sys.argv:
        compute()
    plot()
