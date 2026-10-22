"""Haar-random fermionic Gaussian states at small N: direct sampling vs the
finite-N Jacobi formula with the corrected weight x^{-1/2}(1-x)^{N-2l} and with
the weight as printed in the submitted version, x^{-1/2}(1-x)^{N-2l-1/2}."""
import json, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from jacobi_finite_N import avg   # avg(N, l, b): exact average of M_2^GNL, weight x^{-1/2}(1-x)^b

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261022)
m2 = lambda x: -np.log2(1 - x + x * x)

def sample_M2(N, l, n_samples, batch=20000):
    """M_2^GNL = sum_m m_2(lambda_m^2) for Gamma = O Gamma_0 O^T, O Haar on O(2N)."""
    J = np.zeros((2 * N, 2 * N))
    for k in range(N):
        J[2 * k, 2 * k + 1], J[2 * k + 1, 2 * k] = 1.0, -1.0
    out = []
    while sum(len(o) for o in out) < n_samples:
        # first 2l rows of a Haar orthogonal matrix = Haar frame from QR with sign fix
        G = rng.standard_normal((batch, 2 * N, 2 * l))
        Q, R = np.linalg.qr(G)
        Q = Q * np.sign(np.diagonal(R, axis1=1, axis2=2))[:, None, :]
        GA = np.einsum('bia,ij,bjc->bac', Q, J, Q)          # Gamma_A = O_A Gamma_0 O_A^T
        ev = np.linalg.eigvalsh(1j * GA)[:, l:]              # the l non-negative eigenvalues
        out.append(m2(np.clip(ev, 0, 1) ** 2).sum(axis=1))
    v = np.concatenate(out)[:n_samples]
    return v.mean(), v.std(ddof=1) / np.sqrt(n_samples)

JSON = os.path.join(HERE, 'small_N_check.json')
series = {r'$\ell=1$': lambda N: 1, r'$\ell=2$': lambda N: 2, r'$\ell=N/2$': lambda N: N // 2}
Ns = np.arange(2, 17, 2)
data = json.load(open(JSON)) if os.path.exists(JSON) else {}
for name, lf in series.items():
    if name in data:
        continue
    rows = []
    for N in Ns:
        l = int(lf(int(N)))
        if 2 * l > N:
            continue
        mc, err = sample_M2(N, l, 100000)
        rows.append(dict(N=int(N), l=l, mc=float(mc), err=float(err),
                         corrected=float(avg(N, l, N - 2 * l)), printed=float(avg(N, l, N - 2 * l - 0.5))))
    data[name] = rows
json.dump(data, open(JSON, 'w'), indent=1)

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 11})
fig, ax = plt.subplots(figsize=(4.3, 3.9))
colors = ['#1f4e79', '#2e8b57', '#b03a2e']
for (name, rows), c in zip(data.items(), colors):
    rows = [r for r in rows if not (name == r'$\ell=N/2$' and r['N'] < 6)]   # N=2,4 coincide with l=1,2
    N = np.array([r['N'] for r in rows]); l = np.array([r['l'] for r in rows])
    ax.plot(N, np.array([r['corrected'] for r in rows]) / l, '-', c=c, lw=1.4)
    ax.plot(N, np.array([r['printed'] for r in rows]) / l, ':', c=c, lw=1.6)
    ax.errorbar(N, np.array([r['mc'] for r in rows]) / l, yerr=np.array([r['err'] for r in rows]) / l,
                fmt='o', ms=4, mfc='none', mec=c, ecolor=c, capsize=0, label=name)
ax.plot([], [], '-', c='k', lw=1.4, label=r'$w(x)=x^{-1/2}(1-x)^{N-2\ell}$')
ax.plot([], [], ':', c='k', lw=1.6, label=r'$w(x)=x^{-1/2}(1-x)^{N-2\ell-1/2}$')
ax.set_xlabel(r'$N$'); ax.set_ylabel(r'$\overline{M}_2^{\mathrm{GNL}}/\ell$')
ax.set_xticks(Ns)
ax.legend(frameon=False, fontsize=8.5, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2,
          columnspacing=1.0, handlelength=1.8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'small_N_check.pdf')); fig.savefig(os.path.join(HERE, 'small_N_check.png'), dpi=200)
for name, rows in data.items():
    for r in rows:
        print(f"{name:12s} N={r['N']:2d} l={r['l']:2d}  MC={r['mc']:.5f}({r['err']:.5f})  "
              f"corr={r['corrected']:.5f} z={(r['mc']-r['corrected'])/r['err']:+.1f}  "
              f"printed={r['printed']:.5f} z={(r['mc']-r['printed'])/r['err']:+.1f}")
