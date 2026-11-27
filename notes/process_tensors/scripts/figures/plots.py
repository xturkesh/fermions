import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from m1core import I2l
d, q = 2, 2
def I2inf(L, p):
    D = float(q)**L; u = (1-p)**2; GI, GF = D, (u*q*q + 1 - u)**L
    a = (D*D-1)/(d*d*D*D-1); b = D*(d*d-1)/(d*d*D*D-1)
    cI = (D*GI - GF)/(D*(D*D-1)); cF = (D*GF - GI)/(D*(D*D-1))
    T = np.array([[a, b], [cI*a, cI*b + cF]])
    w, R = np.linalg.eig(T); i = np.argmax(w.real); r = R[:, i].real
    wl, Lv = np.linalg.eig(T.T); j = np.argmax(wl.real); l = Lv[:, j].real
    l = l/(l @ r); ys = (cI + cF*b)/(1 - cF*a); wv = np.array([1.0, ys])
    return -np.log2((r[0])*(l @ wv))

OUT = '../../src/figs/'
TXT, TXT2, GRID, SPINE = '#191919', '#666666', '#e8e8e8', '#bdbdbd'
CAT = ['#257de0', '#f47230', '#00b894', '#6a4ce2']                 # validated categorical subset
ORD = ['#84b9f3', '#338ced', '#2370cd', '#1757ac', '#0e3d7e']      # validated ordinal blue ramp
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 8.5,
    'axes.edgecolor': SPINE, 'axes.labelcolor': TXT, 'xtick.color': TXT2, 'ytick.color': TXT2,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
    'lines.linewidth': 1.4, 'lines.markersize': 4.5, 'legend.frameon': False, 'legend.fontsize': 7.5,
    'axes.spines.top': False, 'axes.spines.right': False, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
def tag(ax, s): ax.text(-0.13, 1.04, s, transform=ax.transAxes, fontsize=9.5, fontweight='bold', color=TXT, va='bottom')
def label_end(ax, x, y, s, c, dx=0.0, dy=0.0, ha='left'):
    ax.text(x + dx, y + dy, s, color=TXT, fontsize=7.5, ha=ha, va='center',
            bbox=dict(boxstyle='round,pad=0.12', fc='white', ec=c, lw=0.8))
m1 = np.load('data_m1_ipr.npz'); m2 = np.load('data_m2.npz'); mc = np.load('data_ipr_mc.npz')
d, q = 2, 2

# ---------------- Figura modello 1 ----------------
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.45))
ax = axs[0]; Ls = np.arange(1, 23)
ax.plot(Ls, 2*Ls, ls='--', color=TXT2, lw=1.1)
label_end(ax, 1.5, 21.0, r'senza rumore: $2\log_2 D$', TXT2)
for i, p in enumerate([0.05, 0.1, 0.3]):
    y = np.array([I2inf(L, p) for L in Ls]); c = ORD[1 + i]
    ax.plot(Ls, y, '-o', color=c, ms=3.2, mec='white', mew=0.5)
    gF = (1-p)**2*q*q + 1 - (1-p)**2; Lstar = np.log(d*d)/np.log(q*q/gF)
    ax.axvline(Lstar, color=c, lw=0.8, ls=':')
    k = np.argmax(y); label_end(ax, Ls[k], y[k] + 2.4, r'$p=%.2g$' % p, c, ha='center')
ax.set_xlabel('qudit dell’ambiente $L$  ($D=2^L$)'); ax.set_ylabel(r'$I_2^{\infty}(\mathrm{passato:futuro})$ [bit]')
ax.set_xlim(0.5, 22.5); ax.set_ylim(-0.5, 31); tag(ax, '(a)')
ax = axs[1]
for i, p in enumerate([0.1, 0.03, 0.01, 0.005]):
    x, y = m1['P2_x_p%.3f' % p], m1['P2_y_p%.3f' % p]
    ax.plot(x, y, color=ORD[1 + i], label=r'$p=%g$' % p)
xs = np.log(d*d)/(2*(1 - q**-2.0)); ax.axvline(xs, color=TXT2, ls='--', lw=1.0)
ax.text(xs + 0.03, 0.55, r'$x^*=%.3f$' % xs, color=TXT, fontsize=7.5)
ax.set_xlabel(r'$x=pL$'); ax.set_ylabel(r'$I_2^{\infty}/(2\log_2 D)$'); ax.set_xlim(0, 1.6); ax.set_ylim(-0.03, 1.05)
ax.legend(loc='lower left', title=r'verso il limite di scala', title_fontsize=7.5); tag(ax, '(b)')
fig.tight_layout(); fig.savefig(OUT + 'fig_model1.pdf'); plt.close(fig)

# ---------------- Figura modello 2: dinamica ----------------
lam = float(m2['lam'])
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.45))
ax = axs[0]; t = m2['A_t'][:21]
for i, p in enumerate([0.0, 0.2, 0.7]):
    yv = np.abs(m2['A_p%.1f' % p][:21])/lam**(2*t)
    env = np.array([yv[max(0, j-1):j+2].max() for j in range(len(yv))])[:-1]; yv = yv[:-1]; tt = t[:-1]     # envelope (complex eigenvalues -> oscillations)
    ax.plot(tt, yv, 'o', color=CAT[i], ms=2.4, alpha=0.35, mew=0)
    ax.plot(tt, env, '-', color=CAT[i], label=r'$p=%g$' % p)
ax.text(12.5, 21.5, r'$\propto t$: blocco di Jordan', color=TXT, fontsize=7.5)
ax.text(13.5, 6.0, r'limitato: puro $\lambda^{2t}$', color=TXT, fontsize=7.5)
ax.set_xlabel(r'tempo $t$'); ax.set_ylabel(r'$|C_2(t)|\,/\,\lambda^{2t}$  (inviluppo)'); ax.set_xlim(0, 20.5); ax.set_ylim(0, 28)
ax.legend(loc='upper left'); tag(ax, '(a)')
ax = axs[1]; pg = m2['B_p']; pf = np.linspace(0, 1, 400); u = (1 - pf)**2
ax.plot(pf, np.maximum(u, lam**2), color=CAT[0], lw=1.2)
ax.plot(pf, np.maximum(u, lam), color=CAT[1], lw=1.2)
ax.plot(pg, m2['B_traceless'], 'o', color=CAT[0], mec='white', mew=0.5, label=r'$A$ a traccia nulla: $\max(u,\lambda^2)$')
ax.plot(pg, m2['B_generic'], 's', color=CAT[1], mec='white', mew=0.5, ms=3.8, label=r'$A$ generico: $\max(u,\lambda)$')
for pc, c in [(1 - lam, CAT[0]), (1 - np.sqrt(lam), CAT[1])]:
    ax.axvline(pc, color=c, ls=':', lw=0.9)
ax.text(1 - lam + 0.015, 0.56, r'$p_c=1-\lambda$', fontsize=7.2, color=TXT)
ax.text(1 - np.sqrt(lam) + 0.015, 0.66, r'$p_c=1-\sqrt{\lambda}$', fontsize=7.2, color=TXT)
ax.set_xlabel(r'rumore sul bagno $p$  ($u=(1-p)^2$)'); ax.set_ylabel('tasso di rilassamento'); ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.legend(loc='upper right', fontsize=7); tag(ax, '(b)')
fig.tight_layout(); fig.savefig(OUT + 'fig_model2_rates.pdf'); plt.close(fig)

# ---------------- Figura modello 2: free probability ----------------
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.45))
ax = axs[0]; t = m2['C_t']
for i, p in enumerate([0.0, 0.01, 0.03, 0.1]):
    ax.semilogx(t, m2['C_p%.2f' % p], color=ORD[1 + i] if p > 0 else CAT[1], label=r'$p=%g$' % p)
ax.axhline(0.5625, color=CAT[1], ls='--', lw=0.9); ax.axhline(0.3125, color=ORD[4], ls='--', lw=0.9)
import matplotlib.transforms as mtr
tr = mtr.blended_transform_factory(ax.transAxes, ax.transData)
ax.text(1.01, 0.5625, 'libero', transform=tr, fontsize=7.5, color=TXT, va='center')
ax.text(1.01, 0.3125, 'monotono', transform=tr, fontsize=7.5, color=TXT, va='center')
ax.set_xlabel(r'tempo $t$'); ax.set_ylabel(r'$C_2(t)$,  $A=Z+\frac{1}{2}$, $B=X+\frac{1}{2}$'); ax.set_ylim(0.0, 0.75)
ax.legend(loc='lower right', ncol=2, fontsize=7); tag(ax, '(a)')
ax = axs[1]
for i, (key, tk, n, lab) in enumerate([('D2', 'D2_t', 2, r'$\kappa_2$ ($k=2$)'), ('D3', 'D3_t', 3, r'$\kappa_3$ ($k=3$, qutrit)'), ('D4', 'D4_t', 4, r'$\kappa_4,\kappa_2^2$ ($k=4$)')]):
    tt, yy = m2[tk], np.abs(m2[key]); c = CAT[i]
    ax.semilogy(tt, yy, 'o', color=c, mec='white', mew=0.5, label=lab)
    ax.semilogy(tt, yy[-1]*(0.9**n)**(tt - tt[-1]), '-', color=c, lw=1.0)
ax.text(0.97, 0.95, r'linee: $(1-p)^{n t}$,  $p=0.1$', transform=ax.transAxes, ha='right', va='top', fontsize=7.5, color=TXT2)
ax.set_xlabel(r'tempo $t$'); ax.set_ylabel(r'termine con $\kappa_n(A)$'); ax.legend(loc='lower left', fontsize=7); tag(ax, '(b)')
fig.tight_layout(); fig.savefig(OUT + 'fig_model2_free.pdf'); plt.close(fig)

# ---------------- Figura IPR temporale ----------------
fig, axs = plt.subplots(2, 2, figsize=(6.3, 4.6))
ax = axs[0, 0]; n = m1['P6_n']
ax.loglog(n, m1['P6_pure_D1'], color=CAT[1], label=r'$D=1$')
for i, D in enumerate([2, 8, 32, 128]):
    ax.loglog(n, m1['P6_pure_D%d' % D], color=ORD[1 + i], label=r'$D=%d$' % D)
ax.axhline(1, color=TXT2, lw=0.9, ls=':'); ax.axhline(2, color=TXT2, lw=0.9, ls='--')
import matplotlib.transforms as mtr
tr = mtr.blended_transform_factory(ax.transAxes, ax.transData)
ax.text(1.01, 1.0, 'uniforme', transform=tr, fontsize=7, color=TXT, va='center'); ax.text(1.01, 2.0, 'PT', transform=tr, fontsize=7, color=TXT, va='center')
ax.set_ylim(0.9, 1e3); ax.set_xlabel(r'numero di misure $n$'); ax.set_ylabel(r'$N\,\mathcal{I}_2$  (bagno puro)')
ax.legend(loc='upper left', ncol=1, fontsize=7); tag(ax, '(a)')
ax = axs[0, 1]; bins = np.linspace(-8, 2.5, 70)
for i, D in enumerate([1, 4, 16, 64]):
    ax.hist(mc['H_D%d' % D], bins=bins, density=True, histtype='step', lw=1.3, color=[CAT[1], ORD[1], ORD[2], ORD[4]][i], label=r'$D=%d$' % D)
yy = np.linspace(-8, 2.5, 300); ax.plot(yy, np.exp(yy)*np.exp(-np.exp(yy)), '--', color=TXT2, lw=1.0, label='Porter–Thomas')
ax.set_xlabel(r'$\ln(Np)$,  $n=12$ misure'); ax.set_ylabel('densità'); ax.set_xlim(-8, 2.5); ax.legend(loc='upper left', fontsize=7); tag(ax, '(b)')
ax = axs[1, 0]; V = mc['V']
for i, D in enumerate([4, 16, 64]):
    sel = V[:, 0] == D; xv = (1 - 1/d)*V[sel, 1]/D
    ax.plot(xv, V[sel, 2], 'o', color=ORD[2 + i], mec='white', mew=0.5, label=r'$D=%d$' % D)
xx = np.linspace(0, 1.6, 10); ax.plot(xx, xx, '--', color=TXT2, lw=1.0); ax.text(1.0, 0.82, r'$\sigma^2=(1-\frac{1}{d})\,\frac{n}{D}$', fontsize=7.5, color=TXT)
ax.set_xlabel(r'$(1-1/d)\,n/D$'); ax.set_ylabel(r'$\sigma^2=\mathrm{Var}\,\ln(Np)$'); ax.legend(loc='upper left', fontsize=7); tag(ax, '(c)')
ax = axs[1, 1]; n1 = m1['P8_n1']
for i, D in enumerate([2, 8, 32, 128]):
    y = m1['P8_D%d' % D]; ok = np.isfinite(y); c = ORD[1 + i]
    ax.semilogx(n1[ok], y[ok], color=c)
    j = np.argmin(np.abs(y[ok] - 0.45)); label_end(ax, n1[ok][j]*0.8, 0.45 + 0.07*(i % 2 == 0) - 0.07*(i % 2 == 1), r'$D=%d$' % D, c, ha='right')
ax.axhline(1, color=TXT2, ls='--', lw=0.9); ax.text(600, 1.035, '1 bit', fontsize=7, color=TXT, ha='right')
ax.set_xlabel(r'$n_1=n_2$'); ax.set_ylabel(r'memoria nei record $\Delta I_{\rm rec}$ [bit]'); ax.set_ylim(0, 1.15); ax.set_xlim(0.8, 800); tag(ax, '(d)')
fig.tight_layout(); fig.savefig(OUT + 'fig_ipr.pdf'); plt.close(fig)
print("figures written")
