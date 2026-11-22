"""
Figures for Section sec:numerics (reads data/numerics_*.npz written by numerics_imbalance.py,
numerics_entropies_small.py and numerics_entropies_large.py):
    figs/imbalance.pdf   (fig:imbalance)
    figs/entropies.pdf   (fig:entropies)
Colors: ordinal one-hue ramps for ordered parameters (lambda, N) and the Okabe-Ito categorical colors for the
four entropies; all palettes were checked with a CVD/contrast validator.  Usage: python3 numerics_plots.py
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from numerics_imbalance import imbalance_bounds, uniform_bound   # noqa: E402

DATA = os.path.join(HERE, "data")
FIGS = os.path.normpath(os.path.join(HERE, "..", "figs"))
LN2 = np.log(2.0)

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
    "axes.unicode_minus": False, "axes.formatter.use_mathtext": True,
    "font.size": 8.5, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 7.3, "legend.frameon": False, "legend.handlelength": 1.8,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "lines.linewidth": 1.25, "pdf.fonttype": 42,
})

INK2 = "#52514e"      # secondary ink
MUTED = "#9a9893"     # muted gray
AXIS0 = "#c3c2b7"     # zero line
LAM_RAMP = {0.5: "#86b6ef", 1.0: "#3987e5", 2.0: "#1c5cab", 4.0: "#0d366b"}    # ordinal (blue)
N_RAMP_IB = {100: "#f2a05f", 200: "#d55e00", 400: "#7f3300"}                   # ordinal (vermillion)
C_SG, C_SN, C_SF, C_SM = "#0072B2", "#D55E00", "#009E73", "#CC79A7"            # Okabe-Ito
SG_RAMP = {64: "#7fb8de", 128: "#0072B2", 256: "#003a5c"}
SN_RAMP = {64: "#f2a05f", 128: "#d55e00", 256: "#7f3300"}


def logx(ax, lo, hi, every=1):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(LogLocator(base=10.0 ** every, numticks=30))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=(1.0,) if every > 1 else np.arange(2, 10) * 0.1,
                                          numticks=100))
    ax.xaxis.set_minor_formatter(NullFormatter())


def panel_label(ax, s, x=0.025, y=0.965, **kw):
    ax.text(x, y, s, transform=ax.transAxes, ha="left", va="top", **kw)


# ============================================================================== imbalance
def fig_imbalance():
    tr = np.load(os.path.join(DATA, "numerics_imbalance_traces.npz"))
    ib = {N: np.load(os.path.join(DATA, f"numerics_imbalance_ibar_N{N}.npz")) for N in (100, 200, 400)}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.7), gridspec_kw=dict(width_ratios=[1.12, 1.0]))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.16, top=0.97, wspace=0.27)

    # (a) I(t), N = 200
    t = tr["times"]
    ax1.axhline(0.0, color=AXIS0, lw=0.6, zorder=0)
    m = t <= 40.0
    ax1.plot(t[m], tr["clean"][m], color=MUTED, lw=0.8, zorder=1)
    ax1.text(0.95, -0.36, r"$q\equiv1$", color=INK2, fontsize=7.5, ha="left", va="center")
    for i, lam in enumerate(tr["lams"]):
        mu, se = tr["mean"][i], tr["se"][i]
        col = LAM_RAMP[float(lam)]
        ax1.fill_between(t, mu - 2 * se, mu + 2 * se, color=col, alpha=0.35, lw=0, zorder=2)
        ax1.plot(t, mu, color=col, zorder=3, label=rf"$\lambda={lam:g}$")
        ax1.plot([1.5e5, 1e6], [tr["ibar"][i]] * 2, color="k", lw=0.7, ls=(0, (1.5, 1.5)), zorder=4)
    for lam, y in ((4.0, 0.285), (2.0, 0.205), (1.0, 0.115), (0.5, 0.045)):
        ax1.text(9e3, y, rf"$\lambda={lam:g}$", color=INK2, fontsize=7.5, ha="right", va="center")
    logx(ax1, 0.1, 1e6)
    ax1.set_ylim(-0.45, 1.03)
    ax1.yaxis.set_major_locator(MultipleLocator(0.25))
    ax1.set_xlabel(r"time $Jt$")
    ax1.set_ylabel(r"imbalance $I(t)$")
    panel_label(ax1, r"(a) $N=200$", x=0.6, y=0.95)

    # (b) Ibar(lambda) and the lower bounds of thm:imbalance
    marks = {100: "o", 200: "s", 400: "D"}
    for N in (100, 200, 400):
        d = ib[N]
        ax2.errorbar(d["lams"], d["ibar"], yerr=2 * d["ibar_se"], fmt=marks[N], ms=3.2, mew=0.8,
                     color=N_RAMP_IB[N], mfc="white" if N != 400 else N_RAMP_IB[N], lw=0.8, capsize=0,
                     label=rf"$\bar I$, $N={N}$", zorder=3)
    lams = ib[400]["lams"]
    hs = lams / 2.0
    k400 = ib[400]["kappa"] + 2 * ib[400]["kappa_se"]
    best3 = np.array([max(imbalance_bounds(400, h, k)[:2]) for h, k in zip(hs, k400)])
    khat = np.max(np.vstack([ib[N]["kappa"] + 2 * ib[N]["kappa_se"] for N in (100, 200, 400)]), axis=0)
    unif = np.array([uniform_bound(h, k)[0] for h, k in zip(hs, khat)])
    ax2.plot(lams, best3, color=INK2, lw=1.0, marker="|", ms=4, zorder=2,
             label=r"bound (iii), $N=400$")
    ax2.plot(lams, unif, color=INK2, lw=1.0, ls=(0, (4, 2)), zorder=2,
             label=r"bound (iv), $\kappa=\hat\kappa$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(0.2, 10)
    ax2.set_ylim(1e-9, 1.0)
    ax2.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    ax2.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8"])
    ax2.xaxis.set_minor_formatter(NullFormatter())
    ax2.yaxis.set_major_locator(LogLocator(base=100.0, numticks=20))
    ax2.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(1.0,), numticks=40))
    ax2.yaxis.set_minor_formatter(NullFormatter())
    ax2.set_xlabel(r"disorder strength $\lambda=2h/J$")
    ax2.set_ylabel(r"infinite-time imbalance $\bar I$")
    ax2.legend(loc="lower right", ncol=1, borderaxespad=0.3)
    panel_label(ax2, "(b)")
    out = os.path.join(FIGS, "imbalance.pdf")
    fig.savefig(out)
    plt.close(fig)
    return out, dict(best3=best3, unif=unif, khat=khat)


# ============================================================================== entropies
def fig_entropies():
    sm = np.load(os.path.join(DATA, "numerics_entropies_small.npz"))
    lg = {lam: np.load(os.path.join(DATA, f"numerics_entropies_large_lam{lam:g}.npz")) for lam in (1.0, 2.0)}
    fig, axs = plt.subplots(2, 2, figsize=(6.6, 4.9))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.085, top=0.985, wspace=0.22, hspace=0.3)
    (axa, axb), (axc, axd) = axs

    # (a) N = 10 exact
    t = sm["times"][1:]
    ell = int(sm["ell"])
    axa.axhline(2 * ell * LN2, color=AXIS0, lw=0.7, zorder=0)
    axa.axhline(ell * LN2, color=AXIS0, lw=0.7, zorder=0)
    axa.text(1.3e-2, 2 * ell * LN2 - 0.08, r"$2\ell\ln2$", color=INK2, fontsize=7.3, va="top")
    axa.text(1.3e-2, ell * LN2 + 0.08, r"$\ell\ln2$", color=INK2, fontsize=7.3, va="bottom")
    series = (("Sfull", C_SF, r"$S_{\mathrm{full}}$", "o"), ("Smix", C_SM, r"$S_{\mathrm{mix}}$", "s"),
              ("SG", C_SG, r"$S_{\mathrm{G}}$", "^"), ("SN", C_SN, r"$S_N$", "v"))
    for key, col, lab, mk in series:
        axa.plot(t, sm[key][1:], color=col, label=lab, marker=mk, ms=3.0, markevery=(4, 8), mew=0.6,
                 mfc="white")
    axa.text(2e5, sm["Sfull"][-1] + 0.25, r"$S_{\mathrm{full}}$", color=INK2, fontsize=8, ha="center")
    axa.text(2e5, sm["Smix"][-1] + 0.22, r"$S_{\mathrm{mix}}$", color=INK2, fontsize=8, ha="center")
    axa.text(2e5, sm["SG"][-1] + 0.2, r"$S_{\mathrm{G}}$", color=INK2, fontsize=8, ha="center")
    axa.text(2e5, sm["SN"][-1] - 0.22, r"$S_N$", color=INK2, fontsize=8, ha="center", va="top")
    logx(axa, 1e-2, 1e6, every=2)
    axa.set_ylim(0, 7.4)
    axa.set_xlabel(r"time $Jt$")
    axa.set_ylabel("entropy (nats)")
    axa.legend(loc="upper left", bbox_to_anchor=(0.0, 0.9), ncol=1)
    panel_label(axa, r"(a) $N=10$, $\ell=5$, $h=J$: exact")

    # (b) S_N - S_G at early times
    axb.axhline(0.0, color=AXIS0, lw=0.7, zorder=0)
    tl = lg[2.0]["times"]
    for lam, ls, lab in ((1.0, (0, (4, 2)), r"$N=256$, $\lambda=1$"), (2.0, "-", r"$N=256$, $\lambda=2$")):
        d = lg[lam]
        y = d["SN_N256"] - d["SG_N256"]
        e = d["DNG_se_N256"]
        col = LAM_RAMP[lam]
        axb.fill_between(tl[1:], (y - 2 * e)[1:], (y + 2 * e)[1:], color=col, alpha=0.3, lw=0)
        axb.plot(tl[1:], y[1:], color=col, ls=ls, label=lab)
    axb.plot(t, (sm["SN"] - sm["SG"])[1:], color="k", lw=1.1, label=r"$N=10$, $\lambda=2$ (exact)")
    logx(axb, 1e-2, 1e2)
    axb.set_ylim(-0.32, 0.17)
    axb.set_xlabel(r"time $Jt$")
    axb.set_ylabel(r"$S_N-S_{\mathrm{G}}$")
    axb.text(0.03, 0.62, r"$S_{\mathrm{G}}<S_N$", transform=axb.transAxes, color=INK2, fontsize=7.5)
    axb.text(0.03, 0.47, r"$S_{\mathrm{G}}>S_N$", transform=axb.transAxes, color=INK2, fontsize=7.5)
    axb.legend(loc="lower left", borderaxespad=0.3)
    panel_label(axb, "(b)")

    # (c), (d) large N
    for ax, lam, lab in ((axc, 1.0, "(c)"), (axd, 2.0, "(d)")):
        d = lg[lam]
        tt = d["times"][1:]
        for N in (64, 128, 256):
            for key, ramp, sym in (("SG", SG_RAMP, r"S_{\mathrm{G}}"), ("SN", SN_RAMP, r"S_N")):
                y = d[f"{key}_N{N}"][1:]
                e = d[f"{key}_se_N{N}"][1:]
                ax.fill_between(tt, y - 2 * e, y + 2 * e, color=ramp[N], alpha=0.3, lw=0)
                ax.plot(tt, y, color=ramp[N], lw=1.0 if N != 256 else 1.25,
                        label=rf"${sym}$, $N={N}$")
        logx(ax, 1e-2, 1e8, every=2)
        ax.set_xlabel(r"time $Jt$")
        ax.set_ylabel("entropy (nats)")
        panel_label(ax, rf"{lab} $\lambda={lam:g}$, sector average")
    axc.set_ylim(0, 3.7)
    axd.set_ylim(0, 1.6)
    handles, labels = axc.get_legend_handles_labels()
    order = [0, 2, 4, 1, 3, 5]
    axc.legend([handles[i] for i in order], [labels[i] for i in order], ncol=2, loc="center left",
               bbox_to_anchor=(0.0, 0.62), columnspacing=1.0)
    out = os.path.join(FIGS, "entropies.pdf")
    fig.savefig(out)
    plt.close(fig)
    return out


def theta_sigma(F):
    """constants of eq:main-constants for a non-increasing F on {0,...,len(F)-1} (truncated)."""
    d = np.arange(len(F))
    S1 = float((d[1:] * F[1:]).sum())
    S1p = float((d[1:] * np.sqrt(F[1:])).sum())
    Th = float(2 * (d[1:] * (2 * d[1:] + 1) * F[1:]).sum())
    M2 = float(((d + 1.0) ** 2 * F).sum())
    return S1, S1p, Th, M2


def summary():
    """numbers quoted in sections/numerics.tex (printed; see numerics_plots.out)."""
    ib = {N: np.load(os.path.join(DATA, f"numerics_imbalance_ibar_N{N}.npz")) for N in (100, 200, 400)}
    lams = ib[400]["lams"]
    print("=== imbalance: Ibar(N) +- SE, kappa_N +- SE (max over all intervals), lower bounds of thm:imbalance ===")
    print("    bounds use kappa_N + 2 SE; khat = max_N (kappa_N + 2 SE); kthm = 2 Sigma_1(F_emp) with F_emp from N=400")
    khat = np.max(np.vstack([ib[N]["kappa"] + 2 * ib[N]["kappa_se"] for N in (100, 200, 400)]), axis=0)
    for i, lam in enumerate(lams):
        h = lam / 2.0
        row = "  ".join(f"{ib[N]['ibar'][i]:.5f}({ib[N]['ibar_se'][i]:.5f})" for N in (100, 200, 400))
        kap = "  ".join(f"{ib[N]['kappa'][i]:.2f}({ib[N]['kappa_se'][i]:.2f})" for N in (100, 200, 400))
        k4 = ib[400]["kappa"][i] + 2 * ib[400]["kappa_se"][i]
        glob, blk, Lb = imbalance_bounds(400, h, k4)
        u, Ls = uniform_bound(h, khat[i])
        S1 = float(ib[400]["S1emp"][i])
        tail = float(ib[400]["Femp"][i][-1])
        ut, Lt = uniform_bound(h, 2 * S1)
        print(f"  lambda={lam:4.2f}: Ibar N=100,200,400: {row} | kappa_N: {kap}")
        print(f"      N=400: global {glob:.2e}, block {blk:.2e} (L={Lb}), best(iii) {max(glob, blk):.2e}, "
              f"ratio Ibar/best = {ib[400]['ibar'][i] / max(glob, blk):.0f} | (iv) khat={khat[i]:.2f}: L*={Ls}, "
              f"{u:.2e}, ratio {ib[400]['ibar'][i] / u:.1e} | Sigma_1(F_emp)={S1:.3g} (F_emp(399)={tail:.1e}), "
              f"(iv) with kthm={2*S1:.3g}: {ut:.2e}")
    print("=== constants of thm:arealaw(b) and thm:numberentropy with F_emp (N=400, truncated at d<=399) ===")
    for lam in (1.0, 2.0):
        i = int(np.argmin(np.abs(lams - lam)))
        F = ib[400]["Femp"][i]
        S1, S1p, Th, M2 = theta_sigma(F)
        print(f"  lambda={lam}: Sigma_1={S1:.3g}, Sigma'_1={S1p:.3g}, Theta={Th:.3g}, M_2={M2:.3g}; "
              f"area-law bound 2ln2*Theta = {2*LN2*Th:.3g} (<= 8ln2*M_2 = {8*LN2*M2:.3g}); "
              f"number-entropy bound 1/2 ln(2 pi e (Theta + Sigma'_1^2 + 1/12)) = "
              f"{0.5*np.log(2*np.pi*np.e*(Th + S1p**2 + 1/12)):.2f}; per-sector bound {0.5*np.log(2*np.pi*np.e*(Th + 1/12)):.2f}")
    for lam in (1.0, 2.0):
        f = os.path.join(DATA, f"numerics_entropies_large_lam{lam:g}.npz")
        if not os.path.exists(f):
            continue
        d = np.load(f)
        t = d["times"]
        late = t >= 1e3
        print(f"=== large N, lambda={lam}: plateaus (time average over grid points t>=1e3; error over sectors) ===")
        for N in (64, 128, 256):
            vm = d[f"var_mix_N{N}"]
            print(f"  N={N}: S_G = {float(d[f'SG_plateau_N{N}']):.4f}({float(d[f'SG_plateau_se_N{N}']):.4f}), "
                  f"S_N = {float(d[f'SN_plateau_N{N}']):.4f}({float(d[f'SN_plateau_se_N{N}']):.4f}), "
                  f"E_qH(p^q) = {d[f'EHq_N{N}'][late].mean():.3f}, Var_mix = {vm[late].mean():.3f} (max_t {vm.max():.3f}), "
                  f"Massey with measured Var (non-rigorous) max_t = {0.5*np.log(2*np.pi*np.e*(vm.max()+1/12)):.3f}, "
                  f"E_q max_grid S(rho^q_A) = {float(d[f'SGsup_N{N}']):.3f}({float(d[f'SGsup_se_N{N}']):.3f}), "
                  f"max_t S_N = {d[f'SN_N{N}'].max():.3f}, max_t S_G = {d[f'SG_N{N}'].max():.3f}")
            dd = d[f"SN_N{N}"] - d[f"SG_N{N}"]
            e = d[f"DNG_se_N{N}"]
            sig = (dd > 2 * e) & (t > 0)
            print(f"      S_N - S_G: max = {dd.max():.4f}({e[np.argmax(dd)]:.4f}) at t = {t[np.argmax(dd)]:.3g}; "
                  f"> 2 SE at grid times: " + (f"[{t[sig].min():.3g}, {t[sig].max():.3g}]" if sig.any() else "none")
                  + "; values at t=1, 2.15, 4.64, 10: " + ", ".join(
                      f"{dd[np.argmin(abs(t - tt))]:+.4f}({e[np.argmin(abs(t - tt))]:.4f})" for tt in (1, 2.154, 4.642, 10)))


if __name__ == "__main__":
    if sys.argv[1:] == ["summary"]:
        summary()
        sys.exit(0)
    which = sys.argv[1:] or ["imbalance", "entropies"]
    if "imbalance" in which:
        out, b = fig_imbalance()
        print("wrote", out)
        for lam, x, y, k in zip(np.load(os.path.join(DATA, "numerics_imbalance_ibar_N400.npz"))["lams"],
                                b["best3"], b["unif"], b["khat"]):
            print(f"  lambda={lam:5.2f}: best bound (iii) at N=400 = {x:.2e}; bound (iv) with kappa_hat={k:.2f}: {y:.2e}")
    if "entropies" in which:
        print("wrote", fig_entropies())
