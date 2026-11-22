import sys
import numpy as np

for tag in sys.argv[1:]:
    d = np.load(f"res_{tag}.npz")
    t = d["tgrid"]
    Ls = sorted({int(k.split('_')[0][1:]) for k in d.files if k.startswith('L')})
    sel = [4, 8, 12, 16, 20, 28, 36, 44, 56]
    print(f"\n===== set {tag} =====  times:", " ".join(f"{t[i]:.0e}" for i in sel))
    for L in Ls:
        SN = d[f"L{L}_SN_t"]; Var = d[f"L{L}_Var_t"]
        Shat = d[f"L{L}_Shat"]; Hbar = d[f"L{L}_Hbar"]; dl = d[f"L{L}_delta"]; Vb = d[f"L{L}_Varbar"]
        n = len(Shat)
        m = SN.mean(0); se = SN.std(0)/np.sqrt(n)
        bound = 0.5*np.log(2*np.pi*np.e*(4*dl**2 + 1/12))
        supS = SN.max(1)
        print(f"L={L:2d} n={n:5d}  <S_N(t)>:", " ".join(f"{m[i]:.3f}" for i in sel),
              f"| sem~{se[sel].max():.3f}")
        print(f"      Shat(inf-time avg)={Shat.mean():.4f}+-{Shat.std()/np.sqrt(n):.4f}  E H(pbar)={Hbar.mean():.4f}+-{Hbar.std()/np.sqrt(n):.4f}"
              f"  frac(Shat<=Hbar+1e-3)={np.mean(Shat <= Hbar + 1e-3):.3f}  <Var_pbar>={Vb.mean():.3f}"
              f"  <Var(t=1e13)>={Var[:, -1].mean():.3f}")
        gap = Hbar - Shat
        print(f"      gap E[H(pbar)-Shat] = {gap.mean():.4f}+-{gap.std()/np.sqrt(n):.4f}")
        print(f"      <delta_circ>={dl.mean():.2f}  <T1 bound 1/2ln(2pie(4d^2+1/12))>={bound.mean():.2f}  min(bound - sup_t S_N)={np.min(bound - supS):.2f}"
              f"  max_t S_N max over samples={supS.max():.3f}  trivial ln(L/2+1)={np.log(L/2+1):.2f}")
        # slope of <S_N> vs ln ln t over windows
        for (a, b) in [(1e1, 1e4), (1e4, 1e8), (1e8, 1e13)]:
            ia, ib = np.argmin(abs(t - a)), np.argmin(abs(t - b))
            x = np.log(np.log(t[ia:ib+1])); y = m[ia:ib+1]
            sl = np.polyfit(x, y, 1)[0]
            print(f"         slope d<S_N>/d lnln t on [{a:.0e},{b:.0e}] = {sl:+.4f}", end="")
        print()
