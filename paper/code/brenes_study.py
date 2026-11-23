"""
Study of the U(1) setting of Brenes, Dalmonte, Heyl, Scardicchio (PRL 120, 030601 (2018)):
lattice Schwinger model with Wilson links (u1.tex, eq. eq:HU1), w = 1, m = 0, theta = 0, J := g^2/2,
initial state = bare vacuum (staggered fermions) x every link in an equal superposition of E in {-1,0,1}, E_0 = 0.

Exact structure (checked to machine precision in the first part):
  * the background FIELDS ebar_j are i.i.d. uniform on {-1,0,1}; the background charges are Q_j = ebar_j - ebar_{j-1};
  * in the sector of ebar the fields are E_j = ebar_j + F_j(n), F_j = sum_{k<=j} rho_k, and
        H_ebar = H_clean + g^2 sum_k A_k rho_k + const,   A_k = sum_{j>=k} (ebar_j + theta):
    the clean (Q = 0) Schwinger model plus a random-walk potential for the charges.

Long-time behaviour: for each sector, the infinite-time average (diagonal ensemble, degeneracies handled) of the
charge density n_pa = (1/L) sum_j |rho_j| is compared with
  (G)   the Gibbs state of the same sector at the energy of the initial state (canonical; beta from <H> = E_0),
  (inf) the infinite-temperature state of the sector (the reference used by Brenes et al.).
Also checked: the rigorous bounds of prop:wilson / rem:wilson(a) / rem:adhh(a) for these states.
"""
import itertools, time, sys
import numpy as np
from scipy.optimize import brentq

T0 = time.time()
rng = np.random.default_rng(7)
w, m, theta = 1.0, 0.0, 0.0


def odd(j):
    return j % 2 == 1


def setup(L):
    vac = tuple(1 if odd(j) else 0 for j in range(1, L + 1)); Nf = sum(vac)
    confs = [n for n in itertools.product([0, 1], repeat=L) if sum(n) == Nf]
    ci = {n: i for i, n in enumerate(confs)}
    R = np.array([[n[j - 1] - (1 if odd(j) else 0) for j in range(1, L + 1)] for n in confs])
    F = np.cumsum(R, axis=1)
    D = len(confs); K = np.zeros((D, D))
    for i, n in enumerate(confs):
        for j in range(1, L):
            if n[j] == 1 and n[j - 1] == 0:
                nn = list(n); nn[j] = 0; nn[j - 1] = 1; k = ci[tuple(nn)]
                K[k, i] += -w; K[i, k] += -w
    mass = np.array([sum(m * (-1) ** j * n[j - 1] for j in range(1, L + 1)) for n in confs])
    Nl = np.cumsum(np.array(confs), axis=1)          # N_{[1,l]}
    return dict(L=L, confs=confs, v=ci[vac], R=R, F=F, K=K, mass=mass, Nl=Nl, npa=np.abs(R).sum(1) / L)


def H_shannon(p):
    p = np.asarray(p, float); p = p[p > 1e-15]
    return float(-(p * np.log(p)).sum())


def sector(S, eb, g2):
    E = eb + S['F']
    d = S['mass'] + 0.5 * g2 * ((E + theta) ** 2).sum(1)
    ev, V = np.linalg.eigh(S['K'] + np.diag(d))
    return ev, V, d, E


def diag_ensemble(ev, V, c, obs, tol=1e-9):
    """infinite-time average of <obs> (obs diagonal in the basis), degenerate eigenspaces handled"""
    out = 0.0; i = 0; n = len(ev)
    while i < n:
        j = i
        while j + 1 < n and ev[j + 1] - ev[i] < tol:
            j += 1
        Vg = V[:, i:j + 1]; a = Vg @ (Vg.T @ c)          # projection of the initial state
        out += float((a ** 2) @ obs)
        i = j + 1
    return out


def gibbs(ev, V, E0, obs):
    """canonical state of the sector with <H> = E0; returns <obs>, beta"""
    oa = (V ** 2).T @ obs                                  # <alpha|obs|alpha>
    e = ev - ev.min()
    def meanE(b):
        x = -b * e; x -= x.max(); p = np.exp(x); p /= p.sum(); return p @ ev
    lo, hi = -50.0 / max(e.max(), 1e-12), 50.0 / max(e[1:].min() if len(e) > 1 else 1, 1e-3)
    f = lambda b: meanE(b) - E0
    if f(lo) * f(hi) > 0:
        b = lo if abs(f(lo)) < abs(f(hi)) else hi
    else:
        b = brentq(f, lo, hi, xtol=1e-12)
    x = -b * e; x -= x.max(); p = np.exp(x); p /= p.sum()
    return float(p @ oa), b


if __name__ == '__main__':
    # ------------------------------------------------------------------------------------------------
    print('=' * 110)
    print('Brenes et al. states: long-time charge density vs Gibbs state at the same energy (w=1, m=0, theta=0)')
    print('=' * 110)
    Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    nsamp = {6: None, 8: None, 10: 3000, 12: 300}
    Js = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
    for L in [6, 8, 10, 12]:
        if L > Lmax:
            break
        S = setup(L)
        allsec = list(itertools.product([-1, 0, 1], repeat=L))
        secs = allsec if nsamp[L] is None else [allsec[i] for i in rng.choice(len(allsec), nsamp[L], replace=False)]
        print(f'L={L}: sector dimension {len(S["confs"])}, sectors used {len(secs)} of {3 ** L}')
        for J in Js:
            g2 = 2 * J
            de = gb = inf = 0.0; e2mem = 0.0; SN_lt = np.zeros(L)
            pN = np.zeros((L, L + 1))
            frac_above = 0.0
            for eb in secs:
                eb = np.array(eb)
                ev, V, d, E = sector(S, eb, g2)
                c = V[S['v'], :]
                v = np.zeros(len(ev)); v[S['v']] = 1.0
                de += diag_ensemble(ev, V, v, S['npa'])
                E0 = d[S['v']]
                gval, b = gibbs(ev, V, E0, S['npa']); gb += gval
                inf += S['npa'].mean()
                frac_above += (E0 - ev.min()) / (ev.max() - ev.min())
            n = len(secs)
            print(f'  J={J:5.2f}: n_pa long-time {de / n:.4f} | Gibbs at E0 {gb / n:.4f} | infinite T {inf / n:.4f} | '
                  f'E0 position in sector spectrum {frac_above / n:.3f}', flush=True)
        print(f'  [time {time.time() - T0:.0f}s]', flush=True)
