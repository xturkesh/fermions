"""Brenes et al. states, part 3: (D) the large-J plateau of the long-time charge density versus L, with the fraction of
sectors where the vacuum mixes appreciably (long-time n_pa > 0.05); (A) L = 12, sampled sectors, J = 0.5 ... 5."""
import itertools, time, sys
import numpy as np
sys.path.insert(0, '.')
from brenes_study import setup, sector, diag_ensemble, gibbs
T0 = time.time(); rng = np.random.default_rng(5)
part = sys.argv[1] if len(sys.argv) > 1 else 'DA'
if 'D' in part:
    print('(D) J = 80: long-time n_pa versus L', flush=True)
    for L, ns in [(6, None), (8, None), (10, 2000), (12, 200)]:
        S = setup(L)
        allsec = list(itertools.product([-1, 0, 1], repeat=L)) if ns is None else None
        secs = [np.array(e) for e in allsec] if ns is None else [rng.integers(-1, 2, size=L) for _ in range(ns)]
        vals = []
        for eb in secs:
            ev, V, d, E = sector(S, eb, 160.0)
            v = np.zeros(len(ev)); v[S['v']] = 1.0
            vals.append(diag_ensemble(ev, V, v, S['npa']))
        vals = np.array(vals)
        print(f'  L={L:2d} ({len(secs)} sectors): mean {vals.mean():.5f}, median {np.median(vals):.2e}, '
              f'fraction > 0.05: {np.mean(vals > 0.05):.4f}, fraction > 0.2: {np.mean(vals > 0.2):.4f}  [{time.time() - T0:.0f}s]', flush=True)
if 'A' in part:
    print('(A) L=12, 300 random sectors: long-time n_pa vs Gibbs at E0', flush=True)
    S = setup(12)
    secs = [rng.integers(-1, 2, size=12) for _ in range(300)]
    for J in [0.5, 1.0, 2.0, 5.0]:
        de = gb = 0.0
        for eb in secs:
            ev, V, d, E = sector(S, eb, 2 * J)
            v = np.zeros(len(ev)); v[S['v']] = 1.0
            de += diag_ensemble(ev, V, v, S['npa']); gb += gibbs(ev, V, d[S['v']], S['npa'])[0]
        print(f'  J={J:5.2f}: n_pa long-time {de / 300:.4f} | Gibbs at E0 {gb / 300:.4f} | gap {(gb - de) / 300:.4f}  [{time.time() - T0:.0f}s]', flush=True)
print(f'[time {time.time() - T0:.0f}s]')
