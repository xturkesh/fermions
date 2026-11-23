"""Brenes et al. states, part 2: (A) L=12 (sampled sectors); (B) large-J scaling of the long-time charge density;
(C) the electric-energy band of the vacuum: band size, in-band infinite-temperature density, and the long-time density
under the second-order (Schrieffer-Wolff) effective in-band Hamiltonian."""
import itertools, time, sys
import numpy as np
sys.path.insert(0, '.')
from brenes_study import setup, sector, diag_ensemble, gibbs, theta, w
T0 = time.time(); rng = np.random.default_rng(3)
part = sys.argv[1] if len(sys.argv) > 1 else 'ABC'

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

if 'B' in part:
    print('(B) large-J scaling, L=8, all sectors: long-time n_pa and J^2 n_pa', flush=True)
    S = setup(8); allsec = [np.array(e) for e in itertools.product([-1, 0, 1], repeat=8)]
    for J in [5.0, 10.0, 20.0, 40.0, 80.0]:
        de = 0.0
        for eb in allsec:
            ev, V, d, E = sector(S, eb, 2 * J)
            v = np.zeros(len(ev)); v[S['v']] = 1.0
            de += diag_ensemble(ev, V, v, S['npa'])
        de /= len(allsec)
        print(f'  J={J:6.1f}: n_pa long-time {de:.5f}   J*n_pa {J * de:.4f}   J^2*n_pa {J * J * de:.3f}', flush=True)

if 'C' in part:
    print('(C) electric-energy band of the vacuum, L=8 (J large): second-order effective in-band dynamics', flush=True)
    S = setup(8); allsec = [np.array(e) for e in itertools.product([-1, 0, 1], repeat=8)]
    K = S['K']; npa = S['npa']; v0 = S['v']
    J = 50.0; g2 = 2 * J
    stats = dict(band=0.0, inband_inf=0.0, eff_de=0.0, iso=0.0, exact_de=0.0, comp=0.0)
    for eb in allsec:
        E = eb + S['F']; Nel = ((E + theta) ** 2).sum(1)           # electric energy in units of g^2/2 (integers)
        B = np.where(np.abs(Nel - Nel[v0]) < 1e-9)[0]; C = np.setdiff1d(np.arange(len(Nel)), B)
        stats['band'] += len(B); stats['inband_inf'] += npa[B].mean()
        # second-order effective Hamiltonian on the band (energies of H0 = (g^2/2) Nel, band energy (g^2/2)Nel_v0)
        E0 = 0.5 * g2 * Nel[v0]; Ec = 0.5 * g2 * Nel[C]
        KBC = K[np.ix_(B, C)]
        Heff = KBC @ np.diag(1.0 / (E0 - Ec)) @ KBC.T + K[np.ix_(B, B)]
        ev, V = np.linalg.eigh(Heff)
        i0 = int(np.where(B == v0)[0][0]); e = np.zeros(len(B)); e[i0] = 1.0
        stats['eff_de'] += diag_ensemble(ev, V, e, npa[B])
        # connected component of the vacuum under Heff
        A = (np.abs(Heff) > 1e-12).astype(int); seen = {i0}; fr = [i0]
        while fr:
            k = fr.pop(); nb = set(np.where(A[k])[0]) - seen; seen |= nb; fr += list(nb)
        stats['comp'] += len(seen); stats['iso'] += (len(seen) == 1)
        ev2, V2, d2, _ = sector(S, eb, g2); vv = np.zeros(len(ev2)); vv[v0] = 1.0
        stats['exact_de'] += diag_ensemble(ev2, V2, vv, npa)
    n = len(allsec)
    print(f'  J={J}: mean band size {stats["band"] / n:.2f} (of {len(npa)}), in-band infinite-T n_pa {stats["inband_inf"] / n:.4f}, '
          f'connected component of the vacuum under H_eff: mean size {stats["comp"] / n:.2f}, isolated in {100 * stats["iso"] / n:.1f}% of sectors; '
          f'long-time n_pa: H_eff {stats["eff_de"] / n:.4f}, exact {stats["exact_de"] / n:.5f}', flush=True)
print(f'[time {time.time() - T0:.0f}s]')
