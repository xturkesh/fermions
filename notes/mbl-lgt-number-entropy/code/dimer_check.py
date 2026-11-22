"""(i) ED: the gauge initial state [odd links |up>_z, even links in sigma^x eigenstates x_{2k}=(-1)^k] (x) CDW
    produces sector weights supported on q_{2k-1}=q_{2k} (random-dimer law), uniform over pairs.
(ii) free-fermion numerics: single particle, random-dimer potential +-2h (pairs) vs iid Bernoulli, h=0.25 J (<J/2):
    <x^2>(t) averaged over disorder.  Random dimer theory (Dunlap-Wu-Phillips; Jitomirskaya-Schulz-Baldes-Stolz):
    superdiffusive growth near critical energies E=+-2h when 4h<2J."""
import numpy as np, itertools
from scipy.linalg import eigh_tridiagonal
from lgt1d_ed import build
# (i)
N = 6; J = 1.0; h = 0.25
m = build(N, J, h, 'open_bl')
# gauge state: link b (0-indexed b=0..N-1 is link (b+1,b+2) 1-indexed); 1-indexed link j=b+1
up = np.array([1, 0.]); plus = np.array([1, 1.]) / np.sqrt(2); minus = np.array([1, -1.]) / np.sqrt(2)
gs = np.array([1.])
for b in range(N):
    j = b + 1
    if j % 2 == 1: v = up
    else: v = plus if (-1) ** (j // 2) == 1 else minus
    gs = np.kron(gs, v)
occ = [1, 3, 5]   # 0-indexed sites 1,3,5 = 1-indexed even sites 2,4,6 (n_{2k}=1, n_{2k-1}=0)
vf = np.zeros(2 ** N); vf[int(''.join('1' if j in occ else '0' for j in range(N)), 2)] = 1
psi0 = np.kron(vf, gs)
Qs = [m['Q'][j].toarray() for j in range(N)]
tot = 0
for q in itertools.product([1, -1], repeat=N):
    P = np.eye(len(psi0))
    for j in range(N): P = P @ (np.eye(len(psi0)) + q[j] * Qs[j]) / 2
    w = np.vdot(psi0, P @ psi0).real
    if w > 1e-12:
        tot += w; print('  sector q =', q, ' weight', round(w, 4))
print('  total weight', round(tot, 6))
# (ii)
rng = np.random.default_rng(3)
Nn = 3000; x0 = Nn // 2; S = 60
times = [10, 30, 100, 300, 1000, 3000]
for law in ['dimer', 'iid']:
    for hh in [0.25, 0.75]:
        m2 = np.zeros(len(times))
        for s in range(S):
            if law == 'dimer':
                q = np.repeat(rng.choice([-1, 1], size=Nn // 2), 2)
            else:
                q = rng.choice([-1, 1], size=Nn)
            E, V = eigh_tridiagonal(2 * hh * q.astype(float), -J * np.ones(Nn - 1))
            c = V[x0]
            for it, t in enumerate(times):
                u = V @ (np.exp(-1j * E * t) * c)
                m2[it] += ((np.arange(Nn) - x0) ** 2 * np.abs(u) ** 2).sum() / S
        sl = np.polyfit(np.log(times[2:]), np.log(m2[2:]), 1)[0]
        print(f'{law:5s} h/J={hh}: <x^2>(t) for t={times}: ' + ' '.join(f'{v:.1f}' for v in m2) + f'   slope(t>=100)={sl:.2f}')
