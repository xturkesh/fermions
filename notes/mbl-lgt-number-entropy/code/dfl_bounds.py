"""
Sector-averaged free-fermion quantities for the 1D Z2 LGT (open chain), h_q = -J Lap + 2h diag(q), q iid +-1.
(1) eigenfunction correlator  E Q(x,y),  Q = sum_a |phi_a(x) phi_a(y)|  >= sup_t |U_t(x,y)|   (bulk vs boundary, several N)
(2) time-averaged CDW imbalance  E Ibar = E (1/N) sum_a a_a^2 ; identity E_a a_a = 2h sum_x q_x(-1)^x phi_a(x)^2
(3) half-chain number entropy (mixture over sectors and sector-averaged), number variance, Gaussian EE, vs t and N
(4) binary entropy inequality h(nu) <= 2 ln2 sqrt(nu(1-nu)),  h(nu) >= 4 ln2 nu(1-nu)
"""
import numpy as np, sys
from scipy.linalg import eigh_tridiagonal

def eig(N, J, h, q):
    return eigh_tridiagonal(2 * h * q.astype(float), -J * np.ones(N - 1))

def hbin(nu):
    nu = np.clip(nu, 1e-300, 1 - 1e-16)
    return -(nu * np.log(nu) + (1 - nu) * np.log1p(-nu))

def part4():
    nu = np.linspace(1e-9, 1 - 1e-9, 200001)
    r1 = hbin(nu) / (2 * np.log(2) * np.sqrt(nu * (1 - nu)))
    r2 = hbin(nu) / (4 * np.log(2) * nu * (1 - nu))
    print(f'(4) max h/(2ln2 sqrt(nu(1-nu))) = {r1.max():.6f} (<=1),  min h/(4ln2 nu(1-nu)) = {r2.min():.6f} (>=1)')

def part1(J=1.0, samples=400, rng=None):
    print('(1) E Q(x,y) [Q >= sup_t|U_t(x,y)|], x = bulk center or boundary site 0')
    for h in [0.25, 0.5, 1.0, 2.0]:
        for N in [100, 200, 400]:
            Qc = np.zeros(N); Qb = np.zeros(N)
            for s in range(samples):
                q = rng.choice([-1, 1], size=N)
                E, V = eig(N, J, h, q)
                A = np.abs(V)
                Qc += A @ A[N // 2] / samples
                Qb += A @ A[0] / samples
            d = np.arange(1, 41)
            # fit decay rate on distances 5..40
            sl = slice(5, 40)
            etac = -np.polyfit(d[sl], np.log(Qc[N // 2 + d[sl]]), 1)[0]
            etab = -np.polyfit(d[sl], np.log(Qb[d[sl]]), 1)[0]
            print(f'   lambda=2h/J={2*h:4.1f} N={N}: bulk EQ(d=1,10,20,40)=' +
                  ' '.join(f'{Qc[N//2+k]:.2e}' for k in (1, 10, 20, 40)) +
                  f' rate={etac:.3f} | boundary EQ(d=1,10,20,40)=' + ' '.join(f'{Qb[k]:.2e}' for k in (1, 10, 20, 40)) + f' rate={etab:.3f}')
        sys.stdout.flush()

def part2(J=1.0, samples=400, rng=None):
    print('(2) time-averaged CDW imbalance E Ibar = E (1/N) sum a_alpha^2  (exact, simple spectrum)')
    maxid = 0.0
    for h in [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]:
        row = []
        for N in [50, 100, 200, 400]:
            acc = 0.0
            for s in range(samples):
                q = rng.choice([-1, 1], size=N)
                E, V = eig(N, J, h, q)
                st = (-1.0) ** np.arange(N)
                a = st @ (V ** 2)
                b = (q * st) @ (V ** 2)
                maxid = max(maxid, np.abs(E * a - 2 * h * b).max())
                acc += (a ** 2).sum() / N / samples
            row.append(acc)
        print(f'   lambda={2*h:4.2f}: E Ibar(N=50,100,200,400) = ' + ' '.join(f'{x:.4f}' for x in row))
    print(f'   identity  E_a a_a = 2h sum_x q_x (-1)^x phi_a(x)^2 : max violation {maxid:.1e}')

def poisson_binomial(nu):
    p = np.array([1.0])
    for v in nu:
        p = np.convolve(p, [1 - v, v])
    return p

def part3(J=1.0, samples=300, rng=None):
    print('(3) half chain, CDW: Gaussian EE (sector avg), number variance, S_N of mixture and sector-avg S_N, vs t')
    times = [0.0, 1.0, 10.0, 1e2, 1e3, 1e4, 1e6, 1e8]
    for h in [0.5, 1.0]:
        for N in [64, 128, 256]:
            l = N // 2
            SG = np.zeros(len(times)); VarS = np.zeros(len(times)); SNs = np.zeros(len(times))
            pmix = [np.zeros(l + 1 + 1) for _ in times]
            occ = np.arange(0, N, 2)
            for s in range(samples):
                q = rng.choice([-1, 1], size=N)
                E, V = eig(N, J, h, q)
                for it, t in enumerate(times):
                    U = (V * np.exp(-1j * E * t)) @ V.T
                    Phi = U[:, occ]
                    CA = Phi[:l] @ Phi[:l].conj().T
                    nu = np.clip(np.linalg.eigvalsh(CA), 0, 1)
                    SG[it] += hbin(nu).sum() / samples
                    VarS[it] += (nu * (1 - nu)).sum() / samples
                    p = poisson_binomial(nu)
                    pp = p[p > 1e-300]
                    SNs[it] += -(pp * np.log(pp)).sum() / samples
                    pmix[it][:len(p)] += p / samples
            print(f'   lambda={2*h}, N={N}:')
            for it, t in enumerate(times):
                pm = pmix[it]; pm = pm[pm > 1e-300]
                k = np.arange(len(pmix[it]))
                mean = (k * pmix[it]).sum(); varmix = (k ** 2 * pmix[it]).sum() - mean ** 2
                Smix = -(pm * np.log(pm)).sum()
                bound = 0.5 * np.log(2 * np.pi * np.e * (varmix + 1 / 12))
                print(f'     t={t:8.0e}: S_G={SG[it]:.3f}  E Var_q N_A={VarS[it]:.3f}  Var_mix N_A={varmix:.3f}  '
                      f'E_q S_N^q={SNs[it]:.3f}  S_N(mixture)={Smix:.3f}  (Massey bound {bound:.3f})')
            sys.stdout.flush()

if __name__ == '__main__':
    rng = np.random.default_rng(12345)
    part4()
    part2(rng=rng)
    part1(rng=rng)
    part3(rng=rng)
