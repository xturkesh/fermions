"""Sanity checks of the scalar inequalities used in T3.

(1) h(nu) >= 4 ln2 nu(1-nu)               (binary entropy, nats)
(2) x <= -ln(1-x) <= 2 ln2 x on [0,1/2]    (Renyi-2 vs variance)
(3) integer max-entropy bound H <= 1/2 ln(2 pi e (Var+1/12)) on random integer laws
(4) Poisson-binomial: 1/2 ln(8 Var/pi) <= H_inf <= H <= 1/2 ln(2 pi e (Var+1/12))
(5) Poisson-binomial: H - 1/2 ln(2 pi e Var) -> 0 as Var -> infinity (worst cases)
(6) max_k p_k <= sqrt(pi/8)/sigma
"""
import numpy as np

rng = np.random.default_rng(1)
ln2 = np.log(2.0)

def h(nu):
    nu = np.clip(nu, 1e-300, 1 - 1e-16)
    return -nu*np.log(nu) - (1-nu)*np.log1p(-nu)

nu = np.linspace(1e-12, 1-1e-12, 2_000_001)
r = h(nu) / (nu*(1-nu))
print("(1) min h/(nu(1-nu)) =", r.min(), " 4ln2 =", 4*ln2, " argmin =", nu[r.argmin()])

x = np.linspace(1e-12, 0.5, 1_000_001)
y = -np.log1p(-x)
print("(2) min(-ln(1-x)/x) =", (y/x).min(), " max =", (y/x).max(), " 2ln2 =", 2*ln2)

def H(p):
    p = p[p > 0]
    return -(p*np.log(p)).sum()

worst = -np.inf
for trial in range(20000):
    K = rng.integers(1, 40)
    p = rng.random(K)**rng.uniform(0.2, 8)
    if rng.random() < 0.3:          # sparse supports
        p *= (rng.random(K) < 0.3)
        if p.sum() == 0: p[0] = 1
    p /= p.sum()
    k = np.arange(K)
    var = (p*k**2).sum() - ((p*k).sum())**2
    worst = max(worst, H(p) - 0.5*np.log(2*np.pi*np.e*(var + 1/12)))
print("(3) max over random integer laws of H - 1/2 ln(2 pi e (Var+1/12)) =", worst, "(must be <= 0)")

def poisson_binomial(nus):
    p = np.array([1.0])
    for v in nus:
        p = np.convolve(p, [1-v, v])
    return p

worst_lo = np.inf; worst_hi = -np.inf; worst_max = -np.inf
for trial in range(5000):
    n = rng.integers(1, 200)
    mode = rng.integers(0, 3)
    if mode == 0:
        nus = rng.random(n)
    elif mode == 1:
        nus = rng.random(n)**6           # many small nu
    else:
        nus = np.where(rng.random(n) < 0.5, 0.5, rng.random(n)*1e-3)
    p = poisson_binomial(nus)
    var = (nus*(1-nus)).sum()
    if var < 1e-6: continue
    Hp = H(p)
    Hinf = -np.log(p.max())
    lo = 0.5*np.log(8*var/np.pi)
    hi = 0.5*np.log(2*np.pi*np.e*(var + 1/12))
    worst_lo = min(worst_lo, Hinf - lo)
    worst_hi = max(worst_hi, Hp - hi)
    worst_max = max(worst_max, p.max()*np.sqrt(var)/np.sqrt(np.pi/8))
print("(4) min(H_inf - 1/2 ln(8Var/pi)) =", worst_lo, "(>=0);  max(H - upper) =", worst_hi, "(<=0)")
print("(6) max of p_max*sigma/sqrt(pi/8) =", worst_max, "(<=1)")

print("(5) asymptotics: H - 1/2 ln(2 pi e Var) for growing Var (worst over shapes)")
for n in [10, 100, 1000, 10000]:
    vals = []
    for shape in ["uniform", "half", "small", "mixed"]:
        if shape == "uniform": nus = rng.random(n)
        elif shape == "half": nus = np.full(n, 0.5)
        elif shape == "small": nus = np.full(n, 0.05)
        else: nus = np.concatenate([np.full(n//2, 0.5), rng.random(n - n//2)**4])
        p = poisson_binomial(nus)
        var = (nus*(1-nus)).sum()
        vals.append((shape, round(var, 2), H(p) - 0.5*np.log(2*np.pi*np.e*var)))
    print("   n =", n, vals)
