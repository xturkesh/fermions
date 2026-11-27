# Temporal IPR of measurement records in model 1 (Haar collision model): S (d) + E (D), Haar on S⊗E each step,
# global depolarizing p on E after each step, S measured in the computational basis after each step
# (measure-and-continue). IPR = sum_x p(x)^2, x = (o_1..o_n). Exact 2-replica transfer vs Monte Carlo.
import numpy as np, itertools as it
rng = np.random.default_rng(21)
def haar(N):
    Z = (rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))/np.sqrt(2)
    Q, R = np.linalg.qr(Z); return Q*(np.diag(R)/np.abs(np.diag(R)))
# ---------------- exact transfer ----------------
def exact(d, D, p, n, PE, recorded):
    """recorded: list of bools per step. Returns E sum_x p(x)^2 for the recorded steps (others marginalized)."""
    q = d*D; w0, w1 = 1/(q*q-1), -1/(q*(q*q-1)); u = (1-p)**2
    cI, cF = (1-u)/D, u                                   # global depolarizing on E in moment space
    ovS = {'pure': (1., 1.), 'I': (d*d, d), 'Delta': (d, d)}     # Tr(tau_S X) for tau=(e,F)
    comps = [('pure', np.array([1.0, PE]))]               # list of (S-input type, E moments (m_I, m_F))
    for t, rec in enumerate(recorded):
        a = np.zeros(2)                                   # Haar twirl coefficients a_e, a_F
        for (st, m) in comps:
            o = ovS[st]
            b = np.array([o[0]*m[0], o[1]*m[1]])          # Tr(tau_{SE} X) for tau=e,F
            a += np.array([w0*b[0] + w1*b[1], w1*b[0] + w0*b[1]])
        E_e, E_F = np.array([D*D, D]), np.array([D, D*D])  # moments of e_E and F_E
        if t == n-1:                                      # final: contract S output, trace E per replica
            cS = (d, d) if rec else (d*d, d)
            return a[0]*cS[0]*D*D + a[1]*cS[1]*D
        if rec:
            m = d*(a[0]*E_e + a[1]*E_F); comps = [('Delta', m)]          # sum_o <oo|sigma|oo> = 1, d outcomes... weight via Delta
            comps = [('Delta', (a[0]*E_e + a[1]*E_F))]
        else:
            comps = [('I', a[0]*E_e), ('Delta', a[1]*E_F)]
        comps = [(st, np.array([m[0], cI*m[0] + cF*m[1]])) for (st, m) in comps]   # noise on E
def ipr_exact(d, D, p, n, PE): return exact(d, D, p, n, PE, [True]*n)
def lam(d, D): return (D+1)/(d*D+1), (D-1)/(d*D-1)
# ---------------- Monte Carlo ----------------
def records(d, D, p, n, pureE):
    rhoE = np.zeros((D, D), complex)
    if pureE: rhoE[0, 0] = 1
    else: rhoE = np.eye(D)/D
    s0 = np.zeros((d, d)); s0[0, 0] = 1
    branches = {(): np.kron(s0, rhoE)}
    for t in range(n):
        U = haar(d*D); new = {}
        for x, r in branches.items():
            r = U @ r @ U.conj().T
            if p > 0:                                      # depolarize E
                T = r.reshape(d, D, d, D); trE = np.einsum('aebe->ab', T)
                r = (1-p)*r + p*np.kron(trE, np.eye(D)/D)
            for o in range(d):
                P = np.zeros((d, d)); P[o, o] = 1; PP = np.kron(P, np.eye(D))
                new[x + (o,)] = PP @ r @ PP
        branches = new
    return {x: np.real(np.trace(r)) for x, r in branches.items()}
def mc(d, D, p, n, pureE, S, n1=None):
    out = []
    for _ in range(S):
        pr = records(d, D, p, n, pureE); ipr = sum(v*v for v in pr.values())
        if n1 is None: out.append([ipr]); continue
        past, fut = {}, {}
        for x, v in pr.items():
            past[x[:n1]] = past.get(x[:n1], 0) + v; fut[x[n1:]] = fut.get(x[n1:], 0) + v
        out.append([ipr, sum(v*v for v in past.values()), sum(v*v for v in fut.values())])
    out = np.array(out); return out.mean(0), out.std(0)/np.sqrt(S)
if __name__ == "__main__":
    d = 2
    print("closed form check (noiseless): IPR_n = (1+PE)/2 lam+^n + (1-PE)/2 lam-^n, lam = (D+-1)/(dD+-1)")
    for (D, pureE, p) in [(2, True, 0), (3, False, 0), (2, False, 0.3), (3, True, 0.3)]:
        PE = 1.0 if pureE else 1/D
        for n in [1, 2, 4]:
            m, e = mc(d, D, p, n, pureE, 1500)
            lp, lm = lam(d, D); cf = (1+PE)/2*lp**n + (1-PE)/2*lm**n if p == 0 else float('nan')
            print("  D=%d %s p=%.1f n=%d  MC=%.5f+-%.5f  transfer=%.5f  closed=%.5f" % (D, 'pure ' if pureE else 'mixed', p, n, m[0], e[0], ipr_exact(d, D, p, n, PE), cf))
    print("past/future marginal IPRs (n1=2, n=4):")
    for (D, pureE, p) in [(2, True, 0), (3, True, 0.3)]:
        PE = 1.0 if pureE else 1/D
        m, e = mc(d, D, p, 4, pureE, 1500, n1=2)
        ex = [exact(d, D, p, 4, PE, [True]*4), exact(d, D, p, 4, PE, [True, True, False, False]), exact(d, D, p, 4, PE, [False, False, True, True])]
        print("  D=%d p=%.1f  MC (all,past,fut)=%s +- %s  exact=%s" % (D, p, np.round(m, 5), np.round(e, 5), np.round(ex, 5)))
