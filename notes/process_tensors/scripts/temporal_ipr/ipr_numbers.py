import numpy as np
from records_m1 import exact
def exact_c(d, D, cI, cF, recorded, PE):
    # same as records_m1.exact but with generic noise constants (local noise)
    q = d*D; w0, w1 = 1/(q*q-1), -1/(q*(q*q-1)); n = len(recorded)
    ovS = {'pure': (1., 1.), 'I': (d*d, d), 'Delta': (d, d)}
    comps = [('pure', np.array([1.0, PE]))]
    E_e, E_F = np.array([D*D, D], float), np.array([D, D*D], float)
    for t, rec in enumerate(recorded):
        a = np.zeros(2)
        for (st, m) in comps:
            o = ovS[st]; b = np.array([o[0]*m[0], o[1]*m[1]]); a += np.array([w0*b[0] + w1*b[1], w1*b[0] + w0*b[1]])
        if t == n-1:
            cS = (d, d) if rec else (d*d, d); return a[0]*cS[0]*D*D + a[1]*cS[1]*D
        comps = [('Delta', a[0]*E_e + a[1]*E_F)] if rec else [('I', a[0]*E_e), ('Delta', a[1]*E_F)]
        comps = [(st, np.array([m[0], cI*m[0] + cF*m[1]])) for (st, m) in comps]
def local_c(q, L, p):
    D = q**L; u = (1-p)**2; GI, GF = float(D), (u*q*q + 1 - u)**L
    return (D*GI - GF)/(D*(D*D-1)), (D*GF - GI)/(D*(D*D-1))
d = 2
print("(a) normalized temporal IPR  N*IPR = d^n IPR  [1=uniform, 2=Porter-Thomas]")
for PEtag in ['pure', 'mixed']:
    for D in [1, 2, 8, 32, 128]:
        row = []
        for n in [1, 4, 16, 64, 256]:
            lp = (D+1)/(d*D+1); lm = (D-1)/(d*D-1) if D > 1 else 0.0
            PE = 1.0 if PEtag == 'pure' or D == 1 else 1/D
            row.append(d**n*((1+PE)/2*lp**n + (1-PE)/2*lm**n))
        print("  %s D=%3d  n=1,4,16,64,256: %s" % (PEtag, D, "  ".join("%.4g" % v for v in row)))
print("\n(b) record memory witness dI_rec = log2[IPR_all/(IPR_past IPR_fut)]  vs quantum Renyi-2 I2 (n1=n2)")
for PEtag in ['pure', 'mixed']:
    for D in [2, 8, 32]:
        PE = 1.0 if PEtag == 'pure' else 1/D
        row = []
        for n1 in [1, 2, 4, 8, 16]:
            R = lambda rec: exact(d, D, 0.0, 2*n1, PE, rec)
            row.append(np.log2(R([True]*(2*n1))/(R([True]*n1+[False]*n1)*R([False]*n1+[True]*n1))))
        print("  %s D=%2d  n1=n2=1,2,4,8,16: %s bits   (quantum I2 -> %.0f bits)" % (PEtag, D, "  ".join("%.2e" % v for v in row), 2*np.log2(D)))
print("\n(c) local depolarizing p on each of L env qubits: record IPR growth rate per step (log2 of N*IPR ratio) and dI_rec, pure bath, n1=n2=12")
for p in [0.0, 0.1, 0.3]:
    out = []
    for L in [1, 2, 4, 6, 8, 10, 12]:
        D = 2**L; cI, cF = local_c(2, L, p)
        r = lambda rec: exact_c(d, D, cI, cF, rec, 1.0)
        g = np.log2(d**25*r([True]*25)/(d**24*r([True]*24)))
        dI = np.log2(r([True]*24)/(r([True]*12+[False]*12)*r([False]*12+[True]*12)))
        out.append("L=%d: g=%.2e dI=%.2e" % (L, g, dI))
    print("  p=%.1f  " % p + " | ".join(out))
