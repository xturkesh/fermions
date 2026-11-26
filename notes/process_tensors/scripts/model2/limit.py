# D -> infinity two-sector form of the k=2 bottleneck OTOC, with noise u=(1-p)^2 on E (copy-wise) and optional noise on C.
import numpy as np
from bottleneck import gate, recursion, perm_C, dep_on_pair, X, Y, Z, I2
def maps(U, d, pC=0.0):
    Uc = U.conj().T; U4 = U.reshape(d, d, d, d); Ud4 = Uc.reshape(d, d, d, d)
    def Ulayer(Xt, tau):
        Pc = perm_C(d, tau)
        full = np.einsum('ikjl,mnop->imknjolp', Xt, Pc).reshape([d]*8)
        full = np.einsum('abxy,xyzwpqrs->abzwpqrs', Ud4, full)
        full = np.einsum('cdzw,abzwpqrs->abcdpqrs', Ud4, full)
        full = np.einsum('abcdxyrs,xypq->abcdpqrs', full, U4)
        full = np.einsum('abcdpqxy,xyrs->abcdpqrs', full, U4)
        return dep_on_pair(full.reshape(d**4, d**4), pC, 'C', d)
    def R(Yop, sig):
        return np.einsum('aibjckdl,klij->abcd', Yop.reshape([d]*8), perm_C(d, sig))
    return Ulayer, R
def otoc_inf(U, A, B, t, u=1.0, pC=0.0):
    d = A.shape[0]; Ulayer, R = maps(U, d, pC)
    Xe = np.einsum('ij,kl->ikjl', A, A); XF = np.zeros_like(Xe)
    BB = np.kron(np.kron(B, np.eye(d)), np.kron(B, np.eye(d)))
    Fac = np.eye(d**4).reshape(d*d, d*d, d*d, d*d).transpose(1, 0, 2, 3).reshape(d**4, d**4)
    for step in range(t):
        Ye, YF = Ulayer(Xe, 0), Ulayer(XF, 1)
        if step == t-1: break
        Xe, XF = R(Ye, 0)/d**2, u*R(YF, 1)/d**2 + R(Ye, 1)/d**2 - R(Ye, 0)/d**3
    return np.real(np.trace(Ye @ BB @ Fac) + np.trace(YF @ BB @ Fac))/d**2
def single(U, A, B, t):                                       # two-point function C1(t) = Tr(A(t)B)/(d^2 D)
    d = A.shape[0]; Xa = A.copy()
    for step in range(t):
        Yop = U.conj().T @ np.kron(Xa, np.eye(d)) @ U
        if step == t-1: break
        Xa = np.einsum('aibi->ab', Yop.reshape(d, d, d, d))/d
    return np.real(np.trace(Yop @ np.kron(B, np.eye(d))))/d**2
if __name__ == "__main__":
    U = gate(0.5, 0.3); A, B = Z, Z
    print("check D->inf limit vs exact recursion at large D:")
    for t in [2, 3, 5]:
        print("  t=%d  inf=%.6f  D=1e3: %.6f  D=1e5: %.6f   | noisy pE=0.2: inf=%.6f D=1e5: %.6f" % (
            t, otoc_inf(U, A, B, t), recursion(U, A, B, t, 1000), recursion(U, A, B, t, 100000),
            otoc_inf(U, A, B, t, u=0.8**2), recursion(U, A, B, t, 100000, pE=0.2)))
