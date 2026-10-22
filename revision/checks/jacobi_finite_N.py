import numpy as np
from scipy.special import roots_jacobi, eval_jacobi, gammaln
from scipy.integrate import quad
m2 = lambda x: -np.log2(1-x+x*x)
def avg(N, l, b):
    a=-0.5; xs, ws = roots_jacobi(400, a, b)   # nodes in y, weight (1-y)^a(1+y)^b
    x=(1-xs)/2
    tot=0.0
    for n in range(l):
        hn=np.exp(gammaln(n+a+1)+gammaln(n+b+1)-gammaln(n+1)-gammaln(n+a+b+2))*(n+a+b+1)/(2*n+a+b+1) if n>0 else np.exp(gammaln(a+1)+gammaln(b+1)-gammaln(a+b+2))
        tot+=np.sum(ws*eval_jacobi(n,a,b,xs)**2*m2(x))/hn/2**(a+b+1)
    return tot
def inf(r):
    c=4*r*(1-r)
    f=lambda t: c/np.pi*np.cos(t)**2/(1-c*np.sin(t)**2)*m2(c*np.sin(t)**2)
    return quad(f,0,np.pi/2)[0]
if __name__ == "__main__":
  print(" N   l   correct/N   written/N   diff/N    N->inf")
  for N,l in [(8,2),(8,4),(32,4),(32,8),(32,16),(64,32),(128,64),(1024,512)]:
      c=avg(N,l,N-2*l); w=avg(N,l,N-2*l-0.5)
      print(f"{N:4d} {l:4d}  {c/N:.5f}    {w/N:.5f}    {(w-c)/N:+.5f}   {inf(l/N):.5f}")
