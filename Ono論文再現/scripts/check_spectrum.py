from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from perturbation import Perturbation
p=Perturbation(nk=32)
w1=.43;w2=.7
# Independently evaluate the density-matrix resolvent expansion to quadratic order.
def ordered(a,b):
 total=0j
 for i in range(p.nk):
  e=p.E[i];q=p.q0[i];q1=p.q1[i];q2=p.q2[i];r0=np.diag(p.n)
  diff=e[:,None]-e[None,:]
  def rho1(B,w):return -(B@r0-r0@B)/(w+1j*p.gamma-diff)
  r1=rho1(q,b)
  r2=-(q@r1-r1@q)/(a+b+1j*p.gamma-diff)
  total+=np.trace(q@r2+.5*q@rho1(q1,a+b)+q1@r1+q2@r0)
 return total/p.nk
assert np.allclose(p.chi2(w1,w2),ordered(w1,w2)+ordered(w2,w1),atol=1e-14)
assert np.allclose(p.chi2(w1,w2),p.chi2(w2,w1),atol=1e-14)
assert np.allclose(p.chi2(-w1,-w2),p.chi2(w1,w2).conj(),atol=1e-14)
assert abs(Perturbation(nk=32,hz=0).chi2(w1,w2))<1e-14
assert abs(Perturbation(nk=32,hy=0).chi2(w1,w2))<1e-14
print('5 perturbation checks passed')
