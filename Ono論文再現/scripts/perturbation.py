"""Supplement S40,S41,S43,S46,S55, spinless Rice-Mele, no RTA source expansion."""
import numpy as np
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.])
class Perturbation:
 def __init__(self,nk=800,hx=1.,hy=.5,hz=.1,gamma=.01):
  k=-np.pi+2*np.pi*np.arange(nk)/nk
  s=np.sin(k/2);c=np.cos(k/2)
  H=-hx*c[:,None,None]*SX-hy*s[:,None,None]*SY-hz*SZ
  self.E,U=np.linalg.eigh(H)
  def band(Q):return U.conj().swapaxes(-1,-2)@Q@U
  self.q0=band(hx/2*s[:,None,None]*SX-hy/2*c[:,None,None]*SY)
  self.q1=band(-hx/4*c[:,None,None]*SX-hy/4*s[:,None,None]*SY)
  self.q2=band(-hx/16*s[:,None,None]*SX+hy/16*c[:,None,None]*SY)
  self.n=np.array([1.,0.]);self.gamma=gamma;self.nk=nk
 def c1(self,A,B,w,broad=None):
  w=np.atleast_1d(w);z=w+1j*(self.gamma if broad is None else broad)
  v=np.zeros(len(w),complex)
  for m in range(2):
   for n in range(2):
    numerator=A[:,m,n]*B[:,n,m]*(self.n[n]-self.n[m])
    v+=np.sum(numerator[:,None]/(z[None,:]-self.E[:,n,None]+self.E[:,m,None]),axis=0)/self.nk
  return v
 def linear(self,w):
  return np.mean(self.q1[:,0,0])+self.c1(self.q0,self.q0,w)
 def chi2(self,w1,w2,sum_gamma=1):
  w1,w2=np.broadcast_arrays(w1,w2);shape=w1.shape
  a=w1.ravel();b=w2.ravel();result=np.empty(a.size,complex)
  for start in range(0,a.size,256):
   aa=a[start:start+256];bb=b[start:start+256];om=aa+bb
   outer=om+1j*self.gamma*sum_gamma
   val=np.zeros(len(aa),complex)
   # Sum both orderings directly; same Q0 for the two inputs.
   for mu in range(2):
    for nu in range(2):
     d0=outer[None,:]-self.E[:,nu,None]+self.E[:,mu,None]
     for lam in range(2):
      prod=self.q0[:,mu,nu]*self.q0[:,nu,lam]*self.q0[:,lam,mu]
      n1=prod*(self.n[mu]-self.n[lam]);n2=prod*(self.n[nu]-self.n[lam])
      for w in (aa,bb):
       z=w+1j*self.gamma
       d1=z[None,:]-self.E[:,lam,None]+self.E[:,mu,None]
       d2=z[None,:]-self.E[:,nu,None]+self.E[:,lam,None]
       val+=np.sum((n1[:,None]/d1+n2[:,None]/d2)/d0,axis=0)/self.nk
   val+=2*np.mean(self.q2[:,0,0])
   val+=self.c1(self.q0,self.q1,om,broad=sum_gamma*self.gamma)
   val+=self.c1(self.q1,self.q0,aa)+self.c1(self.q1,self.q0,bb)
   result[start:start+len(aa)]=val
  return result.reshape(shape)
