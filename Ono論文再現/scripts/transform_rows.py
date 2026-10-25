"""Fourier transform measured rows only; FT convention exp(+i omega t)."""
from pathlib import Path
import numpy as np
from scipy.signal import ZoomFFT
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fig2_float64_stronger'
DT=.01;TMIN=-1000.;N=200001
F=.01;EPS=.001;TAUF=.1

def load(w):
 p=OUT/f'row_{round(w*1e6):07d}.bin'
 v=np.fromfile(p,dtype=np.float64)
 if v.size!=2*N or not np.all(np.isfinite(v)):raise ValueError(f'Bad row {p}')
 return v[:N],v[N:]

def transform(w1,w2):
 dc,ds=load(w1);q=dc-1j*ds
 q[[0,-1]]*=.5
 w2=np.atleast_1d(w2);om=w1+w2
 if len(om)>1 and np.allclose(np.diff(om),np.diff(om)[0],atol=1e-12):
  f=ZoomFFT(N,[-om[0]/(2*np.pi),-om[-1]/(2*np.pi)],m=len(om),fs=1/DT,endpoint=True)(q)
  f*=DT*np.exp(1j*om*TMIN)
 else:
  t=TMIN+DT*np.arange(N)
  f=np.array([DT*np.sum(q*np.exp(1j*w*t)) for w in om])
 fhat=F*TAUF/np.sqrt(2*np.pi)*np.exp(-w2*w2*TAUF*TAUF/2)
 return 2*np.pi*f/(EPS*F*fhat)

if __name__=='__main__':
 import sys
 w=float(sys.argv[1]) if len(sys.argv)>1 else .5
 vals=transform(w,np.array([.5,1.,1.5]))
 # Independent direct quadrature verifies ZoomFFT frequency sign, origin and normalization.
 from perturbation import Perturbation
 exact=Perturbation().chi2(w,np.array([.5,1.,1.5]))
 print('omega2, measured, perturbation')
 for w2,num,ref in zip([.5,1.,1.5],vals,exact):print(w2,num,ref)
 dc,ds=load(w);t=TMIN+DT*np.arange(N);q=dc-1j*ds;q[[0,-1]]*=.5
 direct=np.array([DT*np.sum(q*np.exp(1j*(w+w2)*t)) for w2 in [.5,1.,1.5]])
 fhat=F*TAUF/np.sqrt(2*np.pi)*np.exp(-np.array([.5,1.,1.5])**2*TAUF**2/2)
 direct*=2*np.pi/(EPS*F*fhat)
 print('ZoomFFT-direct error',max(abs(vals-direct)))
 assert np.allclose(vals,direct,rtol=1e-8,atol=1e-10)
