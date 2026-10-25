#!/usr/bin/env python3
"""Extract central mixed variations and symmetrized second-order response.
FT convention: integral dt exp(+i*w*t) Q(t). No damping by default.
"""
import argparse
import json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('jobs',type=Path)
p.add_argument('--omega2',nargs=3,type=float,default=[-5,5,201],metavar=('MIN','MAX','COUNT'))
p.add_argument('--window-tau',type=float,default=None)
a=p.parse_args()
if a.window_tau is not None and a.window_tau<=0:p.error('window-tau must be positive')
w2=np.linspace(a.omega2[0],a.omega2[1],int(a.omega2[2]))
records={}
for d in sorted(a.jobs.iterdir()):
 if not (d/'metadata.json').exists():continue
 m=json.loads((d/'metadata.json').read_text())
 x=np.loadtxt(d/'obs_t.dat')
 t=x[:,0]+m['tmin']
 # Delta positive real in the initial gauge: -U_el * Phi_sc.
 q=-m['U_el']*(x[:,3]+1j*x[:,4])
 if not np.all(np.isfinite(q)):raise ValueError(f'nonfinite output in {d}')
 records[(m['F0'],m['omega'],m['kind'])]=(t,q,m)
for f0,omega,kind in list(records):
 if kind!='cos_plus':continue
 t,qp,m=records[(f0,omega,kind)]
 _,qm,_=records[(f0,omega,'cos_minus')]
 _,sp,_=records[(f0,omega,'sin_plus')]
 _,sm,_=records[(f0,omega,'sin_minus')]
 dc=(qp-qm)/2
 ds=(sp-sm)/2
 # dc includes epsilon, matching the variation in Ono's Eq.14.
 window=1 if a.window_tau is None else np.exp(-t*t/(2*a.window_tau**2))
 weights=np.ones(t.size)*m['dt'];weights[[0,-1]]*=0.5
 out=[]
 for w in w2:
  phase=np.exp(1j*(omega+w)*t)*weights*window
  mixed=np.sum(phase*(dc-1j*ds))
  fhat=f0*m['tau_f']/np.sqrt(2*np.pi)*np.exp(-w*w*m['tau_f']**2/2)
  out.append(2*np.pi*mixed/(m['epsilon']*f0*fhat))
 out=np.asarray(out)
 name=a.jobs/f'chi_bar_F{f0:.8g}_w1_{omega:.8g}.dat'
 np.savetxt(name,np.column_stack([w2,out.real,out.imag,abs(out)]),header='omega2 Re_chi_bar Im_chi_bar abs_chi_bar; central epsilon difference; assumes vanishing linear response')
 np.savetxt(a.jobs/f'variation_F{f0:.8g}_w1_{omega:.8g}.dat',np.column_stack([t,dc.real,dc.imag,ds.real,ds.imag]),header='time Re_delta_cos Im_delta_cos Re_delta_sin Im_delta_sin')
 print(name)
