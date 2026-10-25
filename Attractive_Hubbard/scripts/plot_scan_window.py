"""Reprocess the full chi^(2) map using a chosen postprocessing window rate."""
from pathlib import Path
import argparse,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/codex_mpl_op')
import numpy as np
from scipy.signal import ZoomFFT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
parser=argparse.ArgumentParser();parser.add_argument('--gamma',type=float,required=True);parser.add_argument('--directory',type=Path,default=Path(__file__).resolve().parents[1]/'results/scan_x_float64');args=parser.parse_args()
gamma=args.gamma
if gamma<0:raise ValueError('gamma must be nonnegative')
out=args.directory.resolve()
try: import tomllib
except ModuleNotFoundError: import tomli as tomllib
config=tomllib.loads((out/'scan_config.toml').read_text());model=config['model'];scan=config['scan']
lo,hi=scan['amplitudes'];fl,fh=float(lo),float(hi);ratio=fh/fl;eps=float(model['epsilon']);tauf=float(model['tau_f']);U=float(model['U'])
pol=np.asarray(model['polarization'],dtype=float);pol/=np.linalg.norm(pol);angle=np.degrees(np.arctan2(pol[1],pol[0]))
a=np.load(out/'spectra.npz');w=a['omega1'];wp=w[w>=0];dt=float(model['dt']);tmin=float(model['tmin']);n=round((float(model['tmax'])-tmin)/dt)+1;t=tmin+dt*np.arange(n)
window=np.exp(-gamma*np.maximum(t,0));window[[0,-1]]*=.5
z=np.empty((len(w),len(w)),complex)
def difference(i,amp,kind):
 v=[]
 for sign in ['plus','minus']:
  path=out/f'omega_{i:03d}_F{amp}_{kind}_{sign}'/'observables.f64'
  v.append(np.memmap(path,dtype='<f8',shape=(n,10),mode='r',order='F')[:,3].copy())
 return -U*(v[0]-v[1])/2
fhat=fh*tauf/np.sqrt(2*np.pi)*np.exp(-w*w*tauf**2/2)
for i,omega in enumerate(wp):
 lc=difference(i,lo,'cos');hc=difference(i,hi,'cos')
 ls=difference(i,lo,'sin');hs=difference(i,hi,'sin')
 q=((hc-ratio*lc)-1j*(hs-ratio*ls))/(1-1/ratio)*window;om=omega+w
 ft=ZoomFFT(n,[-om[0]/(2*np.pi),-om[-1]/(2*np.pi)],m=len(w),fs=1/dt,endpoint=True)(q)*dt*np.exp(1j*om*tmin)
 row=ft*np.pi/(eps*fh*fhat) # chi = chi_bar/2
 z[len(wp)-1+i]=row
 if i:z[len(wp)-1-i]=row[::-1].conj()
 if i%20==0:print('transformed',i+1,'/',len(wp),flush=True)
tag=f'{gamma:g}'.replace('.','p');stem=f'chi2_map_gamma{tag}'
np.savez_compressed(out/f'{stem}.npz',omega1=w,omega2=w,chi=z,chi_bar=2*z,window_gamma=gamma)
fig,axs=plt.subplots(1,3,figsize=(16,5),layout='constrained');halfstep=(w[1]-w[0])/2;extent=(w[0]-halfstep,w[-1]+halfstep,w[0]-halfstep,w[-1]+halfstep)
for ax,data,title in zip(axs,[abs(z),z.real,z.imag],[r'$|\chi^{(2)}|$',r'$\mathrm{Re}\,\chi^{(2)}$',r'$\mathrm{Im}\,\chi^{(2)}$']):
 opts=dict(cmap='magma') if title.startswith('$|') else dict(cmap='RdBu_r',vmin=-np.max(abs(data)),vmax=np.max(abs(data)))
 im=ax.imshow(data.T,origin='lower',extent=extent,aspect='equal',**opts)
 ax.set(xlabel=r'$\omega_1$',ylabel=r'$\omega_2$',title=title,xlim=(-5,5),ylim=(-5,5));fig.colorbar(im,ax=ax,shrink=.8)
fig.suptitle(rf'Attractive Hubbard: $\Delta$ response, polarization ${angle:g}^\circ$, $\chi^{{(2)}}=\bar\chi^{{(2)}}/2$'+'\n'+rf'Postprocessing window $W(t)=e^{{-{gamma:g}\max(t,0)}}$ (dynamics unchanged)')
fig.savefig(out/f'{stem}.png',dpi=180);fig.savefig(out/f'{stem}.svg');plt.close(fig)
print('SAVED',stem, 'max abs chi',np.max(abs(z)),flush=True)
