"""Reprocess saved traces at omega1=omega2 for several causal exponential windows."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/codex_mpl_op')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'results/scan_x_float64'
a=np.load(out/'spectra.npz');w=a['omega1'];wp=w[w>=0]
n=140001;dt=.01;t=-700+dt*np.arange(n)
etas=np.array([0.,.01,.02,.05])
windows=np.exp(-etas[:,None]*np.maximum(t,0)[None,:]);windows[:,[0,-1]]*=.5
z=np.empty((len(etas),len(wp)),complex)
def delta_diff(i,amp,kind):
 arrays=[]
 for sign in ['plus','minus']:
  p=out/f'omega_{i:03d}_F{amp}_{kind}_{sign}'/'observables.f64'
  arrays.append(np.memmap(p,dtype='<f8',mode='r',shape=(n,10),order='F')[:,3].copy())
 return 4*(arrays[0]-arrays[1])/2
for i,omega in enumerate(wp):
 lc=delta_diff(i,'0.001','cos');hc=delta_diff(i,'0.01','cos')
 ls=delta_diff(i,'0.001','sin');hs=delta_diff(i,'0.01','sin')
 q=((hc-10*lc)-1j*(hs-10*ls))/.9
 fhat=.01*.1/np.sqrt(2*np.pi)*np.exp(-omega**2*.1**2/2)
 # chi = chi_bar/2. All four windows use identical dynamics and normalization.
 z[:,i]=dt*(windows@(q*np.exp(2j*omega*t)))*np.pi/(.001*.01*fhat)
assert np.allclose(z[0],np.diag(a['chi'])[w>=0],rtol=1e-7,atol=1e-8)
assert np.allclose(z[1],np.diag(a['chi_eta001'])[w>=0],rtol=1e-7,atol=1e-8)
np.savez_compressed(out/'diagonal_window_comparison.npz',omega=wp,etas=etas,chi=z)
np.savetxt(out/'diagonal_window_comparison.tsv',np.column_stack([wp,abs(z).T]),header='omega abs_chi_eta0 abs_chi_eta001 abs_chi_eta002 abs_chi_eta005',fmt='%.17e')
fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
colors=['.65','tab:blue','tab:orange','tab:red']
for ax in axs:
 for j,eta in enumerate(etas):
  if eta==0:continue
  ax.plot(wp,abs(z[j]),color=colors[j],lw=1.6,label=rf'$\eta={eta:g}$ ($\tau={1/eta:g}$)')
 ax.axvline(1.3813078088431883,color='.45',ls='--',lw=1,label=r'$\Delta_{\rm eq}$')
 ax.set(xlabel=r'$\omega_1=\omega_2=\omega$',ylabel=r'$|\chi^{(2)}(\omega,\omega)|$',ylim=(0,None));ax.grid(alpha=.2)
axs[0].set(xlim=(0,5),title='Full frequency range');axs[0].legend(fontsize=9)
axs[1].set(xlim=(1.1,2.2),title='Near the gap / continuum edge')
fig.suptitle(r'Stronger postprocessing windows: $W_\eta(t)=e^{-\eta\max(t,0)}$, $\chi^{(2)}=\bar\chi^{(2)}/2$')
fig.savefig(out/'chi2_diagonal_windows.png',dpi=180);fig.savefig(out/'chi2_diagonal_windows.svg');plt.close(fig)
for eta,v in zip(etas,z):
 i=np.argmax(abs(v));print(f'eta={eta:g}: peak at omega={wp[i]:.3f}, abs_chi={abs(v[i]):.8g}')
