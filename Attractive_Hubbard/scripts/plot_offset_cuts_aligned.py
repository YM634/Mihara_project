"""Align offset cuts by using the half-sum frequency as the horizontal axis."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/codex_mpl_op')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parents[1]/'results/scan_x_float64'
a=np.load(out/'chi2_offset_cuts_gamma0p05.npz')
fig,axs=plt.subplots(3,1,figsize=(10,9),sharex=True,layout='constrained')
colors=['#3164ce','#269a68','#222222','#e29416','#d73545']
styles=['-', '--', ':', '-.', '-']
curves=[]
for offset,color,style in zip(a['offsets'],colors,styles):
 tag=f'{offset:+g}'.replace('+','plus').replace('-','minus').replace('.','p')
 v=a[tag];x=(v[:,0]+v[:,1])/2;z=v[:,2]+1j*v[:,3];curves.append((offset,x,z))
 for ax,fn in zip(axs,[np.abs,np.real,np.imag]):
  ax.plot(x,fn(z),color=color,ls=style,lw=1.5,label=rf'$XX={offset:g}$')
for ax,label in zip(axs,[r'$|\chi^{(2)}|$',r'$\mathrm{Re}\,\chi^{(2)}$',r'$\mathrm{Im}\,\chi^{(2)}$']):
 ax.set_ylabel(label);ax.set_xlim(-4,4);ax.grid(alpha=.2);ax.axhline(0,color='.5',lw=.5)
axs[0].set_ylim(bottom=0);axs[0].legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=10)
axs[-1].set_xlabel(r'$\bar\omega=(\omega_1+\omega_2)/2=\omega_2+XX/2$')
fig.suptitle(r'Aligned cuts: $\chi^{(2)}(\bar\omega+XX/2,\bar\omega-XX/2)$'+'\n'+r'Window $\Gamma=0.05$; $\chi^{(2)}=\bar\chi^{(2)}/2$')
fig.savefig(out/'chi2_offset_cuts_aligned_gamma0p05.png',dpi=180);fig.savefig(out/'chi2_offset_cuts_aligned_gamma0p05.svg');plt.close(fig)
# Compare only the common half-sum grid, without interpolation.
base=next((x,z) for off,x,z in curves if off==0)
for off,x,z in curves:
 common=np.round(x,8); bgrid=np.round(base[0],8); xs,i,j=np.intersect1d(common,bgrid,return_indices=True)
 rel=np.linalg.norm(z[i]-base[1][j])/np.linalg.norm(base[1][j]);print(f'XX={off:g}: complex relative L2 vs XX=0 on shared grid = {rel:.8g}')
