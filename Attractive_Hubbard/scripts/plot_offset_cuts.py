"""Plot chi(w2+offset,w2) from the gamma=0.05 map without interpolation."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/codex_mpl_op')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parents[1]/'results/scan_x_float64'
a=np.load(out/'chi2_map_gamma0p05.npz');w1=a['omega1'];w2=a['omega2'];z=a['chi']
offsets=np.array([-2.,-1.,0.,1.,2.]);step=w1[1]-w1[0]
fig,axs=plt.subplots(3,1,figsize=(10,10),sharex=True,layout='constrained')
colors=['#3164ce','#269a68','#222222','#e29416','#d73545']
saved={}
for offset,color in zip(offsets,colors):
 shift=round(offset/step); assert np.isclose(shift*step,offset)
 j=np.arange(len(w2));i=j+shift;mask=(i>=0)&(i<len(w1));i=i[mask];j=j[mask]
 assert np.allclose(w1[i],w2[j]+offset)
 vals=z[i,j];x=w2[j]
 for ax,component in zip(axs,[np.abs,np.real,np.imag]):
  ax.plot(x,component(vals),color=color,lw=1.4,label=rf'$XX={offset:g}$')
 tag=f'{offset:+g}'.replace('+','plus').replace('-','minus').replace('.','p')
 np.savetxt(out/f'offset_cut_XX_{tag}_gamma0p05.tsv',np.c_[x,w1[i],vals.real,vals.imag,abs(vals)],header='omega2 omega1 Re_chi Im_chi abs_chi; omega1=omega2+XX; gamma=0.05; chi=chi_bar/2',fmt='%.17e')
 saved[tag]=np.c_[x,w1[i],vals.real,vals.imag]
for ax,label in zip(axs,[r'$|\chi^{(2)}|$',r'$\mathrm{Re}\,\chi^{(2)}$',r'$\mathrm{Im}\,\chi^{(2)}$']):
 ax.set_ylabel(label);ax.grid(alpha=.2);ax.axhline(0,color='.5',lw=.5);ax.set_xlim(-5,5)
axs[0].legend(fontsize=10,loc='upper left',bbox_to_anchor=(1.01,1));axs[0].set_ylim(bottom=0)
axs[-1].set_xlabel(r'$\omega_2\quad(\omega_1=\omega_2+XX)$')
fig.suptitle(r'Offset cuts: $\chi^{(2)}(\omega_2+XX,\omega_2)$, window $\Gamma=0.05$'+'\n'+r'$\chi^{(2)}=\bar\chi^{(2)}/2$; only points within $[-5,5]^2$')
fig.savefig(out/'chi2_offset_cuts_gamma0p05.png',dpi=180);fig.savefig(out/'chi2_offset_cuts_gamma0p05.svg');plt.close(fig)
np.savez_compressed(out/'chi2_offset_cuts_gamma0p05.npz',offsets=offsets,**saved)
print('Saved 5 offset cuts; no frequency interpolation or new time evolution.')
