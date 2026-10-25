"""Create Fig.2(d,e)-style scientific plots from the actual Float64 dynamics."""
from pathlib import Path
import json,time
try: import tomllib
except ModuleNotFoundError: import tomli as tomllib
import numpy as np
from scipy.signal import ZoomFFT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LinearSegmentedColormap
CMAP=LinearSegmentedColormap.from_list("ono",["#050008","#26005a","#3512cf","#0045ff","#00c9f8","#25ee3c","#f9f32c","#ff5b0c","#8c0010"])
from perturbation import Perturbation
from transform_rows import OUT,DT,TMIN,N,F,EPS,TAUF,load,transform
params=tomllib.loads((OUT/'parameters.toml').read_text())
wp=np.array(params['map_positive_frequencies']); wd=np.array(params['diagonal_frequencies'])
ws=sorted(set(wp)|set(wd))
missing=[w for w in ws if not (OUT/f'row_{round(w*1e6):07d}.bin').exists()]
if missing:raise RuntimeError(f'{len(missing)} rows unfinished')
start=time.perf_counter()
w=np.linspace(-4,4,201)
z=np.empty((201,201),complex)
for a in wp:
 i=round((a+4)/.04)
 z[i]=transform(a,w)
 if a>0:z[round((-a+4)/.04)]=z[i,::-1].conj()
measured=np.array([transform(a,np.array([a]))[0] for a in wd])
p=Perturbation()
wf=np.linspace(0,3,1501)
ref1=p.chi2(wf,wf,sum_gamma=1);ref2=p.chi2(wf,wf,sum_gamma=2)
ref_diag=p.chi2(wd,wd,sum_gamma=1)
refcache=OUT/'reference_map.npy'
if refcache.exists():
 refmap=np.load(refcache)
else:
 refmap=p.chi2(w[:,None],w[None,:],sum_gamma=1)
 np.save(refcache,refmap)
np.savez_compressed(OUT/'spectra.npz',omega1=w,omega2=w,chi_bar=z,diagonal_omega=wd,diagonal_chi_bar=measured,reference_map=refmap,reference_diagonal=ref_diag,reference_fine_omega=wf,reference_fine_gamma1=ref1,reference_fine_gamma2=ref2)
np.savetxt(OUT/'diagonal.tsv',np.c_[wd,measured.real,measured.imag,abs(measured),ref_diag.real,ref_diag.imag],header='omega Re_chi_bar Im_chi_bar abs_chi_bar Re_reference Im_reference',fmt='%.17e')
plt.rcParams.update({'font.size':11,'axes.grid':False,'savefig.dpi':180})
fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
ax.plot(wd,measured.real,color='crimson',lw=1.5,label='Time evolution (Float64, central difference)')
ax.plot(wf,ref1.real,color='darkorange',ls='--',lw=1.1,label=r'Perturbation: $\omega_1+\omega_2+i\Gamma$')
ax.plot(wf,ref2.real,color='royalblue',ls='-.',lw=1.1,label=r'Perturbation: $\omega_1+\omega_2+2i\Gamma$')
ax.set(xlim=(0,3),ylim=(-.21,.21),xlabel=r'$\omega_1=\omega_2$',ylabel=r'$\mathrm{Re}\,\bar\chi^{(2)}(\omega_1,\omega_1)$',title='Rice-Mele second-harmonic response | Fig. 2(d) comparison')
ax.axhline(0,color='.7',lw=.5);ax.legend(loc='lower right',fontsize=9)
ins=inset_axes(ax,width='34%',height='36%',loc='upper right',borderpad=1.2)
ins.plot(wd,measured.real,color='crimson',lw=1.2)
ins.plot(wf,ref1.real,color='darkorange',ls='--',lw=1)
ins.plot(wf,ref2.real,color='royalblue',ls='-.',lw=1)
ins.set(xlim=(.48,.54),ylim=(0,.095));ins.tick_params(labelsize=8)
fig.savefig(OUT/'fig2d_second_harmonic.png');fig.savefig(OUT/'fig2d_second_harmonic.svg')
plt.close(fig)
fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
im=ax.imshow(abs(z).T,origin='lower',extent=(-4.02,4.02,-4.02,4.02),cmap=CMAP,vmin=0,vmax=.2,interpolation='nearest',aspect='equal')
ax.set(xlim=(-4,4),ylim=(-4,4),xlabel=r'$\omega_1$',ylabel=r'$\omega_2$',title=r'$|\bar\chi^{(2)}(\omega_1,\omega_2)|$ | Fig. 2(e) comparison')
fig.colorbar(im,ax=ax,label=r'$|\bar\chi^{(2)}|$')
fig.savefig(OUT/'fig2e_second_order_map.png');fig.savefig(OUT/'fig2e_second_order_map.svg')
plt.close(fig)
# Both components retained so imaginary response can also be inspected.
fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=True,layout='constrained')
for ax,component,label in zip(axs,[np.real,np.imag],['Real','Imaginary']):
 ax.plot(wd,component(measured),color='crimson',label='Time evolution')
 ax.plot(wf,component(ref1),'--',color='darkorange',label='Perturbation (sum broadening Gamma)')
 ax.set_ylabel(label+r' $\bar\chi^{(2)}$');ax.grid(alpha=.2)
axs[0].legend(fontsize=9);axs[1].set_xlabel(r'$\omega_1=\omega_2$');axs[1].set_xlim(0,3)
fig.savefig(OUT/'diagonal_real_imag.png');plt.close(fig)
rel=abs(z-refmap)/np.maximum(abs(refmap),1e-12)
fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
a=axs[0].imshow(abs(refmap).T,origin='lower',extent=(-4.02,4.02,-4.02,4.02),vmin=0,vmax=.2,cmap=CMAP)
axs[0].set_title('Perturbation reference');fig.colorbar(a,ax=axs[0])
b=axs[1].imshow(np.log10(np.maximum(rel,1e-6)).T,origin='lower',extent=(-4.02,4.02,-4.02,4.02),vmin=-4,vmax=1,cmap='viridis')
axs[1].set_title('log10(relative difference)');fig.colorbar(b,ax=axs[1])
for ax in axs:ax.set(xlim=(-4,4),ylim=(-4,4),xlabel=r'$\omega_1$',ylabel=r'$\omega_2$')
fig.savefig(OUT/'reference_comparison.png');plt.close(fig)
# Parameter checks at omega1=0.5, normalize each by its own F and epsilon.
t=TMIN+DT*np.arange(N);testw2=np.array([.5,1.,1.5])
def check_values(path,F,eps):
 v=np.fromfile(path,dtype=np.float64);q=v[:N]-1j*v[N:];q[[0,-1]]*=.5
 ft=np.array([DT*np.sum(q*np.exp(1j*(.5+b)*t)) for b in testw2])
 fhat=F*TAUF/np.sqrt(2*np.pi)*np.exp(-testw2**2*TAUF**2/2)
 return 2*np.pi*ft/(eps*F*fhat)
original=transform(.5,testw2)
checks={}
for name,amp,eps in [('half_amplitude',.005,.001),('half_epsilon',.01,.0005)]:
 path=OUT/f'check_{name}.bin'
 if path.exists():
  checked=check_values(path,amp,eps)
  checks[name]={'omega2':testw2.tolist(),'relative_change':(abs(checked-original)/np.maximum(abs(original),1e-12)).tolist()}
mask=abs(refmap)>.01
stats={'map_points':201,'actual_positive_rows':len(ws),'negative_omega1':'reconstructed using chi(-w1,-w2)=conj(chi(w1,w2))','diagonal_complex_relative_L2_vs_S55':float(np.linalg.norm(measured-ref_diag)/np.linalg.norm(ref_diag)),'map_complex_relative_L2_vs_S55':float(np.linalg.norm(z-refmap)/np.linalg.norm(refmap)),'map_median_relative_difference_abs_reference_gt_0.01':float(np.median(rel[mask])),'permutation_relative_L2':float(np.linalg.norm(z-z.T)/np.linalg.norm(z)),'parameter_checks':checks,'analysis_seconds':time.perf_counter()-start}
(OUT/'validation.json').write_text(json.dumps(stats,indent=2)+'\n')
print(json.dumps(stats,indent=2))
