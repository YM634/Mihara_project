"""Watch the 22.5-degree scan, then save spectra and compare against x polarization."""
from pathlib import Path
import subprocess,sys,os,json
import numpy as np
root=Path(__file__).resolve().parents[1];out=root/'results/scan_angle22p5_float64';scripts=root/'scripts'
subprocess.run([sys.executable,str(scripts/'analyze_scan.py'),str(out),'--watch'],check=True)
subprocess.run([sys.executable,str(scripts/'plot_scan_window.py'),'--directory',str(out),'--gamma','0.05'],check=True)
x=np.load(root/'results/scan_x_float64/chi2_map_gamma0p05.npz');b=np.load(out/'chi2_map_gamma0p05.npz')
assert np.array_equal(x['omega1'],b['omega1']) and np.array_equal(x['omega2'],b['omega2'])
z0=x['chi'];z1=b['chi'];w=x['omega1'];diff=z1-z0
stats=dict(relative_L2=float(np.linalg.norm(diff)/np.linalg.norm(z0)),max_absolute_difference=float(np.max(abs(diff))),max_abs_chi_x=float(np.max(abs(z0))),max_abs_chi_angle22p5=float(np.max(abs(z1))),complete_cases=len(list(out.glob('*/complete.toml'))))
for q in out.glob('*/complete.toml'):
 assert (q.parent/'observables.f64').stat().st_size==140001*10*8
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/codex_mpl_op')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(1,3,figsize=(16,5),layout='constrained');extent=(-5.025,5.025,-5.025,5.025)
maxv=max(abs(z0).max(),abs(z1).max())
for ax,z,title in [(axs[0],abs(z0),r'$0^\circ$'),(axs[1],abs(z1),r'$22.5^\circ$'),(axs[2],abs(diff),r'$|\chi_{22.5}^{(2)}-\chi_0^{(2)}|$')]:
 opts=dict(vmin=0,vmax=maxv) if ax is not axs[2] else dict(vmin=0)
 im=ax.imshow(z.T,origin='lower',extent=extent,aspect='equal',cmap='magma',**opts);ax.set(title=title,xlabel=r'$\omega_1$',ylabel=r'$\omega_2$',xlim=(-5,5),ylim=(-5,5));fig.colorbar(im,ax=ax,shrink=.8)
fig.suptitle(r'Polarization comparison: $\chi^{(2)}=\bar\chi^{(2)}/2$, window $\Gamma=0.05$')
fig.savefig(out/'polarization_comparison.png',dpi=180);plt.close(fig)
(out/'comparison_x.json').write_text(json.dumps(stats,indent=2)+'\n')
size=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
with (out/'結果.md').open('a') as f:
 f.write(f'\n## Gamma=0.05でのx方向との比較\n\n複素応答の相対L2差: {stats["relative_L2"]:.6g}。最大絶対差: {stats["max_absolute_difference"]:.6g}。\n全{stats["complete_cases"]}条件のファイルサイズを確認。保存容量は約{size/1e9:.3f} GB。Gammaは後処理の窓であり、時間発展の緩和ではない。\n')
print('PIPELINE COMPLETE',json.dumps(stats),flush=True)
