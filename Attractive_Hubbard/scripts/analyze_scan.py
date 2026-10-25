"""Extract symmetrized second-order Delta response to Ax; reusable raw traces remain intact."""
import argparse, json, time, os
from pathlib import Path
import numpy as np
from scipy.signal import ZoomFFT
try: import tomllib
except ModuleNotFoundError: import tomli as tomllib
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/codex_mpl_op')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

parser=argparse.ArgumentParser()
parser.add_argument('directory',type=Path)
parser.add_argument('--watch',action='store_true')
args=parser.parse_args();out=args.directory.resolve()
raw=tomllib.loads((out/'scan_config.toml').read_text());c=raw['model'];scan=raw['scan']
pol=np.asarray(c['polarization'],dtype=float);pol/=np.linalg.norm(pol);angle=float(np.degrees(np.arctan2(pol[1],pol[0])))
dt=float(c['dt']);tmin=float(c['tmin']);n=round((float(c['tmax'])-tmin)/dt)+1
t=tmin+dt*np.arange(n);wpos=np.linspace(scan['omega_min'],scan['omega_max'],scan['omega_count']);w=np.r_[-wpos[:0:-1],wpos]
lo,hi=scan['amplitudes'];fl,fh=float(lo),float(hi);ratio=fh/fl;eps=float(c['epsilon']);tauf=float(c['tau_f']);U=float(c['U'])
fhhat=fh*tauf/np.sqrt(2*np.pi)*np.exp(-w*w*tauf*tauf/2)
eta=.01
z=np.full((len(w),len(w)),np.nan+1j*np.nan);zd=z.copy();done=set();checks=[]

def read(i,amp,kind):
    p=out/f'omega_{i:03d}_F{amp}_{kind}'
    meta=tomllib.loads((p/'complete.toml').read_text())
    if meta['rows']!=n or meta['columns']!=10: raise ValueError(str(p))
    a=np.memmap(p/'observables.f64',dtype='<f8',mode='r',shape=(n,10),order='F')
    if not np.all(np.isfinite(a)):raise ValueError('Nonfinite '+str(p))
    return a

def central(i,amp,kind):
    p=read(i,amp,kind+'_plus');m=read(i,amp,kind+'_minus')
    # Observable is real Delta=-U Re(Phi), so real-signal conjugation applies.
    return -U*(p[:,3]-m[:,3])/2

def ft(q,om):
    q=q.copy();q[[0,-1]]*=.5
    v=ZoomFFT(n,[-om[0]/(2*np.pi),-om[-1]/(2*np.pi)],m=len(om),fs=1/dt,endpoint=True)(q)
    return dt*np.exp(1j*om*tmin)*v

while len(done)<len(wpos):
    changed=False
    for i,a in enumerate(wpos):
        if i in done:continue
        if not all((out/f'omega_{i:03d}_F{amp}_{kind}_{sign}'/'complete.toml').exists() for amp in (lo,hi) for kind in ('cos','sin') for sign in ('plus','minus')):continue
        lc,hc=central(i,lo,'cos'),central(i,hi,'cos')
        ls,hs=central(i,lo,'sin'),central(i,hi,'sin')
        qc=(hc-ratio*lc)/(1-1/ratio);qs=(hs-ratio*ls)/(1-1/ratio)
        q=qc-1j*qs;om=a+w;factor=2*np.pi/(eps*fh*fhhat)
        vals=ft(q,om)*factor
        damped=ft(q*np.exp(-eta*np.maximum(t,0)),om)*factor
        j=len(wpos)-1+i;z[j]=vals;zd[j]=damped
        if i>0:z[len(wpos)-1-i]=vals[::-1].conj();zd[len(wpos)-1-i]=damped[::-1].conj()
        # Fourier implementation checked independently at selected map columns.
        selected=[0,len(w)//2,len(w)-1]
        qt=q.copy();qt[[0,-1]]*=.5
        direct=np.array([dt*np.sum(qt*np.exp(1j*om[k]*t))*factor[k] for k in selected])
        directerr=float(np.max(abs(direct-vals[selected]))/max(np.max(abs(direct)),1e-30))
        high=hc-1j*hs;low=lc-1j*ls
        scaleerr=float(np.linalg.norm(high-ratio**2*low)/max(np.linalg.norm(high),1e-30))
        checks.append(dict(omega1=float(a),fourier_relative_error=directerr,amplitude_scaling_relative_L2=scaleerr))
        done.add(i);changed=True
        print(f'{time.strftime("%H:%M:%S")} transformed {len(done)}/{len(wpos)} omega1={a:.2f} scaling_error={scaleerr:.3g}',flush=True)
    if changed:
        np.savez_compressed(out/'spectra_progress.npz',omega1=w,omega2=w,chi_bar=z,chi_bar_eta001=zd,completed_positive_indices=sorted(done))
    if len(done)==len(wpos):break
    if not args.watch:raise SystemExit(f'Only {len(done)}/{len(wpos)} rows complete; use --watch')
    time.sleep(5)

# Common pulse-only cases must also complete before declaring the scan finished.
while not all((out/f'common_F{amp}_f'/'complete.toml').exists() for amp in (lo,hi)):
    if not args.watch:raise SystemExit('Common cases missing')
    time.sleep(5)
np.savez_compressed(out/'spectra.npz',omega1=w,omega2=w,chi_bar=z,chi_bar_eta001=zd,chi=z/2,chi_eta001=zd/2,eta=eta)
plt.rcParams.update({'font.size':11,'savefig.dpi':180})
extent=(w[0]-.025,w[-1]+.025,w[0]-.025,w[-1]+.025)
for name,v,label in [('chi2_map',z,'Finite-time transform (no window)'),('chi2_map_eta001',zd,r'Postprocessing window: $e^{-0.01\max(t,0)}$')]:
    fig,axs=plt.subplots(1,3,figsize=(17,5),layout='constrained')
    for ax,data,title in zip(axs,[abs(v),v.real,v.imag],[r'$|\bar\chi^{(2)}|$',r'$\mathrm{Re}\,\bar\chi^{(2)}$',r'$\mathrm{Im}\,\bar\chi^{(2)}$']):
        if title.startswith('$|'):
            im=ax.imshow(data.T,origin='lower',extent=extent,aspect='equal',cmap='magma')
        else:
            limit=float(np.max(abs(data)));im=ax.imshow(data.T,origin='lower',extent=extent,aspect='equal',cmap='RdBu_r',vmin=-limit,vmax=limit)
        ax.set(xlabel=r'$\omega_1$',ylabel=r'$\omega_2$',title=title,xlim=(-5,5),ylim=(-5,5));fig.colorbar(im,ax=ax,shrink=.8)
    fig.suptitle(rf'Attractive Hubbard: $\Delta=-U\,\mathrm{{Re}}\,\Phi$, polarization ${angle:g}^\circ$'+'\n'+label)
    fig.savefig(out/f'{name}.png');fig.savefig(out/f'{name}.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
im=ax.imshow(abs(zd).T,origin='lower',extent=extent,aspect='equal',cmap='magma');ax.set(xlabel=r'$\omega_1$',ylabel=r'$\omega_2$',title=rf'$|\bar\chi^{{(2)}}|$ (window $\eta=0.01$, polarization ${angle:g}^\circ$)',xlim=(-5,5),ylim=(-5,5));fig.colorbar(im,ax=ax);fig.savefig(out/'chi2_abs.png');plt.close(fig)
fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True,layout='constrained')
for ax,func,label in [(axs[0],np.real,'Real'),(axs[1],np.imag,'Imaginary')]:
    ax.plot(w,func(np.diag(z)),label='No window',alpha=.7);ax.plot(w,func(np.diag(zd)),label='Window eta=0.01');ax.set_ylabel(label+r' $\bar\chi^{(2)}(\omega,\omega)$');ax.grid(alpha=.2)
axs[0].legend();axs[1].set_xlabel(r'$\omega_1=\omega_2$');fig.savefig(out/'chi2_diagonal.png');plt.close(fig)
metrics=dict(complete_cases=len(list(out.glob('*/complete.toml'))),expected_cases=810,grid_points=len(w),observable='-U Re(Phi)',input='vector potential along polarization (NOT electric E)',polarization=pol.tolist(),angle_degrees=angle,symmetrization='chi_bar(w1,w2)=chi(w1,w2)+chi(w2,w1); chi=chi_bar/2',window_eta=eta,permutation_relative_L2_raw=float(np.linalg.norm(z-z.T)/np.linalg.norm(z)),permutation_relative_L2_eta001=float(np.linalg.norm(zd-zd.T)/np.linalg.norm(zd)),max_fourier_relative_error=max(x['fourier_relative_error'] for x in checks),amplitude_scaling_relative_L2_median=float(np.median([x['amplitude_scaling_relative_L2'] for x in checks])),checks=checks)
(out/'analysis.json').write_text(json.dumps(metrics,indent=2)+'\n')
report=f'''# x軸から{angle:g}度方向の二次応答スキャンの結果

全{metrics['complete_cases']}条件完了。Float64、L=50、U=-4、T=0.02、mu=0、平均場更新あり。dt=0.01、t=-700..700。omega1とomega2は-5..5、0.05刻みの201点。非負omega1の101行を実計算し、負側は実観測量の共役対称性から復元した。

観測量はDelta=-U Re(Phi)、入力は偏光方向のベクトルポテンシャルA（偏光ベクトル={pol.tolist()}）。電場Eに対する応答とは異なる。図は対称化した和chi_barを示す。対称な展開係数chiはその1/2であり、spectra.npzに両方を保存した。

F_low=0.001、F_high=0.01、epsilon=0.001。cos/sinの正負中央差分を取り、(D_high-10 D_low)/0.9で二次推定を構成した。フーリエ変換はexp(+i omega t)、台形則、omega=omega1+omega2。有限幅のテストパルスを無限長の単色波で近似した抽出であり、完全な二周波数応答とは有限幅誤差がある。

chi2_map.pngは窓なしの有限時間結果、chi2_map_eta001.pngとchi2_abs.pngは後処理でexp[-0.01 max(t,0)]を掛けた結果。運動方程式には緩和を加えていない。この窓は見え方・有限時間感度の比較用で、物理的緩和を入れた再計算ではない。

独立直接積分との最大相対差: {metrics['max_fourier_relative_error']:.3g}。
交換対称性の相対L2差: 窓なし{metrics['permutation_relative_L2_raw']:.3g}、窓あり{metrics['permutation_relative_L2_eta001']:.3g}。
高低外場のF^2スケーリング相対L2差（周波数ごとの中央値）: {metrics['amplitude_scaling_relative_L2_median']:.3g}。
これらは診断であり、dt、格子数、時間範囲、epsilonを変えた収束試験は未実施。特に無減衰の有限時間変換では共鳴ピークの高さが時間窓に依存する。

元時系列は全て保存済み。解析を変えて再計算できる。
'''
(out/'結果.md').write_text(report)
print('ANALYSIS COMPLETE '+json.dumps({k:v for k,v in metrics.items() if k!='checks'}),flush=True)
