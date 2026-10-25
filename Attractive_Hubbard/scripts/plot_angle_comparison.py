"""Compare C++ and Julia for x-axis +22 degree polarization, in both MF modes."""
from pathlib import Path
try:import tomllib
except ModuleNotFoundError:import tomli as tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/angle22_comparison';OUT.mkdir(exist_ok=True)
modes=[('self_consistent','Self-consistent MF'),('fixed','Fixed equilibrium MF')]
sets={}
for key,label in modes:
 d=ROOT/'results'/f'angle22_{key}'
 s=tomllib.loads((d/'summary.toml').read_text())
 assert s['threads']==1 and s['fixed_mean_field']==(key=='fixed') and len(s['cases'])==3
 sets[key]=(d,s)
 for case in s['cases']:
  c=case['config'];p=np.array(c['polarization'])
  assert np.isclose(np.degrees(np.arctan2(p[1],p[0])),22,atol=1e-12)
  assert abs(np.linalg.norm(p)-1)<1e-14
  f=d/case['name']/'cpp_1';ax=np.loadtxt(f/'Ax_t.txt');ay=np.loadtxt(f/'Ay_t.txt')
  mask=abs(ax)>1e-15
  assert np.max(abs(ay[mask]/ax[mask]-np.tan(np.radians(22))))<1e-14
  assert case['same_initial_max_Phi_complex_error']<1e-12
  assert max(case['same_initial_max_errors'][1:])<1e-11
plt.rcParams.update({'savefig.dpi':170,'font.size':10})
fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True,layout='constrained')
for col,(mode,label) in enumerate(modes):
 d,_=sets[mode]
 for row,name in enumerate(['gaussian','cos_plus','sin_plus']):
  a=np.loadtxt(d/name/'julia.tsv');b=np.loadtxt(d/name/'cpp_observables.tsv')
  axes[row,col].plot(a[:,0],4*(b[:,3]-b[0,3]),color='k',lw=1.5,label='C++')
  axes[row,col].plot(a[:,0],4*(a[:,3]-a[0,3]),color='tab:orange',ls='--',lw=1.2,label='Julia')
  axes[row,col].set_title(f'{label} | {name}');axes[row,col].set_ylabel(r'Re $[\Delta(t)-\Delta(0)]$');axes[row,col].grid(alpha=.25)
axes[0,0].legend();axes[0,1].legend()
for ax in axes[-1]:ax.set_xlabel('Physical time')
fig.suptitle('Polarization: +22 degrees from x | L=50 | dt=0.01 | Float64 | single thread')
fig.savefig(OUT/'angle22_OP_comparison.png');fig.savefig(OUT/'angle22_OP_comparison.svg');plt.close(fig)
fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True,layout='constrained')
for col,(mode,label) in enumerate(modes):
 d,_=sets[mode]
 for row,name in enumerate(['gaussian','cos_plus','sin_plus']):
  a=np.loadtxt(d/name/'julia.tsv');b=np.loadtxt(d/name/'cpp_observables.tsv');cc=np.loadtxt(d/name/'julia_cpp_initial.tsv')
  phi=lambda v:v[:,3]+1j*v[:,4]
  axes[row,col].semilogy(a[:,0],np.maximum(abs(phi(a)-phi(b)),1e-18),label='Independent equilibrium')
  axes[row,col].semilogy(a[:,0],np.maximum(abs(phi(cc)-phi(b)),1e-18),label='Shared C++ initial state')
  axes[row,col].set_ylim(1e-17,1e-11);axes[row,col].set_title(f'{label} | {name}');axes[row,col].set_ylabel(r'$|\Phi_{Julia}-\Phi_{C++}|$');axes[row,col].grid(alpha=.25)
axes[0,0].legend(fontsize=8);axes[0,1].legend(fontsize=8)
for ax in axes[-1]:ax.set_xlabel('Physical time')
fig.suptitle('OP consistency at +22 degrees')
fig.savefig(OUT/'angle22_OP_errors.png');plt.close(fig)
# Current in both axes: a direct check of the oblique Peierls coupling.
fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,layout='constrained')
for col,(mode,label) in enumerate(modes):
 d,_=sets[mode];a=np.loadtxt(d/'cos_plus/julia.tsv');b=np.loadtxt(d/'cos_plus/cpp_observables.tsv')
 for row,index,yl in [(0,5,'Jx'),(1,6,'Jy')]:
  axes[row,col].plot(a[:,0],b[:,index],color='k',label='C++')
  axes[row,col].plot(a[:,0],a[:,index],color='tab:orange',ls='--',label='Julia')
  axes[row,col].set_title(label);axes[row,col].set_ylabel(yl);axes[row,col].grid(alpha=.25)
axes[0,0].legend();axes[0,1].legend()
for ax in axes[-1]:ax.set_xlabel('Physical time')
fig.suptitle('Oblique currents: f + epsilon g_cos | +22 degrees')
fig.savefig(OUT/'angle22_currents.png');plt.close(fig)
text='''# 外場をx軸から22度にしたC++/Julia比較

## 条件

偏光は単位ベクトル(cos22°, sin22°)=(0.9271838545667874, 0.374606593415912)。入力ファイルから実際の角度、単位長、およびAy/Ax=tan22°を確認した（極小の裾を除く）。

両方Float64、1プロセス・1スレッド、L=50、J=1、U=−4、T=0.02、mu=0、S波、g_elph=0、dt=0.01、t=−7..7、Nt=1400。平衡の収束閾値1e−12。F0=0.3、epsilon=0.01、omega=1、tau_f=0.1、tau_g=1。

fのみ、f+epsilon g_cos、f+epsilon g_sinの3外場を、平均場更新あり・なしでそれぞれ比較した。両コードの外場配列と5点4次補間を揃える。極小値のC++読込停止を防ぐ共通の1e−300カットオフも前回と同じ。

C++は元ソースから作成済みの自己無撞着版と固定版の参照実行ファイルを使用。Juliaではupdate_mean_field=true/falseを切り替えた。独立平衡とC++初期状態を共有する2通りを比較し、位相合わせやOP・電流の補正は行っていない。

## 最大絶対差

Phiは対期待値、Delta=−U Phi。全時刻で比較。

| 平均場 | 外場 | 独立平衡の複素Phi差 | 同じ初期状態の複素Phi差 | 同じ初期状態の電流差 | 同じ初期状態の全エネルギー差 |
|---|---|---:|---:|---:|---:|
'''
for mode,label in modes:
 d,s=sets[mode]
 for r in s['cases']:
  e=r['same_initial_max_errors'];text+=f"| {label} | {r['name']} | {r['max_Phi_complex_error']:.3e} | {r['same_initial_max_Phi_complex_error']:.3e} | {max(e[5:7]):.3e} | {e[9]:.3e} |\n"
text+='''
平均場を固定する場合のC++のEint/Etotだけ、既知の定数−U|Phi_eq|²を加え、Juliaのエネルギー原点と揃えている。定数は運動方程式には影響しない。生のC++エネルギーも各条件のcpp_observables_raw.tsvに保存する。

## 判定

今回の22度偏光でも、平均場更新あり・固定の両方で、C++とJuliaが倍精度の丸め誤差程度まで一致した。初期状態を独立に求めた場合の差には、収束残差の違いが含まれる。これは今回の時間・格子条件での実装間の一致であり、全パラメータでの誤差や格子・時間刻みの収束を保証するものではない。2D走査は実施していない。

## 再実行と保存先

```bash
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/compare_cpp.jl --angle=22
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/compare_cpp.jl --angle=22 --fixed
MPLCONFIGDIR=/private/tmp/codex_mpl_op python3 Attractive_Hubbard/scripts/plot_angle_comparison.py
```

生データと誤差集計はresults/angle22_self_consistent、results/angle22_fixed。各summary.tomlには条件・誤差と各1回の参考実行時間を保存する。

このフォルダのangle22_OP_comparison.pngはOP波形、angle22_OP_errors.pngは複素OPの差、angle22_currents.pngはx,y方向の電流の比較。
'''
(OUT/'比較結果.md').write_text(text)
print(text)
