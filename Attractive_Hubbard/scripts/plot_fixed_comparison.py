from pathlib import Path
try:import tomllib
except ModuleNotFoundError:import tomli as tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cpp_fixed_comparison'
s=tomllib.loads((OUT/'summary.toml').read_text());rows=s['cases']
assert s['fixed_mean_field'] and s['threads']==1 and len(rows)==5
plt.rcParams.update({'savefig.dpi':170,'font.size':10})
fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True,layout='constrained')
for i,name in enumerate(['gaussian_x','ono_diagonal','doped_sin']):
 d=OUT/name;a=np.loadtxt(d/'julia.tsv');b=np.loadtxt(d/'cpp_observables.tsv');cc=np.loadtxt(d/'julia_cpp_initial.tsv')
 ja=4*(a[:,3]+1j*a[:,4]);cb=4*(b[:,3]+1j*b[:,4]);cj=4*(cc[:,3]+1j*cc[:,4]);t=a[:,0]
 axes[i,0].plot(t,(cb-cb[0]).real,label='C++ fixed MF',color='k',lw=1.5)
 axes[i,0].plot(t,(ja-ja[0]).real,label='Julia fixed MF',color='tab:orange',ls='--')
 axes[i,0].set_title(name);axes[i,0].set_ylabel(r'Re $[\Delta(t)-\Delta(0)]$')
 axes[i,1].semilogy(t,np.maximum(abs(ja-cb),1e-18),label='Independent equilibrium')
 axes[i,1].semilogy(t,np.maximum(abs(cj-cb),1e-18),label='Shared C++ initial state')
 axes[i,1].set_ylim(1e-17,1e-10);axes[i,1].set_ylabel(r'$|\Delta_{Julia}-\Delta_{C++}|$')
 axes[i,1].set_title('Absolute complex gap difference')
 for ax in axes[i]:ax.grid(alpha=.25)
axes[0,0].legend();axes[0,1].legend(fontsize=8)
for ax in axes[-1]:ax.set_xlabel('Physical time')
fig.suptitle('Fixed equilibrium mean field | L=50 | dt=0.01 | Float64 | single thread')
fig.savefig(OUT/'fixed_comparison_dynamics.png');fig.savefig(OUT/'fixed_comparison_dynamics.svg');plt.close(fig)
# Difference in all measured non-OP quantities, using the same initial state.
fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,layout='constrained')
for name in ['gaussian_x','ono_diagonal','doped_sin']:
 d=OUT/name;b=np.loadtxt(d/'cpp_observables.tsv');a=np.loadtxt(d/'julia_cpp_initial.tsv')
 for ax,col,label in zip(axes.flat,[1,5,6,9],['n_up','Jx','Jy','Etot (constant aligned)']):
  ax.plot(a[:,0],a[:,col]-b[:,col],label=name,lw=.8);ax.set_ylabel(f'Julia - C++: {label}');ax.grid(alpha=.25)
axes[0,0].legend(fontsize=8)
for ax in axes[-1]:ax.set_xlabel('Physical time')
fig.savefig(OUT/'fixed_comparison_observable_errors.png');plt.close(fig)
text='''# C++固定平均場版とJuliaの整合性確認

## C++側の実装

元のSC_phonon/sc_dynamics_main_fixMF.cppは、time_evoluation_tstp_fixMFとevaluate_observales_fixMFを呼ぶ。時間発展の各段階でrho_el_k_eq（平衡時の密度行列）を平均場構築に使い、A(t)は各時刻で更新する。これはJuliaのupdate_mean_field=falseと同じ固定方法。

元の固定平均場ソースをそのままコンパイルし、mainだけに17桁出力・初期状態出力・計測を追加した。元ファイルは変更していない。reference/build_fixed.jsonに元ソースのハッシュとビルドコマンドを保存。元のsc_model_fixMF.exはx86_64版であるため、今回の比較には現環境向けに再コンパイルしたsc_model_fixed_reference.exを使用した。

## 比較条件

双方Float64、単一スレッド。L=50（追加でL=20）、J=1、U=−4、T=0.02、S波、g_elph=0、dt=0.01、t=−7..7、Nt=1400。通常mu=0、doped_sinだけmu=0.5。平衡収束閾値1e−12。

外場F0=0.3、epsilon=0.01、omega=1、tau_f=0.1、tau_g=1。zeroはゼロ場、gaussian_xはx方向のf、ono_diagonalは単位対角方向のf+epsilon g_cos、doped_sinは単位化した偏光(1,0.5)のf+epsilon g_sin。外場配列を共通にし、5点4次補間を使用。非正規化数でのC++読込停止を避け、abs(value)<1e−300を両実装共通で0にする。

独立に平衡計算した比較と、C++から出力した同じ初期状態をJuliaに与えた比較を実施。位相合わせ等のフィットは行っていない。速度は今回1回の参考値で、前の自己無撞着版の3回中央値とは区別する。

## OPの比較

Phiは対期待値で、Delta=−U Phi。以下は全時刻での複素数の最大絶対差。

| 条件 | 独立平衡のPhi差 | 同じ初期状態のPhi差 | 同じ初期状態のDelta差 |
|---|---:|---:|---:|
'''
for r in rows:
 e=r['same_initial_max_Phi_complex_error'];text+=f"| {r['name']} | {r['max_Phi_complex_error']:.3e} | {e:.3e} | {4*e:.3e} |\n"
text+='\n## その他の観測量\n\n同じ初期状態を使った成分ごとの最大絶対差。\n\n| 条件 | 密度 | 電流 | Ekin | Eint（定数補正後） | Etot（定数補正後） |\n|---|---:|---:|---:|---:|---:|\n'
for r in rows:
 e=r['same_initial_max_errors'];text+=f"| {r['name']} | {max(e[1:3]):.3e} | {max(e[5:7]):.3e} | {e[7]:.3e} | {e[8]:.3e} | {e[9]:.3e} |\n"
text+=r'''
## エネルギーの定数差

C++固定版は

$$
E_{\mathrm{int,C++}}=2U\operatorname{Re}[\Phi_{\mathrm{eq}}^*\Phi(t)]
$$

を出力する。Julia版は平衡時の自己無撞着エネルギーに合わせた定数を含めて

$$
E_{\mathrm{int,Julia}}=2U\operatorname{Re}[\Phi_{\mathrm{eq}}^*\Phi(t)]-U|\Phi_{\mathrm{eq}}|^2
$$

とする。そのため、C++のEintとEtotに既知の定数−U|Phi_eq|²を加えて比較した。定数は時間発展には影響せず、エネルギーの時間変化には補正不要。OP・電流・密度には何も補正していない。生のC++エネルギーはcpp_observables_raw.tsv、補正後はcpp_observables.tsvに保存。

## 実行時間（参考・各1回）

| 条件 | C++時間発展 [秒] | Julia時間発展 [秒] |
|---|---:|---:|
'''
for r in rows:text+=f"| {r['name']} | {r['cpp_runs'][0]['evolution_seconds']:.5f} | {r['julia_runs'][0]['evolution_seconds']:.5f} |\n"
text+='''
## 結論と範囲

今回の条件では、固定平均場でのOP・密度・電流・エネルギーが倍精度の丸め誤差程度まで一致し、Juliaの固定モードとC++の固定平均場版の整合性を確認した。独立平衡時の差には初期反復・収束残差の違いを含む。この比較は実装の一致を確認するもので、格子・時間刻み・長時間の収束を保証しない。double-double版C++との比較や2D走査は行っていない。

## 再実行

```bash
python3 Attractive_Hubbard/scripts/build_reference.py --fixed
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/compare_cpp.jl --fixed
MPLCONFIGDIR=/private/tmp/codex_mpl_op python3 Attractive_Hubbard/scripts/plot_fixed_comparison.py
```

同一条件の完了ケースは再利用する。結果・誤差・計測はsummary.toml、OP比較図はfixed_comparison_dynamics.png、他の観測量の差はfixed_comparison_observable_errors.png。
'''
(OUT/'比較結果.md').write_text(text)
# Check all quantities after aligning the fixed energy reference.
for r in rows:
 assert r['same_initial_max_Phi_complex_error']<1e-12
 assert max(r['same_initial_max_errors'][1:])<1e-11
print(text[:2400])
