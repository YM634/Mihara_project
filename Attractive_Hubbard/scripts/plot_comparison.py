from pathlib import Path
try: import tomllib
except ModuleNotFoundError: import tomli as tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/cpp_comparison'
s=tomllib.loads((OUT/'summary.toml').read_text())
rows=s['cases'];names={r['name'] for r in rows}
required={'zero','gaussian_x','ono_diagonal','doped_sin','L20_gaussian'}
if names!=required:raise RuntimeError('Benchmark suite is incomplete')
plt.rcParams.update({'font.size':10,'savefig.dpi':170})
fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True,layout='constrained')
for i,(name,label) in enumerate([('gaussian_x','Gaussian / x polarization'),('ono_diagonal','f + epsilon g_cos / diagonal polarization'),('doped_sin','f + epsilon g_sin / mu=0.5')]):
 d=OUT/name;a=np.loadtxt(d/'julia.tsv');b=np.loadtxt(d/'cpp_observables.tsv');cc=np.loadtxt(d/'julia_cpp_initial.tsv')
 ja=4*(a[:,3]+1j*a[:,4]);cb=4*(b[:,3]+1j*b[:,4]);cj=4*(cc[:,3]+1j*cc[:,4]);t=a[:,0]
 axes[i,0].plot(t,(cb-cb[0]).real,label='C++',color='k',lw=1.5)
 axes[i,0].plot(t,(ja-ja[0]).real,label='Julia',color='tab:orange',lw=1.1,ls='--')
 axes[i,0].set_title(label);axes[i,0].set_ylabel(r'Re $[\Delta(t)-\Delta(0)]$')
 axes[i,1].semilogy(t,np.maximum(abs(ja-cb),1e-18),label='Independent equilibrium',lw=1)
 axes[i,1].semilogy(t,np.maximum(abs(cj-cb),1e-18),label='Same C++ initial state',lw=1)
 axes[i,1].set_title('Absolute complex gap difference');axes[i,1].set_ylabel(r'$|\Delta_{Julia}-\Delta_{C++}|$')
 axes[i,1].set_ylim(1e-17,1e-10)
 for ax in axes[i]:ax.grid(alpha=.25)
axes[0,0].legend();axes[0,1].legend(fontsize=8)
for ax in axes[-1]:ax.set_xlabel('Physical time')
fig.suptitle('2D Attractive Hubbard | L=50 | dt=0.01 | Float64 | no parallelism')
fig.savefig(OUT/'comparison_dynamics.png');fig.savefig(OUT/'comparison_dynamics.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
x=np.arange(len(rows));cmed=[];jmed=[];clo=[];chi=[];jlo=[];jhi=[]
for r in rows:
 c=np.array([v['evolution_seconds'] for v in r['cpp_runs']]);j=np.array([v['evolution_seconds'] for v in r['julia_runs']]);cm=np.median(c);jm=np.median(j)
 cmed.append(cm);jmed.append(jm);clo.append(cm-min(c));chi.append(max(c)-cm);jlo.append(jm-min(j));jhi.append(max(j)-jm)
ax.bar(x-.18,cmed,.36,label='C++ reference',yerr=[clo,chi],capsize=3)
ax.bar(x+.18,jmed,.36,label='Julia specialized',yerr=[jlo,jhi],capsize=3)
ax.set_yscale('log');ax.set_xticks(x,[r['name'] for r in rows]);ax.set_ylabel('Evolution + shared observables (seconds)')
ax.set_title('Single-thread comparison | 1,400 steps | median of 3 runs | JIT / I/O excluded')
for i,r in enumerate(rows):ax.text(i,cmed[i]*1.25,f"{r['speedup_dynamics']:.1f}x",ha='center')
ax.set_ylim(min(jmed)/2,max(cmed)*2);ax.legend();ax.grid(axis='y',alpha=.25)
fig.savefig(OUT/'comparison_speed.png');fig.savefig(OUT/'comparison_speed.svg');plt.close(fig)
# Plot current and total energy to check more than the order parameter.
d=OUT/'ono_diagonal';a=np.loadtxt(d/'julia.tsv');b=np.loadtxt(d/'cpp_observables.tsv')
fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,layout='constrained')
for ax,col,label in zip(axes.flat,[4,5,6,9],['Im Phi','Jx','Jy','Total energy']):
 ax.plot(a[:,0],b[:,col],label='C++',lw=1.5,color='k');ax.plot(a[:,0],a[:,col],label='Julia',lw=1,color='tab:orange',ls='--');ax.set_ylabel(label);ax.grid(alpha=.25)
axes[0,0].legend()
for ax in axes[-1]:ax.set_xlabel('Physical time')
fig.savefig(OUT/'comparison_observables.png');plt.close(fig)

header='''# C++とJuliaの結果・速度比較

2D Attractive Hubbardの専用Juliaコードを作成し、既存C++と比較した。2D走査・スペクトル計算は行っていない。

## 比較条件

両方とも通常の倍精度、1プロセス・1スレッド。dt=0.01、t=−7..7、Nt=1400（1401点）、温度0.02、U=−4、J=1、S波、g_elph=0。平衡の平均密度行列変化Frobeniusノルムの閾値1e−12。mu=0が基本で、doped_sinのみmu=0.5。

外場はF0=0.3（AのGaussian係数はF0/(2pi)）、epsilon=0.01、omega=1、tau_f=0.1、tau_g=1。zeroは外場ゼロ、gaussian_xはx方向fのみ、ono_diagonalは単位対角方向のf+epsilon g_cos、doped_sinは偏光(1,0.5)を単位化したf+epsilon g_sin。L20_gaussianはL=20、それ以外L=50。

外場テーブルを両実装で共有し、RK4の中間時刻は5点4次補間。Juliaでは補間と三角関数のうち全k共通の部分を一度だけ評価する。両方とも全時刻でOP・密度・電流・エネルギーを評価する。

## 一致度

以下は全時刻の複素Phiの最大絶対差。Delta=−U Phiなので、ギャップの差はこの4倍。

| 条件 | 最大Phi差 | 最大Delta差 | C++のPhi変化の最大値 |
|---|---:|---:|---:|
'''
body=[]
for r in rows:body.append(f"| {r['name']} | {r['max_Phi_complex_error']:.3e} | {r['max_Delta_complex_error']:.3e} | {r['cpp_max_Phi_change']:.3e} |")
text=header+'\n'.join(body)+'\n\n'
text+='''位相を後からフィット・回転していない。独立な自己無撞着平衡解をそのまま比較した。Juliaの反復初期値はPhi=0.35、C++は乱数であり、同じ閾値でも微小な残差は異なる。C++の初期密度行列をJuliaに読み込む比較も追加した。

| 条件 | 同じ初期状態でのPhi差（実部・虚部ごとの最大） | 最大電流差 | 最大全エネルギー差 |
|---|---:|---:|---:|
'''
for r in rows:
 e=r['max_errors'];cc=r['same_initial_max_errors'];text+=f"| {r['name']} | {max(cc[3:5]):.3e} | {max(e[5:7]):.3e} | {e[9]:.3e} |\n"
text+='\n## 速度\n\n各条件3回の中央値。計測区間はRK4時間発展と共通の観測量評価。平衡計算、Julia JIT、入力生成、ファイル保存を除外した。最小値・最大値を図のエラーバーで表示。\n\n| 条件 | C++ [秒] | Julia [秒] | 速度比 C++/Julia |\n|---|---:|---:|---:|\n'
for r,cm,jm in zip(rows,cmed,jmed):text+=f"| {r['name']} | {cm:.5f} | {jm:.5f} | {cm/jm:.1f}倍 |\n"
text+=f"\n環境: {s['cpu']}、Julia {s['julia_version']}、Julia Threads.nthreads()={s['threads']}。初回の小規模ウォームアップは{s['warmup_seconds']:.6f}秒。C++のコンパイラ・コマンド・参照ソースのSHA256はreference/build.json。\n"
text+='''
これは言語だけを入れ替えた速度差ではない。Julia版はS波・フォノンなしに専用化し、密度行列を3実成分で保持し、作業領域を再利用する。C++は汎用の動的2×2行列・フォノン計算を経由し、外場補間をkごとに繰り返す。主にそのアルゴリズム・メモリ配置・共通計算の削減を含む高速化である。

## 元の実行ファイルとの照合

'''
g=next(r for r in rows if r['name']=='gaussian_x')
text+=f"元のExecute_sc/sc_model.ex（12桁出力）ともgaussian_xで照合し、17桁出力にした参照版との差はPhi成分で最大{max(g['original_vs_instrumented_max_errors'][3:5]):.3e}。出力の丸めを含む。\n\n"
text+=f"Gaussianの中間時刻を解析式で評価するJuliaモードと、C++に揃えたtable4モードとのPhiの最大差は{g['analytic_vs_table4_max_Phi_error']:.3e}。比較の主結果はtable4を使用。\n"
text+='''
## 発見したC++の入力読込問題

純Gaussianの裾にある非正規化数（例5e−324）で、この環境のC++のistreamがfailbitを立てた。元コードは読込失敗を確認せず、以後の値がゼロのままになる。元のGaussian入力ではAxの315値目までしか読めず、パルス本体に到達していなかった。これを物理・数値解法の不一致と混同しないよう切り分けた。

今回の比較入力は、abs(value)<1e−300の値を両実装共通でゼロにした。パルス本体は一切変えず、C++の運動方程式や読込コードも変更していない。診断データはresults/input_reader_diagnostic/、単独確認プログラムはreference/check_field_reader.cpp。古いrun_job_Onoのコード・元実行ファイルは変更していない。

## 精度と今回の範囲

独立な2×2複素密度行列によるRK4、4次補間、ゼロ場、既知のギャップ等の16テストと、Double64/106-bit BigFloat比較の2テストに合格した。

この一致は実装間の照合であり、dt・格子・観測時間の収束を保証するものではない。今回の短時間ベンチマークから、長時間や微小二次応答の精度を断定しない。2D scanは実行していない。

## ファイル

- comparison_dynamics.png / .svg: OP変化と実装間の差。
- comparison_speed.png / .svg: 時間発展の速度比較。
- comparison_observables.png: 電流・虚部OP・エネルギー。
- summary.toml: 各反復の生の計測値と誤差。
- 各条件のjulia.tsv、cpp_observables.tsv、julia_cpp_initial.tsv: 全時系列。
- cpp_1..3: 共通の外場入力、C++出力、初期状態、区間別の計測値。
'''
(OUT/'比較結果.md').write_text(text)
print(text[:2500])
