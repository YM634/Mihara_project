from pathlib import Path
try:import tomllib
except ModuleNotFoundError:import tomli as tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/thread_scaling'
s=tomllib.loads((OUT/'comparison.toml').read_text());rows=s['cases']
fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
for ax,precision in zip(axes,['Float64_L50','Float64_L200','Double64_L50']):
 for update,label,marker in [(True,'Self-consistent MF','o'),(False,'Fixed MF','s')]:
  r=sorted([v for v in rows if v['precision']==precision and v['update_mean_field']==update],key=lambda v:v['threads'])
  x=[v['threads'] for v in r];y=[v['median_seconds'] for v in r]
  lo=[v['median_seconds']-min(v['evolution_seconds']) for v in r];hi=[max(v['evolution_seconds'])-v['median_seconds'] for v in r]
  ax.errorbar(x,y,yerr=[lo,hi],fmt=marker+'-',capsize=3,label=label)
 ax.set_title(precision);ax.set_xticks([1,2,4,8]);ax.set_xlabel('Julia threads');ax.set_ylabel('Evolution + observables (s)');ax.grid(alpha=.25)
axes[0].legend(fontsize=9)
fig.suptitle('k-point threading | dt=0.01 | 1,400 steps | median of 3 warmed runs')
fig.savefig(OUT/'thread_scaling_time.png',dpi=170);fig.savefig(OUT/'thread_scaling_time.svg');plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
for ax,precision in zip(axes,['Float64_L50','Float64_L200','Double64_L50']):
 for update,label,marker in [(True,'Self-consistent MF','o'),(False,'Fixed MF','s')]:
  r=sorted([v for v in rows if v['precision']==precision and v['update_mean_field']==update],key=lambda v:v['threads'])
  ax.plot([v['threads'] for v in r],[v['speedup'] for v in r],marker+'-',label=label)
 ax.axhline(1,color='.6',ls='--');ax.set_title(precision);ax.set_xticks([1,2,4,8]);ax.set_xlabel('Julia threads');ax.set_ylabel('Speedup over 1 thread');ax.grid(alpha=.25)
axes[0].legend(fontsize=9)
fig.savefig(OUT/'thread_scaling_speedup.png',dpi=170);plt.close(fig)
text='''# Juliaスレッド並列化の検証（2026-10-23）

## 実装

RK4の各段階で、異なるk点の密度行列更新をJulia Threadsで並列化した。各段階の完了を待ってから次段階のOPを計算する。各スレッドは重複しないk点だけを書き換える。

OPの総和、観測量の総和、初期平衡計算は従来の逐次処理とし、足し合わせの順番を保持した。これはスレッド数による丸め誤差の差を抑えるため。平衡・観測量の総和は並列化していないので、それらの費用は残る。

1スレッドでは従来の逐次RK4経路を使用する。複数スレッドでは各段階の微分と中間状態の更新を融合して実行する。TOMLのparallel_k=falseで複数スレッドのJuliaを起動しても逐次経路に戻せる。既定はparallel_k=true。

## 計測条件

外場はx軸から22度、f+epsilon g_cos。U=−4、J=1、T=0.02、mu=0、dt=0.01、t=−7..7、Nt=1400、F0=0.3、epsilon=0.01、omega=1、tau_f=0.1、tau_g=1。平均場更新あり・固定の両方を測定。

Float64のL=50,200、Double64のL=50。各スレッド数を独立したJuliaプロセスで、1→2→4→8の順に実行した。プロセス同士は並列実行していない。型と平均場モードごとに小規模なウォームアップを行い、各条件3回の中央値を採用。対象は時間発展と全時刻のOP・密度・電流・エネルギー評価。JIT、外場テーブル準備、平衡計算、保存は除外。

| 精度・格子 | 平均場 | 1スレッド [秒] | 2スレッド [秒] | 4スレッド [秒] | 8スレッド [秒] |
|---|---|---:|---:|---:|---:|
'''
for precision in ['Float64_L50','Float64_L200','Double64_L50']:
 for update in [True,False]:
  r=sorted([v for v in rows if v['precision']==precision and v['update_mean_field']==update],key=lambda v:v['threads'])
  text+=f"| {precision} | {'更新' if update else '固定'} | "+' | '.join(f"{v['median_seconds']:.5f} ({v['speedup']:.2f}倍)" for v in r)+' |\n'
exact=all(r['bitwise_equal'] for r in rows)
text+=f"\n環境: {s['cpu']}、Juliaが認識したCPU数{s['hardware_threads']}。全体のビット単位一致: {exact}。\n"
text+='''
## 整合性

スレッド数1の結果を基準として、全時刻の全10観測列と最終時刻の全kの3成分を比較した。Double64はFloat64に変換せず低位成分を含めて比較する。1スレッドの物理解は以前C++と検証した経路であり、今回の並列経路も独立な2×2密度行列のRK4と検証した。

| 精度・格子 | 平均場 | 全スレッド数の最大OP差 | 最大電流差 | 最大全エネルギー差 | 最終状態の最大差 |
|---|---|---:|---:|---:|---:|
'''
for precision in ['Float64_L50','Float64_L200','Double64_L50']:
 for update in [True,False]:
  r=[v for v in rows if v['precision']==precision and v['update_mean_field']==update]
  errors=[max(float(v[k]) for v in r) for k in ['max_OP_error','max_current_error','max_energy_error','max_state_error']]
  text+=f"| {precision} | {'更新' if update else '固定'} | "+' | '.join(f'{v:.3e}' for v in errors)+' |\n'
text+='''
## 速度の読み方

スレッド起動・待合せが各RK段階で必要なため、小さいFloat64計算では速度が落ちる場合がある。k点数や1点あたりの演算が大きいDouble64では、並列化の効果が出やすい。ただしメモリ帯域、逐次のOP・観測量の総和、CPUコアの性能差、温度・実行順序にも依存する。測定は順番を入れ替えたランダム試験ではないので、小さな速度差を一般化しない。

今回の表で最小時間となる設定を選ぶ。8スレッドが常に最速ではない。周波数ごとの独立ジョブ並列とは別の方式で、今回は単一時系列のk点並列のみ。2D走査は行っていない。

## 実行例

```bash
JULIA_NUM_THREADS=4 sh Attractive_Hubbard/scripts/run.sh Attractive_Hubbard/configs/fixed_mean_field_double64.toml Attractive_Hubbard/results/my_threaded_dd
```

逐次に戻す場合はJULIA_NUM_THREADS=1、または設定ファイルにparallel_k=falseを指定する。run.jlのtiming.tomlに起動スレッド数、parallel_k、active_k_threadsを保存する。active_k_threadsはk更新に利用可能な数であり、観測量や平衡計算が並列化される意味ではない。

生の3回計測はthreadsN/timing.toml、全精度を保持した時系列・状態は.jls、比較結果はcomparison.toml。再集計はcompare_threads.jl、作図はplot_thread_scaling.py。
'''
(OUT/'検証結果.md').write_text(text)
print(text[:3400])
