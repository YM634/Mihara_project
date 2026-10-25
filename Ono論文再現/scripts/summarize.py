from pathlib import Path
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
results=root/'results'
f=tomllib.loads((results/'benchmark_float64/timing.toml').read_text())
d=tomllib.loads((results/'benchmark_double64/timing.toml').read_text())
l=tomllib.loads((results/'benchmark_long_double64/timing.toml').read_text())
rate=l['evolution_total_seconds']/8000/6
rows=[]
for width in [1400,2000]:
 steps=round(width/.01);one=rate*steps
 rows.append(f'| {width} | {steps:,} | {one/60:.2f}分 | {6*one/60:.2f}分 | {806*one/3600:.2f}時間 |')
report=f'''# 計算時間の実測と見積もり

計測日: 2026-10-22。{d['cpu']}、Julia {d['julia_version']}、{d['threads']}スレッド。スペクトルは未計算。

## 実測

Nk=800、dt=0.01、hx=1、hy=0.5、hz=0.1、Gamma=0.01、tau_f=0.1、omega1=0.5、epsilon=1e−6。F0=1e−4,1e−3の各々でf、f+epsilon g_cos、f+epsilon g_sinを計算。1セット=6本。

| 精度 | 時間区間 | ステップ/本 | tau_g | 6本の時間発展 | 6本の時系列保存 |
|---|---|---:|---:|---:|---:|
| Float64 | −20..20 | 4,000 | 3 | {f['evolution_total_seconds']:.3f}秒 | {f['write_total_seconds']:.3f}秒 |
| Double64 | −20..20 | 4,000 | 3 | {d['evolution_total_seconds']:.3f}秒 | {d['write_total_seconds']:.3f}秒 |
| Double64 | −40..40 | 8,000 | 6 | {l['evolution_total_seconds']:.3f}秒 | {l['write_total_seconds']:.3f}秒 |

Double64の1本は短い方で{min(r['evolution_seconds'] for r in d['runs']):.2f}–{max(r['evolution_seconds'] for r in d['runs']):.2f}秒、長い方で{min(r['evolution_seconds'] for r in l['runs']):.2f}–{max(r['evolution_seconds'] for r in l['runs']):.2f}秒。区間2倍の計算時間比は{l['evolution_total_seconds']/d['evolution_total_seconds']:.2f}倍で、概ねステップ数に比例する。

型ごとの初回JITウォームアップはFloat64で{f['warmup_seconds']:.2f}秒、Double64で{d['warmup_seconds']:.2f}秒（短い計測時）。これは上の時間発展から除外した。Julia起動・パッケージロード・次数分離関数の初回コンパイル・差分ファイルの保存はこの表の合計に含まない。DoubleFloatsの初回依存導入・プリコンパイルも除外。従って表はコア計算と6時系列の保存の計測で、シェルからの全経過時間ではない。

## 論文スケールの観測時間への外挿

長い方のDouble64実測から外挿。tau_g=sqrt(2)/Gamma≈141.421に戻す。カーネルの演算回数はtau_gにほぼ依存しないが、長時間本番は未実行なので概算。

| 観測時間幅 | ステップ/本 | 1本 | 1周波数6本 | 201周波数（806本） |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}

201は将来の走査点数の例で、論文が201点と指定している意味ではない。1つのomega1についてomega2全域を同じ時系列から得るので、201×201回の時間発展は不要。fの2振幅はomega1に依存せず共有でき、ケース数は2+4M、M=201なら806本。毎回6本を再実行すると1206本になり約1.5倍かかる。

`paper.toml` の時間幅2000なら、1セットは約{rate*200000*6/60:.1f}分、201周波数は逐次約{rate*200000*806/3600:.1f}時間。独立ジョブを8個同時に動かせる同一速度の環境なら理想的には約{rate*200000*806/3600/8:.1f}時間だが、CPU速度、熱、コア競合、保存負荷で変わる。今回のコード自体は単一スレッド。dt半減は概ね2倍、Nk倍増は概ね2倍。

観測区間とtau_gはスペクトル収束で最終決定する。今回の短い区間は演算速度の実測が目的であり、論文の周波数分解能を再現する設定ではない。FT、摂動論、作図、追加の収束計算の費用は上表に含まない。

## 精度と保存量

短い計測の二次cos混合変分の最大値は、Double64で{d['max_abs_delta_cos']}、Float64で{f['max_abs_delta_cos']}。同じ条件でもFloat64は微小な差分の丸め誤差が大きい。時間が速いことを理由に論文パラメータの二次抽出をFloat64で代用すべきではない。これは未収束の短時間試験であり、Double64の値を物理的な正解と認定する比較ではない。

保存はtime,A,Jの3列。Double64は有効桁を保持した文字列で保存し、Float64へ変換しない。Nk×Ntの全密度行列は保存せず、現在のk状態と観測時系列を保持する。時間幅2000では約200,001行/本。現在の実ファイルサイズに比例させた容量目安は後述の自動算出値を参照。timing.tomlのallocated_bytesは累積割当量でピークRAMではない。
'''
# actual text size extrapolation; pulse and value digit lengths can vary
single_bytes=(results/'benchmark_long_double64/high_cos.tsv').stat().st_size
bytes_full=single_bytes*200001/8001
report+=f'\n長い計測のhigh_cos.tsvは{single_bytes:,} bytes。時間幅2000の1本は約{bytes_full/1e6:.1f}MB、806本なら時系列だけで約{bytes_full*806/1e9:.1f}GB。変分、複数条件、バックアップは別。\n'
(root/'計算時間.md').write_text(report)
# No FFT or spectrum: only the already computed time-domain traces.
q=results/'benchmark_double64'
x=np.loadtxt(q/'high_cos.tsv'); y=np.loadtxt(q/'variations.tsv')
fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
ax[0].plot(x[:,0],x[:,1]);ax[0].set_ylabel('Vector potential A(t)')
ax[1].plot(x[:,0],x[:,2]);ax[1].set_ylabel('Current J(t)')
ax[2].plot(y[:,0],y[:,1],label='cos');ax[2].plot(y[:,0],y[:,2],label='sin');ax[2].set_ylabel('Quadratic variation');ax[2].legend()
ax[2].set_xlabel('Time')
for a in ax:a.grid(alpha=.25)
fig.suptitle('Rice-Mele timing test | Double64 | Nk=800 | dt=0.01 | tau_g=3')
fig.savefig(results/'benchmark_timeseries.png',dpi=160)
print(report[:1600])
