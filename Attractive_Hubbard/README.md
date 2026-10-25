# Attractive Hubbard 専用Juliaコード

既存の `run_job_Ono/run_job.py` が起動する `Execute_sc/sc_model.ex` の元ソースを参照し、2次元正方格子の時間依存BCS平均場をJuliaで実装した。フォルダは既存の `Attractive_Hubbard/` を使用。

`Ono論文再現/` と同様にTOML設定、精度の切替、JITウォームアップを除いた時間計測、時系列保存、検証スクリプトを用意。**今回の範囲はC++との結果・速度比較で、2D走査もスペクトル計算も行わない。比較結果は [結果レポート](results/cpp_comparison/比較結果.md) を参照。**

## 対象とする模型

2D attractive Hubbard、S波、電子フォノン結合なし。C++のSC_symmetry=S、g_elph=0に対応する。D波・フォノン自由度は移植対象外。固定平均場モードも選択できる。RTAや減衰は加えていない。Uは負の相互作用で、化学ポテンシャルの約束もC++に合わせる（Hartree項を追加しない）。

$$
h_k(t)=\begin{pmatrix}
\epsilon_{k-A(t)}-\mu&U\Phi(t)\\
U\Phi^*(t)&-\epsilon_{k+A(t)}+\mu
\end{pmatrix},\qquad
\epsilon_k=-2J(\cos k_x+\cos k_y).
$$

$$
\rho_k=\frac{I+\mathbf r_k\cdot\boldsymbol\tau}{2},\qquad
\Phi=\frac{1}{2N_k}\sum_k(r_{kx}-ir_{ky}),\qquad
\Delta=-U\Phi.
$$

$$
\dot{\mathbf r}_k=2\mathbf h_k\times\mathbf r_k,\quad
\mathbf h_k=\left(U\mathrm{Re}\Phi,-U\mathrm{Im}\Phi,
-2J[\cos k_x\cos A_x+\cos k_y\cos A_y]-\mu\right).
$$

単位行列部分は密度行列の交換子には寄与しないが、エネルギーの評価には含める。RK4の各段階でPhiを自己無撞着に更新する。k格子は[0,2pi)のL×L、端点重複なし。k配列の順序もC++と同じ。

平衡はA=0での有限温度の自己無撞着解。C++との比較では開始時のAが無視できるパルスを使う。初期Phiは実数正0.35から始め、各kの密度行列変化のFrobeniusノルムの平均がtolerance未満で終了する。C++の乱数初期値とは反復の開始点が異なるため、収束直後の微小残差は一致しない場合がある。比較スクリプトでは独立な平衡解の比較と、C++初期状態を共有した比較の両方を行う。

## 外場と出力

外場はベクトルポテンシャル。Gaussianのf、f±epsilon g_cos/sinを選択可能。偏光は指定ベクトルを単位長に正規化する。デフォルトの`table4`は整数時刻の外場を生成後、5点4次補間でRK4の中間時刻を評価し、C++のファイル入力条件を再現する。`analytic`は中間時刻の解析式を直接使う。二つは有限dtで完全には同じではない。

`observables.tsv` の列:

```text
time n_up n_down Phi_re Phi_im Jx Jy Ekin Eint Etot
```

時刻はtminからの物理時刻。C++の内部時刻sはs=t−tminに対応。Phiはギャップそのものではなく、Delta=−U Phi。Eint=U|Phi|²。C++と同じエネルギーの定義で比較する。

## 実行

研究補助フォルダから:

```bash
sh Attractive_Hubbard/scripts/run.sh Attractive_Hubbard/configs/comparison.toml Attractive_Hubbard/results/my_run Float64
```

`Double64`ならDoubleFloatsを使用。`BigFloat`は106ビットで動作する。小数入力は文字列から対象の型に変換し、高精度出力をFloat64に落とさない。既存の出力フォルダは上書きしない。`timing.toml`には平衡・時間発展・保存・初回ウォームアップを分けて記録する。

依存パッケージは `Ono論文再現/.julia_depot` の導入済み環境を共有する。Manifest.tomlを保存している。別環境ではJULIA_DEPOT_PATHを設定してPkg.instantiate()を実行する。

1周波数・2振幅・5外場のセットも生成・実行できる（**今回この本番セットは未実行**）:

```bash
JULIA_DEPOT_PATH="$PWD/Ono論文再現/.julia_depot" julia --project=Attractive_Hubbard Attractive_Hubbard/scripts/run_set.jl Attractive_Hubbard/configs/one_frequency.toml Attractive_Hubbard/results/one_frequency Float64
```

## C++との比較

```bash
python3 Attractive_Hubbard/scripts/build_reference.py
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/compare_cpp.jl
MPLCONFIGDIR=/private/tmp/codex_mpl_op python3 Attractive_Hubbard/scripts/plot_comparison.py
```

比較先は `results/cpp_comparison/`。同一条件の完了済みケースは再利用する。再計算時は新しい出力先に変更するか、既存結果を別名へ保存する。ビルドの参照ルートは引数で変更可能。コンパイラはclang++、Eigenは/opt/local/include/eigen3を使用。元ソース・元実行ファイルは変更しない。

比較用mainだけに17桁出力、初期状態のダンプ、時間計測を追加し、ハミルトニアン・平衡・RK4の元クラスはそのままコンパイルする。`reference/build.json` にコンパイルコマンドとソースのハッシュを保存。元の12桁出力の実行ファイルとも1条件で照合する。

## テスト

```bash
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/tests/runtests.jl
JULIA_DEPOT_PATH="$PWD/Ono論文再現/.julia_depot" julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/tests/precision.jl
```

独立な2×2密度行列RK4との一致、4次補間、既知の平衡ギャップ、ゼロ場、密度・電流、入力検証、高精度型間の比較を含む。

## 平均場固定とdouble-doubleのオプション

通常の自己無撞着計算（従来どおり、既定値）:

```toml
precision = "Float64"
update_mean_field = true
```

平衡OPで平均場を固定し、double-doubleで計算:

```toml
precision = "Double64"
update_mean_field = false
tolerance = "1e-24"
```

`precision`はrun.jlが読み取るトップレベルの設定。実行コマンド末尾の精度指定があれば、そちらが優先される。Double64を選んでも平衡の閾値は自動変更しない。用意したDouble64例では明示的に1e−24に設定した。dtや格子の精度は別途収束確認が必要。

```bash
# 平均場固定、倍精度
sh Attractive_Hubbard/scripts/run.sh Attractive_Hubbard/configs/fixed_mean_field.toml Attractive_Hubbard/results/my_fixed

# 平均場固定、double-double
sh Attractive_Hubbard/scripts/run.sh Attractive_Hubbard/configs/fixed_mean_field_double64.toml Attractive_Hubbard/results/my_fixed_dd

# 自己無撞着更新、double-double
sh Attractive_Hubbard/scripts/run.sh Attractive_Hubbard/configs/self_consistent_double64.toml Attractive_Hubbard/results/my_self_consistent_dd
```

固定するのはHamiltonianに使うPhi_eq。外場A(t)による運動項の時間依存は残す。時間発展後のPhi(t)自体は出力するので、固定平均場でもOPの変化を観察できる。固定モードではRK4の全段階に同じPhi_eqを使う。

固定モードのエネルギーは、固定された平均場Hamiltonianの期待値に、平衡エネルギーと一致させる定数を含める:

$$
E_{\mathrm{int,fixed}}(t)=2U\,\mathrm{Re}[\Phi_{\mathrm{eq}}^*\Phi(t)]-U|\Phi_{\mathrm{eq}}|^2.
$$

自己無撞着モードのEint=U|Phi(t)|²とは区別する。定数項は時間発展には影響しない。元のC++固定平均場版のEintはこの定数を引いていないので、絶対エネルギーには定数差がある。

両モードで、Float64、Double64、106-bit BigFloat（単独run.jl）が使用できる。Double64で生成・積分・保存した値をFloat64へ変換しない。`timing.toml`に実際の精度、update_mean_field、固定時のHamiltonianのOPを記録する。

1周波数セット用のrun_set.jlでは、トップレベルにprecision、[model]内にupdate_mean_fieldを指定する。2D走査は実行していない。

追加検証: 固定平均場と独立な2×2行列RK4の一致、OP期待値が実際に変化すること、固定Hamiltonianのエネルギー定義、両モードのDouble64/106-bit BigFloat一致を確認した。

## C++固定平均場版との照合

元の `sc_dynamics_main_fixMF.cpp` から固定平均場用C++参照実行ファイルをビルドし、Juliaの `update_mean_field=false` と比較できる。

```bash
python3 Attractive_Hubbard/scripts/build_reference.py --fixed
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/compare_cpp.jl --fixed
MPLCONFIGDIR=/private/tmp/codex_mpl_op python3 Attractive_Hubbard/scripts/plot_fixed_comparison.py
```

結果は `results/cpp_fixed_comparison/比較結果.md`。OP・密度・電流と、既知の定数差を合わせたエネルギーを検証。生のC++出力も保存する。このモードでは速度計測は各1回の参考値。自己無撞着版の3回計測結果とは別フォルダに保存し、既存結果を変更しない。

## 22度方向の外場での照合

```bash
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/compare_cpp.jl --angle=22
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/compare_cpp.jl --angle=22 --fixed
MPLCONFIGDIR=/private/tmp/codex_mpl_op python3 Attractive_Hubbard/scripts/plot_angle_comparison.py
```

任意の角度は--angle=角度（度）で指定できる。偏光を(cos角度,sin角度)とし、Gaussian、cos試験場、sin試験場を比較する。結果は角度・平均場モードごとに保存。22度の比較レポートと図はresults/angle22_comparison/。

## k点のスレッド並列化

2026-10-23追加。`JULIA_NUM_THREADS`でワーカ数を指定するとRK4のk点更新を並列化する。

```bash
JULIA_NUM_THREADS=4 sh Attractive_Hubbard/scripts/run.sh Attractive_Hubbard/configs/fixed_mean_field_double64.toml Attractive_Hubbard/results/my_parallel_dd
```

`parallel_k = false`をTOMLに指定すると、複数スレッド起動時も逐次計算する。既定はtrue。平衡計算とOP・観測量の総和は従来順序の逐次処理を維持し、RK段階間で待ち合わせる。

1・2・4・8スレッドで自己無撞着/固定、Float64のL=50,200とDouble64のL=50を各3回測定した。全時刻の全観測値と最終全k状態は、スレッド数に依存せずビット単位で一致した。改修前の22度の逐次結果も維持される。

今回の環境では4スレッドが概ね良好。小さいFloat64計算では待合せの費用が支配する場合がある。8スレッドは4より遅かった。詳しい時間と比較図は `results/thread_scaling/検証結果.md` を参照。

計測を再現する際は `scripts/benchmark_threads.jl` をスレッド数ごとに別プロセスで実行し、`scripts/compare_threads.jl` と `scripts/plot_thread_scaling.py` で集計する。既存出力先は上書きしない。
