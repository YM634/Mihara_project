# Ono論文再現

小野論文の非相互作用Rice–Mele模型をJuliaで実装。初期ベンチマークに加え、指定された強めの外場・倍精度で図2(d),(e)に対応する二次応答を計算するコードを追加した。最新結果は `results/fig2_float64_stronger/` を参照。

- `src/RiceMele.jl`: 外場、瞬間基底状態へのRTA、RK4、電流、次数分離、単一周波数FT関数
- `configs/benchmark.toml`: 800 k点、dt=0.01、t=−20..20、tau_g=3
- `configs/benchmark_long.toml`: 時間区間2倍のスケーリング確認
- `configs/paper.toml`: 本文図2のスケール（観測区間は本実装で選択、未実行）
- `scripts/benchmark.jl`: コンパイル、時間発展、保存を分けて測定
- `tests/`: 行列実装との照合、ゼロ場、次数分離、RK4収束、高精度比較
- `results/`: 実測データと時間集計
- `再現条件.md`: 数式・出典・論文と計測条件の違い
- `計算時間.md`: 実測値と長時間・周波数走査の見積もり

## 実行

研究補助フォルダから実行。Julia 1.11.2、DoubleFloats 1.11.2で確認済み。依存パッケージはこのフォルダの `.julia_depot/` に導入済み。

```bash
sh Ono論文再現/scripts/run.sh Ono論文再現/configs/benchmark.toml Double64 Ono論文再現/results/my_benchmark
```

出力先は新しい名前を指定する（既存結果は上書きしない）。Float64ならDouble64をFloat64に変更する。`timing.toml`に6本の所要時間、保存時間、割当メモリ累積量（ピークRAMではない）、JITウォームアップ、条件を記録する。`*.tsv`は時刻・A・J、`variations.tsv`は二次の混合変分。スペクトルへの変換は行わない。

別環境に移した場合の依存関係復元:

```bash
JULIA_DEPOT_PATH="$PWD/Ono論文再現/.julia_depot" julia --project=Ono論文再現 -e 'using Pkg; Pkg.instantiate()'
```

## 検証

```bash
julia --startup-file=no --project=Ono論文再現 Ono論文再現/tests/runtests.jl
JULIA_DEPOT_PATH="$PWD/Ono論文再現/.julia_depot" julia --startup-file=no --project=Ono論文再現 Ono論文再現/tests/precision.jl
```

12個+2個の検証が成功。直接の2×2複素密度行列によるRK4との電流一致、ゼロ場の電流、時間刻みの4次収束、既知の多項式の次数分離、106-bit BigFloatとDouble64の比較を含む。初期の基本テストに加え、摂動論S55を独立のリゾルベント展開と照合する `scripts/check_spectrum.py` を追加した。時間発展との比較結果は最新の結果フォルダに記録する。

## 指定振幅・倍精度での図2(d),(e)計算

`results/fig2_float64_stronger/` にF_low=1e−3、F_high=1e−2、epsilon=1e−3の結果を保存。epsilonの二次成分を除く中央差分を使用する。元の計測コードと条件は残している。

```bash
JULIA_NUM_THREADS=4 julia --startup-file=no --project=Ono論文再現 Ono論文再現/scripts/calculate_fig2.jl
python3 Ono論文再現/scripts/check_spectrum.py
python3 Ono論文再現/scripts/plot_fig2.py
python3 Ono論文再現/scripts/write_fig2_report.py
```

再実行時は完了済みの行データを再利用する。`results/fig2_float64_stronger/結果.md` に条件・検証値・制限を記録する。作図にはNumPy、SciPy、Matplotlib、Python 3.10ではtomliが必要。
