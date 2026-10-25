# 平均場固定・精度切替の検証

- 基本検証16項目、固定平均場7項目、Double64と106-bit BigFloatの比較8項目、計31項目に成功。
- 固定平均場では独立な2×2密度行列RK4と全状態を比較し、最大差2e−14未満を確認。観測OPは変化するが、HamiltonianのOPは平衡値を使用する。
- 自己無撞着・固定の両モードでDouble64と106-bit BigFloatの最終状態の最大差1e−26未満を確認。
- Float64/Double64 × 自己無撞着/固定の4通りを、L=50、dt=0.01、Nt=1400の設定ファイルから実行。全て1401行・10列で有限値、timing.tomlの精度とモードを確認。
- 既定の自己無撞着Float64は改修前のfirst_juliaの時系列と1e−14未満で一致。
- 出力先: options_fixed_float64、options_fixed_double64、options_self_consistent_float64、options_self_consistent_double64。
- Double64例では平衡閾値を明示的に1e−24とした。高精度化しても時間刻みや格子の誤差は別に残る。
- 2D走査は実施していない。
