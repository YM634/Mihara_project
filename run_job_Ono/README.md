# 小野法のテスト用ジョブ

既存の `Execute_sc/sc_model.ex` を直接使用。元ソースと元ジョブは変更しない。Aを直接与え、E=−dA/dtを解析微分して1列17桁で保存する。Gaussianの定義は参照論文に合わせてexp(−t²/(2tau²))。

## 用意する外場

ゼロ場、f、f±epsilon g_cos、f±epsilon g_sin。±は中央差分とepsilon²成分の検査用。偏光は既定でx方向。ノートの偏光は未確定なので、これを図5の確定条件とは扱わない。

`metadata.json` に時刻原点・規格化・実行ファイルを保存。C++の時刻s=0から開始し、物理時刻t=s+tminと読み替える。k格子は[0,2pi)のL×L。U_el=−4、T=0.02、mu=0、g_elph=0、S波。

## 実行例

作業ディレクトリはこのREADMEの親（研究補助）。Pythonの生成コードは標準ライブラリのみ、解析はNumPyが必要。

短時間の動作テスト（約40秒、出力先が存在したら別名にする）:

```bash
python3 run_job_Ono/run_job.py --output run_job_Ono/test_new --L 50 --tmin -7 --tmax 7 --tau-g 1 --dt 0.02 --F0 0.001 0.01 --omega 1 --run
python3 run_job_Ono/analyze.py run_job_Ono/test_new
```

短時間テストはtau_g=1に短縮しており、図5のスペクトルを再現する設定ではない。

論文に沿うGaussian幅、学生と同じ観測時間で1周波数の入力を生成（計算は実行しない）:

```bash
python3 run_job_Ono/run_job.py --output run_job_Ono/production_w1 --omega 1 --dt 0.02
```

生成された各フォルダの `sh job.sh` で実行する。全ケースの生成と逐次実行には `--run` を追加する。Nt=70000区間、70001点となり、1ケース約6分、既定の11ケースで約1時間の外挿。時間刻みの収束は別途必要。

偏光を単位対角方向にする場合は `--polarization 1 1`。自動的にsqrt(2)で正規化する。Ax=Ay=aを意味する指定ではない。

クラスターではLinux用に同じソースを再コンパイルして `--exe /path/to/sc_model.ex` を指定し、各job.shをそのクラスターのスケジューラから実行する。現在のMac arm64実行ファイルはLinuxでは使えない。スケジューラ種別は未指定のため、特定キューへの投入はしていない。

## 応答の解析

```bash
python3 run_job_Ono/analyze.py run_job_Ono/production_w1 --omega2 -5 5 201
```

`chi_bar_*.dat` は小野論文の対称化応答bar chi。通常の対称な展開係数chiを比較するなら1/2にする。Delta=−U_el Phi_scで換算。FTはexp(+i omega t)、台形積分、時刻原点補正、omega=omega1+omega2を直接評価する。既定は窓なし。`--window-tau VALUE` でGaussian窓を指定できるが、値は学生ノートに未記載。

この解析は反転対称性で一次応答が消えることを利用した二次の近似推定。有限F0の四次以上、有限tau_g、有限時間の誤差を含む。一般の一次応答がある模型向けの次数分離は実装していない。2つのF0の結果が一致するか、epsilonを変えて安定かを確認する。

`benchmark_L10`, `benchmark_L50`, `scaling_L50` は実測・検証用の短時間計算。長時間図5再現データではない。生成時に既存ディレクトリを上書きしない。
