# クラスター移行のための引き継ぎ：研究補助

更新日：2026-10-23。パスは特記しない限り、このファイルを置いた「研究補助」ディレクトリからの相対パス。

## 1. 目的とクラスターでこれから行いたいこと

三原さんが計算した、2次元引力Hubbard模型の超伝導秩序変数の外場に対する二次応答を検証する研究である。小野さんの論文の、時間発展から二周波数の非線形応答を抽出する方法を使用している。

**今後クラスターで行いたいことは、Julia版を実際に動かし、パラメータや計算条件を変えてさまざまな計算を行うこと。** コードを置くだけの移行ではなく、Julia実行環境の構築、動作確認、ジョブ投入、計算結果の保存・解析までを進めたい。

Macでは既に、元のC++コードとの整合性確認、Juliaの精度・並列化の確認、x方向外場による二次応答の2Dスキャンまで完了している。クラスターでは、まずこの既存結果と整合することを確認し、それから計算を拡張する。

現時点で候補となる計算は次のとおり。以下は今後の候補であり、すべて実行済みという意味ではない。

1. 格子数L、dt、時間範囲を変え、二次応答の収束を調べる。
2. 周波数刻みを細かくし、特に対角線上の共鳴付近を詳しく計算する。
3. 平均場更新あり／固定を比較し、集団的なフィードバックの寄与を調べる。
4. 外場振幅F、epsilonを変え、摂動領域と丸め誤差を確認する。
5. 必要な代表点でDouble64とFloat64を比較する。
6. 外場の方向、U、温度、化学ポテンシャルなどを変える。
7. 保存した時系列に異なる窓を使い、解析方法の依存性を調べる。

具体的な次の本番条件、クラスターのキュー・資源指定はまだ決めていない。計算候補を一括で無条件に投入する段階ではない。

## 2. 最初に押さえるべき現在の状態

- 現行の主計算コードは `Attractive_Hubbard/src/AttractiveHubbard.jl`。
- Julia版は元C++の物理模型を踏襲し、s波・電子フォノン結合なしの場合に特化して最適化している。
- 平均場更新あり／固定、Float64／Double64、k点のスレッド並列に対応する。
- 最新の2Dスキャンは **Float64、平均場更新あり、x方向、L=50、2プロセス×各1スレッド** で実行済み。
- 最新の生時系列は `Attractive_Hubbard/results/scan_x_float64/`。計810条件、約9.18 GB。
- 解析図のGammaまたはetaは、Hubbard計算では**後処理の窓関数の減衰率**。運動方程式の緩和率ではない。
- 図によって対称化した和chi_barと、その半分のchiを使い分けている。後述の規約を必ず確認する。
- 学生に渡す1周波数分は `For_Mihara/` に整理済み。
- この引き継ぎ作成時点で、ルートはGitリポジトリではない。GitHubリポジトリ作成、commit、push、クラスターへのcloneはまだ行っていない。

`Attractive_Hubbard/README.md` には開発初期の「2D走査はしていない」等の記述が残る。これは過去の作業範囲の説明であり、現状は上記のとおり2Dスキャン完了済み。最新条件は `configs/scan_x.toml` と結果ディレクトリの設定を参照する。

## 3. ディレクトリと文献

| 場所 | 内容・位置づけ |
|---|---|
| `状況説明.md` | 2026-10-17時点の当初依頼、元C++コードの所在 |
| `ノート分析.md` | 学生ノートの不整合、符号・規格化・不足条件の整理 |
| `onoproject_tex-2.pdf` | 学生のノートと結果 |
| `関連論文/qsxr-c2pq.pdf` | 小野さんの参照論文本文 |
| `関連論文/supplemental.pdf` | 同論文のsupplement |
| `関連論文/Higgs_mode.pdf` | 関連資料 |
| `run_job_Ono/` | 元C++を起動するPythonコード、初期テスト、時間見積もり |
| `Ono論文再現/` | 参照論文のRice–Mele模型をJuliaで再現したコード・結果 |
| `Attractive_Hubbard/` | 今後クラスターで主に動かしたいHubbard専用Juliaコード |
| `For_Mihara/` | 学生との照合用データ、共通README、条件別Params.md |
| `For_Mihara.zip` | 学生配布用の既存zip。README整理前に作成されている可能性があり、最新配布時は内容を確認して作り直す |
| `tmp/` | PDF抽出テキスト等の作業ファイル。計算コードそのものではない |

Mac上の元C++ソースは以下にある。この作業ディレクトリの外にあるため、研究補助全体をcloneしても元ライブラリ一式は付いてこない。

```text
/Users/yutamurakami/Documents/Softwares/Lib_TDDM_ver2/Programs_MF/SC_phonon/
```

`run_job_Ono/run_job.py` が参照する `Execute_sc/sc_model.ex` の元ソースをたどってJulia化した。C++の通常版mainは `sc_dynamics_main.cpp`、固定平均場版mainは `sc_dynamics_main_fixMF.cpp`。

## 4. ここまでの経緯

### 学生ノートと元C++の確認

学生ノートのハミルトニアン・秩序変数の符号、複素共役、Gaussian幅、chiとchi_barの因子2などを検討した。詳細は `ノート分析.md`。ノートの時刻・格子・周波数の記述には不一致もあり、ノートの図と全条件が一意に対応するとまでは確認できていない。

元C++では、外場ファイルに非常に小さな非正規化数（例：5e-324）が入ると `ifstream >> double` が失敗し、その後の読み込みが止まる問題を確認した。C++比較時は両実装でabs(field)<1e-300を0にした共通外場を使用した。元C++ソースを一括修正したわけではない。

診断：`Attractive_Hubbard/results/input_reader_diagnostic/`。

### 小野論文の再現

`Ono論文再現/` で非相互作用Rice–Mele模型を実装し、図2(d),(e)に対応する二次応答を計算した。これはHubbard模型とは別の模型。

代表条件はNk=800、dt=0.01、時間-1000..1000、hx=1、hy=0.5、hz=0.1、RTAのGamma=0.01、Float64。F_low=0.001、F_high=0.01、epsilon=0.001を使い、cos/sinの正負中央差分と2振幅の次数分離を行った。

このRice–Mele実装には瞬間基底状態へのRTAがある。一方、現在のHubbard実装にはないので混同しない。

結果：`Ono論文再現/results/fig2_float64_stronger/`。摂動論との相対L2差は2Dで約3.5%、対角線で約5.7%。弱い応答の領域、特に原点では相対誤差が大きく、完全一致ではない。

### Hubbard専用Julia版とC++比較

同じ模型をBlochベクトル3実数で記述し、共通する外場・三角関数計算をまとめ、配列を再利用した。汎用C++の2×2複素行列処理などとは実装上の仕事量が異なるので、速度差を言語だけの差として解釈しない。

平均場更新あり／固定の両方で、C++の初期状態を共有したJuliaとの一致を確認した。x方向に加え、x軸から22度の外場でもGaussian、cos、sinの3ケースを比較した。

22度の比較では、同じ初期状態からのOP差は約1e-15、電流差は約2e-15、エネルギー差は約2e-14以下だった。独立に平衡を求めた場合は、初期反復の残差が追加される。

結果の所在：

- `Attractive_Hubbard/results/cpp_comparison/比較結果.md`
- `Attractive_Hubbard/results/cpp_fixed_comparison/比較結果.md`
- `Attractive_Hubbard/results/angle22_comparison/比較結果.md`

固定平均場のC++出力とJulia出力ではエネルギーの定数項の規約が異なる。比較では既知の定数差を補正している。

### 精度と並列化

Float64、DoubleFloats.Double64、単独計算用の106-bit BigFloatを実装した。高精度では小数文字列から型を生成し、Double64の結果をFloat64へ落とさずに保存する。

k点のRK4更新はスレッド並列化済み。平均場と観測量の総和は逐次の同じ順序を維持する。1、2、4、8スレッドで比較し、検証した全時刻の観測量と最終k状態はビット単位で一致した。独立の2×2行列RK4との比較もある。

ただし、ビット単位一致は同一実装・同一条件での検証結果であって、別CPU・Juliaバージョンでも必ず保証するものではない。

## 5. 現在のHubbard模型と規約

2次元正方格子、s波、Hartree項なし、フォノンなし、緩和なし。相互作用Uは負。

$$
\xi_k=-2J(\cos k_x+\cos k_y)-\mu,
\qquad
h_k(t)=\begin{pmatrix}
\xi_{k-A(t)} & U\Phi(t)\\
U\Phi(t)^* & -\xi_{k+A(t)}
\end{pmatrix}.
$$

$$
\rho_k=\frac{I+\mathbf r_k\cdot\boldsymbol\tau}{2},
\qquad
\Phi=\frac1{2N_k}\sum_k(r_{kx}-ir_{ky}),
\qquad \Delta=-U\Phi.
$$

$$
\dot{\mathbf r}_k=2\mathbf h_k\times\mathbf r_k,
\quad
\mathbf h_k=\left(U\operatorname{Re}\Phi,-U\operatorname{Im}\Phi,
-2J[\cos k_x\cos A_x+\cos k_y\cos A_y]-\mu\right).
$$

最新条件のU=-4ではDelta=4 Phi。保存列のPhiをDeltaと取り違えない。二次応答の現在の観測量は実スカラーQ=-U Re(Phi)で、|Phi|をそのまま変換しているわけではない。

k格子は各軸2*pi*j/L（j=0..L-1）、Nk=L²。2*pi端点は含めない。初期状態はA=0の自己無撞着熱平衡で、正の実数のPhiを選ぶ。平衡反復の停止条件は連続する密度行列のFrobeniusノルム差のk平均<tolerance。

Float64、L=50、U=-4、T=0.02、mu=0、tolerance=1e-12では、最新スキャンの初期値はおおよそ

$$
\Phi_{\rm eq}=0.3453269522107971,
\qquad \Delta_{\rm eq}=1.381307808843188.
$$

時間積分は陽的古典RK4。更新ありでは各RK段階でPhiを再計算する。`update_mean_field=false` は、HamiltonianのPhiだけを平衡値に固定する。A(t)による運動項の時間依存と、観測するPhi(t)の時間変化は残る。

固定時のエネルギーは以下の規約。

$$
E_{\rm int,fixed}=2U\operatorname{Re}[\Phi_{\rm eq}^*\Phi(t)]-U|\Phi_{\rm eq}|^2.
$$

## 6. 完了した最新2Dスキャン

設定：`Attractive_Hubbard/configs/scan_x.toml`。実行時にコピーした設定は結果ディレクトリの `scan_config.toml`。

| 項目 | 値 |
|---|---|
| J, U, T, mu | 1, -4, 0.02, 0 |
| L, Nk | 50, 2500 |
| dt | 0.01 |
| 時間範囲 | -700..700 |
| 時間ステップ数Nt | 140000 |
| 保存点数 | 140001（両端を含む） |
| 平衡tolerance | 1e-12 |
| 精度・平均場 | Float64、自己無撞着更新あり |
| F_low, F_high, epsilon | 0.001, 0.01, 0.001 |
| tau_f, tau_g | 0.1, 100（いずれも標準偏差） |
| 偏光 | x方向、(1,0) |
| 外場補間 | table4：整数時刻の表から5点4次Lagrange補間 |
| omega1 | 0..5、0.05刻み、101点を時間発展 |
| omega2 | -5..5、0.05刻み、201点を後処理 |
| 完成マップ | omega1,omega2とも-5..5の201×201点 |
| 実行方法 | 2個の独立Juliaプロセス、各1スレッド |

外場は電場Eではなく、ベクトルポテンシャルAを直接与える。

$$
f_F(t)=\frac{F}{2\pi}e^{-t^2/(2\tau_f^2)},
\qquad g_F(t)=\frac{F}{2\pi}e^{-t^2/(2\tau_g^2)},
$$

$$
A_x(t)=f_F(t)\pm\epsilon g_F(t)C(\omega_1t),
\quad C=\cos\text{ または }\sin,
\qquad A_y(t)=0,\quad E=-\partial_t A.
$$

学生ノートはexp[-(t/tau_note)²]という幅の定義。今回のtau_g=100はtau_note=sqrt(2)/0.01に対応する。`configs/one_frequency.toml` は以前の標準偏差141.421...の設定なので、そのまま使うと最新scan_xとは同じ条件にならない。

各omega1につき2振幅×cos/sin×正負=8本。これに周波数に依存しないf単独2本を加えて810本。omega1=0のsin条件もそのまま保存した。負のomega1の時間発展はしておらず、実観測量の共役対称性から復元した。

Macでの実行は2026-10-23 14:33:20から15:25:53、約52分33秒。初期に約7.4秒／1本、その後概ね7.6～8秒／1本。これはApple M2上の2プロセス同時実行での実測であり、クラスターの所要時間は改めて測る。

### 保存データの読み方

出力ルート：`Attractive_Hubbard/results/scan_x_float64/`。

例：omega1=1.40は `omega_028_F0.001_cos_plus/`。番号は0始まりでomega1=番号×0.05。`common_F0.001_f/` 等はf単独の共通条件。

各条件には以下を保存した。

- `observables.f64`：ヘッダなしlittle-endian Float64、140001行×10列、Julia列優先配列。
- `initial_state.f64`, `final_state.f64`：同形式の3行×2500列。
- `config.toml`：個別条件。
- `complete.toml`：形状、列名、完了時刻、計算時間等。これがあれば保存完了扱い。

```text
1 time
2 n_up
3 n_down
4 Phi_re
5 Phi_im
6 Jx
7 Jy
8 Ekin
9 Eint
10 Etot
```

Pythonでの読み込み：

```python
import numpy as np
v = np.fromfile('observables.f64', dtype='<f8')
obs = v.reshape((140001, 10), order='F')
t = obs[:, 0]
phi = obs[:, 3] + 1j * obs[:, 4]
delta_real = 4 * obs[:, 3]
```

全観測量は約9.07 GB、全条件の初期・最終k状態を加えて約9.17 GB。全k点の全時刻の状態は保存していない。これにより全体を約10 GB以内に収めた。解析結果を増やすと追加容量は必要。

## 7. 二次応答の後処理

解析の主コード：`Attractive_Hubbard/scripts/analyze_scan.py`。

まずQ=-U Re(Phi)を使い、外場の正負中央差分を作る。

$$
D_F^a(t)=\frac{Q[A_F^{a,+}](t)-Q[A_F^{a,-}](t)}2,
\qquad a=c,s.
$$

epsilonではまだ割らない。振幅比r=F_high/F_low=10として

$$
d_a(t)=\frac{D_{\rm high}^a(t)-rD_{\rm low}^a(t)}{1-1/r}
=\frac{D_{\rm high}^a(t)-10D_{\rm low}^a(t)}{0.9}.
$$

低次数を消した二次推定であり、高次の混入が完全に消えるわけではない。

$$
W_\eta(t)=e^{-\eta\max(t,0)},
\qquad q_\eta(t)=W_\eta(t)[d_c(t)-i d_s(t)],
$$

$$
\widetilde q_\eta(\Omega)
=\int_{-700}^{700}q_\eta(t)e^{+i\Omega t}dt,
\qquad \Omega=\omega_1+\omega_2.
$$

実装はdt=0.01の台形則（両端の重み1/2）。等間隔omega2の変換にSciPyのZoomFFTを使用する。時間原点tmin=-700の位相も掛けている。非等間隔点には直接積分も使用可能。

$$
\widetilde f_{\rm high}(\omega_2)
=\frac{F_{\rm high}\tau_f}{\sqrt{2\pi}}
 e^{-\omega_2^2\tau_f^2/2},
$$

$$
\bar\chi^{(2)}_\eta(\omega_1,\omega_2)
\simeq\frac{2\pi\widetilde q_\eta(\omega_1+\omega_2)}
{\epsilon F_{\rm high}\widetilde f_{\rm high}(\omega_2)}.
$$

有限幅のgを単色波として近似しており、有限幅の誤差は残る。

$$
\bar\chi^{(2)}(\omega_1,\omega_2)
=\chi^{(2)}(\omega_1,\omega_2)+\chi^{(2)}(\omega_2,\omega_1).
$$

対称な展開係数を選ぶ規約ならchi=chi_bar/2。chi_barは平均ではなく**和**。入力はAなので、ここで得たものをそのまま電場Eに対するchiと呼ばない。

$$
\chi^{(2)}(-\omega_1,-\omega_2)
=\chi^{(2)}(\omega_1,\omega_2)^*.
$$

これで負のomega1側を復元する。対角線では実部が偶、虚部が奇になる。

### Gamma・窓と作成済みの図

HubbardでGamma=0.01、0.02、0.05、0.5と呼んだものは、上記etaと同じ**後処理の減衰率**。t<0はW=1、t>0は指数減衰。物理的RTAの追加とは違う。Gammaを大きくすると遅い時刻の寄与を抑え、細かな構造が平滑になる一方、ピーク幅・高さが変わる。

`results/scan_x_float64/` の主な解析成果物：

| ファイル | 規約・内容 |
|---|---|
| `spectra.npz` | 窓なし／eta=0.01、chi_barとchiの両方、周波数軸 |
| `chi2_map.png` | 窓なしのchi_barの絶対値・実部・虚部 |
| `chi2_abs.png` | eta=0.01のabs(chi_bar) |
| `chi2_map_eta001.png` | eta=0.01のchi_bar、3成分 |
| `chi2_diagonal_abs.png` | 対角線のabs(chi)、窓なしとeta=0.01 |
| `chi2_diagonal_windows.png` | 対角線のabs(chi)、eta=0.01,0.02,0.05の比較 |
| `diagonal_window_comparison.npz` | 上記の複素chi、eta=0も含む |
| `chi2_diagonal_real_imag_eta005_full.png` | eta=0.05、正負の対角線、chiの実部・虚部 |
| `chi2_map_gamma0p05.png`, `.npz` | Gamma=0.05の2D、chiの絶対値・実部・虚部 |
| `chi2_map_gamma0p5.png`, `.npz` | Gamma=0.5の2D、同じ規約 |
| `analysis.json`, `結果.md` | 診断値と条件 |
| `source_snapshot/` | スキャン・主解析のソース控え。一部は後から追加したため、全後処理を含むわけではない |

結果はsum frequencyが約±2 Delta_eq（約±2.763）の近傍で強い応答を示す。これは計算結果の特徴であって、全ての細かなピークが物理的に収束したと確定したわけではない。

### 検証済みの数値と未完了の確認

- 810条件すべての保存完了とバイナリファイルサイズを確認済み。
- ZoomFFTと独立の直接積分の最大相対差：約3.97e-10。
- omega1,omega2交換に対する相対L2差：窓なし約0.64%、eta=0.01では約0.041%。
- 高低振幅のF²スケーリング相対L2差の周波数中央値：約0.41%。
- 代表パルス単独データでパルス後のエネルギー変動は約3e-14。代表例の密度変化は約4e-16。
- 最新2D全体について、L、dt、時間区間、epsilonを変えた系統的収束試験はまだない。
- 後処理の強い窓で曲線が滑らかになったことは、格子・時間積分の収束証明にはならない。

## 8. 学生向けデータ：For_Mihara

まずはomega1=1.40の1セットを渡す方針。対角線では2omega1が2Delta_eqに近く、強い応答を確認できる。必要なら次にomega1=0.50を追加するが、現時点で配布ディレクトリを作成したのは1.40のみ。

8ディレクトリそれぞれに以下がある。

- `observables.dat`：時間、Phi実部、Phi虚部、Jx、Jyの5列、140001行、有効数字17桁。
- `Params.md`：F、cos/sin、正負、外場係数の数値など個別条件。

共通の外場形・Hamiltonian・模型・dt・格子・観測量の定義は `For_Mihara/README.md` に集約した。現在は各Params.mdからREADMEを参照する構成。生時系列であり、窓やFTは掛けていない。

全8条件は元Float64の値とテキスト再読み込み値が一致することを確認した。合計約132 MB。`scripts/export_mihara.py` はこの形式に合わせて更新済みだが、既存出力を上書きしない。再書き出し時は出力先・バックアップ方針を先に決める。

## 9. クラスターでの環境構築・最小実行

### 環境と依存関係

Macで確認した環境はJulia 1.11.2、DoubleFloats 1.11.2。Project.tomlとManifest.tomlを各Juliaプロジェクトに保存してある。Python解析はNumPy、SciPy、Matplotlibを使用する。Python 3.10以下ではtomliも必要。クラスターではまず同じJulia系列で再現確認する。

Juliaだけの実行にMacのC++バイナリや元C++ライブラリは不要。Double64を使う場合はDoubleFloatsの依存関係が必要。

Macにある `Ono論文再現/.julia_depot/` はパッケージキャッシュ・コンパイル済み成果物を含む。これをLinuxへそのまま移すのではなく、クラスターでinstantiateする。クラスターが外部ネットワークへ出られない場合は、ログインノード等で事前に依存関係を用意する必要がある。

以下はルートで実行する例。moduleの名前やJuliaの導入方法はサイトに合わせる。

```bash
# サイトに応じて Julia 1.11 を使えるようにしてから実行
julia --version
export JULIA_DEPOT_PATH="${SCRATCH:-$HOME}/.julia_mihara"
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=Attractive_Hubbard -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=Ono論文再現 -e 'using Pkg; Pkg.instantiate()'
```

Pythonの仮想環境例：

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install numpy scipy matplotlib tomli
```

環境を構築したら、まず小さいテストと短い1条件を実行する。

```bash
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/tests/runtests.jl
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/tests/precision.jl
sh Attractive_Hubbard/scripts/run.sh Attractive_Hubbard/configs/comparison.toml Attractive_Hubbard/results/cluster_smoke Float64
```

`comparison.toml` は短時間・強い外場の検証用で、本番のscan_xと同じ条件ではない。単独 `run.jl` はフラットなTOML、スキャンは `[model]` と `[scan]` のあるTOMLを読む。`scan_x.toml` をそのまま単独run.jlへ渡すことはできない。

単独計算の出力先は既存でないものを指定する。Double64の例は `configs/self_consistent_double64.toml`、`configs/fixed_mean_field_double64.toml`。Double64を選んでも平衡toleranceは自動で厳しくならない。用意した高精度例では1e-24を明記している。

### スキャンをクラスターで実行する方法

`Attractive_Hubbard/scripts/launch_scan.py` はMac向けであり、Juliaの実行パス `/opt/local/bin/julia` が直書きされている。クラスターでそのまま実行しない。また、デポの指定をローカルの共有キャッシュに上書きする実装である。

クラスターではPATH上のjuliaで `scan_worker.jl` を直接起動するか、移植版ランチャーを用意する。workerのCLIは以下。

```text
scan_worker.jl CONFIG OUTPUT WORKER_ID WORKER_COUNT
```

workerはomega1のインデックスをWORKER_ID:WORKER_COUNT:101に分配する。worker1だけがf単独2本も担当する。各workerは起動を維持して複数の条件を順番に処理する。現在のworkerはFloat64・1スレッド専用で、全スキャンをDouble64へ自動切替する実装ではない。Double64スキャンにはworker側の精度対応が別途必要。

例えば、ジョブスケジューラーで確保した同一ノードの2 CPUに対し、以下のように実行できる。これは計算ノード内のシェル用の例であり、ログインノードでの長時間実行を意図しない。CPU binding等はスケジューラーの設定に合わせる。

```bash
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# JULIA_DEPOT_PATHは上で設定したクラスター用のものを使う
cluster_out="$PWD/Attractive_Hubbard/results/cluster_scan_x"
mkdir -p "$cluster_out"
cp Attractive_Hubbard/configs/scan_x.toml "$cluster_out/scan_config.toml"
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/scan_worker.jl "$cluster_out/scan_config.toml" "$cluster_out" 1 2 > "$cluster_out/worker1.log" 2>&1 &
worker1_pid=$!
julia --startup-file=no --project=Attractive_Hubbard Attractive_Hubbard/scripts/scan_worker.jl "$cluster_out/scan_config.toml" "$cluster_out" 2 2 > "$cluster_out/worker2.log" 2>&1 &
worker2_pid=$!
scan_status=0
wait "$worker1_pid" || scan_status=1
wait "$worker2_pid" || scan_status=1
if [ "$scan_status" -ne 0 ]; then
    echo 'Scan failed; inspect worker logs.' >&2
    exit "$scan_status"
fi
python Attractive_Hubbard/scripts/analyze_scan.py "$cluster_out"
```

スキャンworkerは各条件のcomplete.tomlがあればスキップする。保存は一時ファイルからrenameしている。ただし、再開時に全個別条件の設定を自動比較する仕組みはないので、**同じ出力先に別の設定を混ぜない**。上のコピー例は新規出力先用。再開時は保存済みscan_config.tomlを使い、変更した設定を上書きしない。計算条件を変える場合は必ず新しい出力ディレクトリを使う。

スケジューラー（Slurm/PBS等）、partition、walltime、CPU数、メモリ指定、配列ジョブの使用可否は未確認。workerはMPI通信を使わないため、別ノードへの独立割り当ても可能だが、同じ出力ルートを参照できる共有ファイルシステムか、集約方法が必要。実際のsbatch/qsubスクリプトはサイトの仕様を確認して作る。

### 後処理の再実行

```bash
python Attractive_Hubbard/scripts/analyze_scan.py Attractive_Hubbard/results/scan_x_float64
python Attractive_Hubbard/scripts/compare_diagonal_windows.py
python Attractive_Hubbard/scripts/plot_scan_window.py --gamma 0.05
python Attractive_Hubbard/scripts/plot_scan_window.py --gamma 0.5
```

`analyze_scan.py` は出力ディレクトリを引数に取る。`--watch` を付けると、保存完了した周波数から順次変換し、残りを待つ。

`compare_diagonal_windows.py` と `plot_scan_window.py` は現在 `results/scan_x_float64/` を固定参照し、格子・振幅・dt等にも今回の値が直書きされている。別条件を計算する前に、TOMLから読み込むようにするか、対象条件に合わせて修正する必要がある。`analyze_scan.py` は多くをTOMLから読むが、列・ファイル名・810条件の想定等も今回のスキャンに依存する。任意模型・任意振幅に対応した汎用解析器と考えない。

PNG/SVG作図にGPUは不要。ここまでのJulia実装はCPU用であり、GPUやMPIによるk点分割は実装していない。

## 10. 並列化と速度の判断材料

Macのスレッド検証はApple M2上、Juliaが認識したCPU数5で行った。時間発展1400ステップ、全時刻の観測量評価込み、JIT・平衡・保存を除く3回中央値。

| 条件・平均場 | 1スレッド | 4スレッド |
|---|---:|---:|
| Float64 L=50 更新 | 0.06784秒 | 0.06337秒 |
| Float64 L=50 固定 | 0.06487秒 | 0.04965秒 |
| Float64 L=200 更新 | 1.12041秒 | 0.51803秒 |
| Float64 L=200 固定 | 0.99808秒 | 0.33904秒 |
| Double64 L=50 更新 | 1.81532秒 | 1.12350秒 |
| Double64 L=50 固定 | 1.61102秒 | 0.85573秒 |

L=50のFloat64ではk点並列の効果が小さいため、本番は周波数ごとのパラメータ並列を採用した。大きいLやDouble64ではk点並列にも効果がある。クラスターではCPU・メモリ帯域・同時ジョブ数が違うので、同じ結論を無条件に適用しない。

詳細：`Attractive_Hubbard/results/thread_scaling/検証結果.md`。

複数プロセスと各プロセス内のスレッドを併用する場合は、プロセス数×スレッド数が確保CPU数を超えないようにする。現在のスキャンworkerは1スレッドを明示的に要求する。

## 11. GitHubで共有する際の範囲と移植上の注意

「研究補助全体を共有する」という目的に対しても、コード・文書・設定・小さい検証データと、再生成可能な約9 GBの生データ・キャッシュは保存経路を整理した方がよい。

現時点の概算容量はAttractive_Hubbard約9.24 GB、Ono論文再現約0.91 GB（デポを含む）、For_Mihara約0.132 GB。そのほかPDFやzip等がある。最新スキャンは1ファイルが巨大というより、多数の約11 MBファイルの合計が大きい。

提案する管理方針：

- Gitで共有：Julia/Pythonソース、TOMLとManifest、READMEとこの引き継ぎ、ノート分析、必要な小さな検証結果・図。
- 生時系列：クラスターの研究用ストレージへ別転送するか、運用を決めてGit LFS等を利用する。再計算するなら、どの結果がcloneに含まれないかをREADMEへ明記する。
- `.julia_depot/`、コンパイルキャッシュ、`__pycache__/`、`.venv/`、`.DS_Store`、一時出力はGit管理から外す。
- Macで作った`.ex`や`.dylib`等はLinuxクラスターでそのまま動くことを期待しない。Juliaソースを主に使い、C++照合が必要なら外部ライブラリを用意して再ビルドする。
- 元C++ビルド用のEigenパス `/opt/local/include/eigen3` と外部ソースパスもMac固有。クラスターでC++再比較するときだけ移植する。
- Unicode名・空白を含むパスがあるため、シェルではパス変数を引用符で囲む。

このファイルの作成では、ルートの.gitignore作成、Git初期化、LFS設定、外部へのpushやデータ転送はまだ行っていない。既存の `Ono論文再現/.gitignore` は `.julia_depot/` の除外のみ。大容量データを含めるか決めてからGitHub登録を進める。

## 12. クラスター側で最初に引き継いでほしい作業

1. このファイルと `For_Mihara/README.md` を読み、Phi/Delta、A/E、chi/chi_bar、窓/物理緩和の規約を把握する。
2. JuliaとPythonの依存関係を構築し、Hubbardのテストと短時間の単独計算を動かす。
3. omega1=1.40の代表条件を、最新scan_xと同じ条件で再計算し、Macの時系列・学生配布データと比較する。
4. スケジューラーに合わせた起動方法を整え、クラスターでの時間とメモリを測る。
5. 条件・出力先・再開方法を明確にし、ユーザーが指定した新しい計算を順次実行する。
6. 時系列を残して解析の再実行を可能にし、物理結果の変化と窓・離散化・丸め誤差の影響を分けて報告する。

本ファイルに記載したクラスターのコマンドは移行手順の例であり、クラスターでの実行・検証済みの記録ではない。今後は実際のサイト名、CPU、Juliaバージョン、投入スクリプト、資源指定、所要時間を追記していく。
