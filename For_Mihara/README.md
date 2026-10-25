# 三原さん確認用：omega1=1.40の1セット

外場2強度（F=0.001,0.01）× cos/sin × 正負の計8条件です。
全条件に共通する定義・設定はこのREADME.mdに、条件ごとの振幅・波形・符号は各ディレクトリの[Params.md](omega1_1.40_F_0.001_cos_plus/Params.md)に記載しています。
まずは同じ条件でPhi(t),Jx(t),Jy(t)が一致するか確認してください。

## データの形式と秩序変数

各条件の`observables.dat`は、空白区切り、有効数字17桁のテキストです。

| 列 | 内容 |
|---|---|
| 1 | 時間 t |
| 2 | 秩序変数 Phi の実部 |
| 3 | 秩序変数 Phi の虚部 |
| 4 | サイト当たり電流 J_x |
| 5 | サイト当たり電流 J_y |

$$
\Phi=\frac1N\sum_k\langle c_{-k\downarrow}c_{k\uparrow}\rangle,
\qquad \Delta=-U\Phi=4\Phi.
$$

保存している秩序変数はDeltaではなくPhiです。
後処理の窓・フーリエ変換・平衡値の減算は適用していない生の時系列です。
全時刻を保存しているため、異なる窓やフーリエ変換を後から試せます。

## 模型・初期状態

2次元正方格子の引力Hubbard模型、s波超伝導平均場。
hopping J = 1
U = -4
temperature T = 0.02
chemical potential mu = 0
電子・フォノン結合なし。Hartree項なし。緩和項なし。
単位: hbar=1、格子定数=1。AはPeierls位相として無次元。
L = 50（各軸のk点数）
N = Nk = L*L = 2500（総k点数、サイト当たり平均の分母）
k_x,k_y = 2*pi*j/L, j=0,...,L-1（端点2*piは含めない）
k点の並び: kxが先に変化。
初期状態: 同じU,T,muの自己無撞着熱平衡状態。
平衡反復の収束閾値 = 1e-12
収束量: 連続する密度行列のFrobenius差のk平均。
初期Phiは正の実数の位相を選ぶ。
平均場Phi(t)はRK4の各段階で更新する（固定平均場ではない）。
初期Phi_re = 3.4532695221079707e-01
初期Phi_im = -0.0000000000000000e+00

## 時間発展

dt = 0.01
t_min = -700
t_max = 700
Nt = 140000（時間ステップ数）
保存点数 = Nt+1 = 140001（両端を含む、全ステップ保存）
積分法 = explicit classical RK4
数値精度 = Float64（倍精度）
1計算あたりのスレッド数 = 1

## 外場の共通定義

入力はベクトルポテンシャルAです。電場はE=-dA/dt。
方向はx軸、polarization=(1,0)、A_y=E_y=0です。

| パラメータ | 全条件で共通の値 |
|---|---|
| epsilon | 0.001 |
| omega1（角周波数） | 1.40 |
| tau_f（標準偏差） | 0.1 |
| tau_g（標準偏差） | 100 |

$$
f_F(t)=\frac{F}{2\pi}\exp[-t^2/(2\tau_f^2)],
\qquad g_F(t)=\frac{F}{2\pi}\exp[-t^2/(2\tau_g^2)].
$$

$$
A_x(t)=f_F(t)+s\,\epsilon g_F(t)C(\omega_1t),
\qquad A_y(t)=0.
$$

F、C（cosまたはsin）、s（+1または-1）は各Params.mdを参照してください。
`cos_plus`はC=cos,s=+1、`cos_minus`はC=cos,s=-1、sinも同様です。
学生ノートのexp[-(t/tau_note)^2]ではtau_note=sqrt(2)*100に相当します。
外場は整数時刻の表を作り、RK4中間時刻では5点4次Lagrange補間を使います（field_mode=table4）。
中間時刻の解析式を直接評価した計算とは、有限dtで小さな差があり得ます。

## ハミルトニアンと観測量の規約

rho_k=(I+r_k dot sigma)/2、Phi=(mean(r_x)-i mean(r_y))/2。
xi_k=-2J(cos(kx)+cos(ky))-mu。

$$
h_k(t)=\begin{pmatrix}
\xi_{k-A(t)} & U\Phi(t)\\
U\Phi(t)^* & -\xi_{k+A(t)}
\end{pmatrix},\qquad i\dot\rho_k=[h_k,\rho_k].
$$

$$
J_x=\frac1N\sum_k 2J[\sin k_x\cos A_x-r_{k,z}\cos k_x\sin A_x],
$$

J_yもxをyに置き換えた同じ定義。
