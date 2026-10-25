# x方向2D応答用の時系列スキャン

2個の独立Juliaプロセス、各1スレッド。平均場更新あり、Float64。
設定の正本はscan_config.toml。omega_000..100はomega1=0..5（0.05刻み）。負側は実時系列の共役対称性を使って後処理する。omega2の変換は後処理で行う。

外場はベクトルポテンシャルAx=f+epsilon*g、Ay=0。物理電場はE=-dA/dt。学生ノートの包絡exp[-(t/(sqrt(2)/0.01))^2]に対応させ、コードのtau_g=100を採用した。

各ディレクトリのobservables.f64はヘッダなしlittle-endian Float64、140001行10列のJulia列優先配列。Pythonではnp.fromfile(path,dtype="<f8").reshape((140001,10),order="F")で読み込む。列はtime,n_up,n_down,Phi_re,Phi_im,Jx,Jy,Ekin,Eint,Etot。initial_state.f64とfinal_state.f64は同形式の3行2500列。

complete.tomlが存在する条件だけ保存完了とみなす。config.tomlに個別条件を保存。全810条件で観測量9,072,064,800バイト、初期・最終状態97,200,000バイト、合計約9.17 GB。ログ・設定を加えても約10 GB以内。

worker1.logとworker2.logで進捗を確認できる。全プロセスが停止した後、scripts/launch_scan.pyを再実行すると完了済み条件を飛ばして再開する。launch.jsonにPIDを保存。負周波数復元、chi抽出、フーリエ変換と図の作成は今後の後処理。
