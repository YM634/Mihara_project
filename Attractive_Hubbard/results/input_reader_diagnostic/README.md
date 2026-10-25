# 入力読込問題の診断記録

このフォルダは、極小Gaussian裾の非正規化数でC++の入力読込が途中停止した際の診断データ。正しい外場付き計算結果として使用しない。

cpp_3/Ax_t.txtの315番目（0始まり）の値は5e−324。reference/check_field_reader.exで調べるとread_count=315、fail=1、eof=0となり、パルス本体に到達していない。元コードのfind_param_tvectorは各値の読込失敗を確認していない。

実際の比較では、両実装の共通外場からabs(value)<1e−300をゼロにして回避した。正しい比較結果は ../cpp_comparison/ にある。
