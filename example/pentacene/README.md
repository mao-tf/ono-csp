# example/pentacene/

大野さんの事前計算結果（論文再現用サンプル）の置き場所です。
GUI の各タブは、ユーザーが自分の CSV をドロップするまでここのファイルを
サンプルとして表示します（spec.md「実行方式」参照）。

ファイルはまだ未受領です。受領したら以下の名前で置いてください
（列名が異なる場合は GUI 側の列選択デフォルトを合わせて調整します）:

| ファイル | 内容 | 対応タブ |
|---------|------|---------|
| `step1_vdw.csv` | α vs S=a×b（vdW 粗探索） | Tab 2 |
| `step1_dft.csv` | α vs E_intra(8)（DFT-D 精密化） | Tab 3 |
| `step2_map.csv` | θ_incl, φ_incl, E_intra(6) マップ | Tab 4 |
| `step3_N1.csv` | x, y, V, E_inter(7)（N1 スタッキング） | Tab 5 |
| `step3_N2.csv` | 同上（N2 スタッキング） | Tab 5 |
| `step4a_twist.csv` | θ_twist vs E_int(near) | Tab 6a |
| `step4b_incl.csv` | θ'_incl vs E_int(near) | Tab 6b |
| `transfer_integrals.csv` | 多形別・接触タイプ別 J | Tab 7 |
