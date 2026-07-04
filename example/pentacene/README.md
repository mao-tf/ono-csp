# example/pentacene/

大野さんの事前計算結果（論文再現用サンプル）の置き場所です。
GUI の各タブは、ユーザーが自分の CSV をドロップするまでここのファイルを
サンプルとして表示します（spec.md「実行方式」参照）。

ファイルはまだ未受領です。legacy スクリプトの実際の出力フォーマット
（spec.md「大野コード対応表」で確認済み）に合わせ、以下の名前・列で置いてください:

| ファイル | 列 | 対応タブ |
|---------|-----|---------|
| `step1_init_params.csv` | a, b, theta, status | Tab 2（vdW 初期候補） |
| `step1.csv` | a, b, theta, E, E_p1, E_p2, E_t, status, file_name | Tab 3（E_intra(8) 山登り履歴） |
| `step2_para.csv` | z, Et, Ep | Tab 4（長軸シフトスキャン） |
| `step3_para.csv` | cx, cy, cz, a, b, theta, Rt, Rp, E, …（per-dimer E） | Tab 5（層間最適化） |
| `step2_twist.csv` | a, b, theta, Rt, A2, E, E_p1, E_t, … | Tab 6a（twist 変種） |
| `step3_twist.csv` | cx, cy, cz, …, E | Tab 6b |
| `transfer_integrals.csv` | 任意（例: a, b, theta, J_t, J_p または polymorph, contact, J） | Tab 7 |

tcal の `result.txt`（スペース区切り）は CSV に変換してから置いてください。
列名が多少違っても GUI 側の列選択で対応できますが、上記に合わせるとデフォルトで正しく表示されます。
