# example/pentacene/

大野さんの事前計算結果（論文再現用サンプル）の置き場所です。
GUI の各タブは、ユーザーが自分の CSV をドロップするまでここのファイルを
サンプルとして表示します（spec.md「実行方式」参照）。

現在の5タブ構成（Tab1 Molecule Setup / Tab2 Step1 Intralayer / Tab3 Step2
（para・twistサブタブ）/ Tab4 Step3（para・twistサブタブ）/ Tab5 Transfer
Integrals）に対応させると、以下の名前・列で置いてください:

| ファイル | 列 | 対応タブ | 状態 |
|---------|-----|---------|------|
| `step1_init_params.csv` | a, b, theta, S, kind, status | Tab 2（vdW 初期候補） | **生成済み（2026-07-06。大野さんから受領していないが、`step1a_scan()`（`step1.py`の`init_process()`と同じvdW粗探索）をペンタセンでローカル実行するだけで再現可能。DFT不要・数秒で終わる。受領済みのstep1.csvと突合し22候補全て一致することを検証済み）** |
| `step1_vdw_curves.csv` | alpha, theta_ab, a, b, S, valid | Tab 2（Fig.S1(c)風 β スイープ全曲線） | **生成済み（2026-07-06。`step1_init_params.csv`と同じ`step1a_scan()`の副産物。これが無いとFig.S1(c)サブタブがスキャン未実行時に表示されないため追加）** |
| `step1.csv` | a, b, theta, E, E_p1, E_p2, E_t, status, file_name | Tab 2（E_intra(8) 山登り履歴） | **受領済み（2026-07-06、ペンタセン実データ、2692行）** |
| `step2_para.csv` | z, Et, Ep | Tab 3 para（長軸シフトスキャン） | **受領済み（2026-07-06、ペンタセン実データ、81行）** |
| `step2_twist.csv` | a, b, theta, Rt, A2, E, E_p1, E_t, … | Tab 3 twist | 未受領 |
| `step3_para.csv` | cx, cy, cz, a, b, theta, Rt, Rp, E, …（per-dimer E） | Tab 4 para（層間最適化） | 未受領 |
| `step3_twist.csv` | cx, cy, cz, …, E | Tab 4 twist | 未受領 |
| `transfer_integrals.csv` | 任意（例: a, b, theta, J_t, J_p または polymorph, contact, J） | Tab 5 | 未受領 |

tcal の `result.txt`（スペース区切り）は CSV に変換してから置いてください。
列名が多少違っても GUI 側の列選択で対応できますが、上記に合わせるとデフォルトで正しく表示されます。
