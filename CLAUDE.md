# csp — Crystal Structure Prediction Package

大野論文（Ono et al., JACS submitted）の結晶構造予測手法を再現・パッケージ化し GitHub + Zenodo で公開するプロジェクト。

**全体仕様書**: `spec.md` を参照（論文対応表・GUI設計・ファイル構成など）

---

## ディレクトリ規則

- `ono_paper/` — このリポジトリ（ソースコードのみ）
- `ono_paper_dir/` — 計算結果・設定ファイル（ローカル・HPC共通）※ `auto_opt_dir/` と同じ規則

## 計算パイプライン概要

| Step | 内容 | 計算手法 |
|------|------|----------|
| 1a | 層内 vdW 粗探索 | vdW剛体球 → S=a×b 最小化 |
| 1b | 層内 DFT-D 精密化 | Gaussian16: B3LYP-D3/6-311G** (`gd3`, zero damping), BSSE |
| 2 | 長軸傾斜マップ | DFT-D → E_intra(6) ヒートマップ |
| 3 | 層間スタッキング | vdW V(x,y) マップ + DFT-D E_inter(7) |
| 4a/4b | ねじれ・不均一傾斜精密化 | DFT-D → E_int(near) 最小化 |
| 5 | Transfer integral | Gaussian16: B3LYP/6-31G*, 松井先生コード |

**使わないもの**: Amber力場、遠距離相互作用外挿（G1〜G5、`legacy/ono_scripts/tot_energy/`）、
粉末XRD比較（`legacy/ono_scripts/XRD_pattern/`、2026-07-06三好さん暫定決定）。
いずれもcsp本体には統合せず、README等での紹介に留める。

## 主要変数（論文の11変数）

- α: 分子面とグライド面のなす角（Step 1, 初期値 ≈25°）
- a, b: 層内格子定数
- θ_incl, φ_incl: 長軸傾斜角・方向（Step 2）
- θ_twist: ねじれ角（Step 4a, G形用）
- θ'_incl, φ'_incl: 不均一傾斜（Step 4b, N形用）
- x, y, z: 層間ズレ・間距離（Step 3）

## UI

Streamlit（英語）、**5タブ構成**（Tab 1: Molecule Setup / Tab 2: Step 1 Intralayer / Tab 3: Step 2（para・twistサブタブ）/ Tab 4: Step 3（para・twistサブタブ）/ Tab 5: Transfer Integrals）。
2026-07-05にpara/twistの扱いの非対称性を解消するため7タブから再編成（詳細は spec.md「大野コード対応表」）。

実行方式（2026-07-04 ミーティング決定、詳細は spec.md「実行方式」）:
- vdW スキャン（Step 1a）のみ GUI 内で実行。DFT 系ステップは CLI 実行とし、GUI はコマンド紹介と結果 CSV の drag & drop 表示のみ
- 結果表示は「ユーザー CSV 未ドロップ時は example/pentacene/ のサンプルを表示、ドロップで切替」方式
- サンプル（論文図のデータ）同梱は**後回し**（2026-07-05 決定）。残タスク一覧は spec.md「残タスク / TODO」を参照
- 各タブの変数・入出力・グラフの読み方は spec.md「タブ別 詳細ガイド」に記載（質問があるたびに追記していく方式）
- 大野さんの legacy コードに見つけた不具合（append()の引数ミス、pandas2.0でのDataFrame.append廃止、
  モノマーCSVパスのハードコード、tcal_csv.pyの無意味な固定パスプレフィックス等）は
  三好さんの許可のもと修正済み。詳細は spec.md「大野さんのコードに見つけた不具合の修正」参照

---

## auto_opt からの流用コード

リポジトリ: `/Users/miyoshimao/Working/auto_opt/`

| 流用元 | 使用用途 | コピー先 |
|--------|---------|---------|
| `src/auto_opt/utils.py` | `place_monomer()`, `vdw_radius()`, `R2atom()`, `Rod()`, `read_xyz()` | `src/csp/structure/` に組み込み |
| `src/auto_opt/cluster.py` | SGE ジョブ投入（qsub/qstat） | `src/csp/dft/job_cluster.py` |
| `src/auto_opt/app.py` | Streamlit UI パターン（HPC設定UI、SCPコマンド表示） | `app.py` の参考 |
| `src/auto_opt/plot/make_cluster_xyz.py` | py3Dmol 3D表示ロジック | `src/csp/plot/viewer3d.py` |
| `src/auto_opt/plot/energy_map.py` | Plotly ヒートマップ | `src/csp/plot/map2d.py` |

### place_monomer() の注意点

`auto_opt/utils.py` の `place_monomer()` はグライド対称・スクリュー対称の配置に対応している。
csp では「グライド対称のみ」を使うので引数の使い方に注意:
- `phi` → α（ヘリンボーン角の半分）
- `alpha` → θ_incl（長軸傾斜）
- `monomer_dir` → 分子 XYZ/CSV のディレクトリ

### vdW スキャンロジックの参考

`src/auto_opt/vdw/sweep_phi.py` に vdW 接触スキャンのロジックがある。
ただし csp では変数体系が異なる（αスキャン vs φスキャン）ので、参考程度に。

---

## HPC 設定

- ホスト: `miyoshi@133.11.68.31`（スパコン）
- SGE キューシステム（qsub/qstat/qacct）
- Gaussian16 環境: HPC 上にある（パスは HPC の `.bashrc` 参照）
- 作業ディレクトリ: HPC 上の `~/Working/ono_paper_dir/`

## 分子データ

`data/molecules/` にポリアセン XYZ（B3LYP-D3/6-311G** 最適化済み）:
- `naphthalene.xyz`, `anthracene.xyz`, `tetracene.xyz`, `pentacene.xyz`, `hexacene.xyz`
- 大野さん・小野さんから取得（まだない場合は Gaussian で最適化が必要）

## 公開予定

- GitHub: `mao-tf/ono-csp`（作成済み・公開中: https://github.com/mao-tf/ono-csp）
- コラボレーター: 大野さん `ryota-ohno`（承認済み・push権限あり）
- 大野さんの既存スクリプト**受領・読解・統合済み**（2026-07-04, `legacy/ono_scripts/`）:
  - `stepwise_optimization/` — Step 1〜3 の山登り法最適化。幾何は `src/csp/structure/intralayer.py`、
    vdWスキャンは `src/csp/vdw/contact.py`（Tab 2 でGUI実行可）、gjf生成 `dft/make_gjf.py`、
    counterpoiseログ解析 `dft/parse_log.py` に移植済み
  - `tcal_csv/` — transfer integral 一括計算（松井研 tcal ベース、https://github.com/matsui-lab-yamagata/tcal）
  - **変数対応・CSVフォーマット・要確認事項は spec.md「大野コード対応表」を正とする**
    （theta=ヘリンボーン半角α・z軸回転、A2=ねじれ、長軸=z、E_intra(8)=4Et+2Ep1+2Ep2）
- Zenodo で DOI 取得（GitHub Release から自動生成）
- 論文引用: Methods section に URL + DOI を記載
