# csp パッケージ 仕様書 / Package Specification

> 対応論文: Ono et al., "Origin of Layered Herringbone Packing and Polymorphism in Polyacenes: A Quantum Chemical Optimization Approach", JACS (submitted)

---

## 目的・背景 / Purpose

本パッケージは大野論文に記述されたポリアセン結晶構造予測の手法を再現し、GitHub + Zenodo で公開するためのものです。`auto_opt`（VdW→Amber→DFT探索）とは別の独立したパッケージです。

**設計方針**:
- Amber 力場は使わない（vdW剛体球粗探索 → DFT-D 直接）
- 遠距離相互作用の外挿（論文 §2.6）は除外
- Transfer integral は必須
- UI は英語（国際公開向け）

---

## パッケージ概要 / Package Overview

| 項目 | 内容 |
|------|------|
| パッケージ名 | `csp` (Crystal Structure Prediction) |
| ディレクトリ | `/Users/miyoshimao/Working/ono_paper/` |
| リポジトリ名 | `MaoMiyoshi/ono-csp`（仮） |
| 公開先 | GitHub + Zenodo (DOI自動生成) |
| ライセンス | MIT |
| Python | 3.10+ |
| UI | Streamlit（英語） |

**DFT計算設定（論文 METHOD section より）**:
- エネルギー計算: `B3LYP empiricaldispersion=gd3 6-311g** counterpoise=2` (Gaussian16)
- Transfer integral: `B3LYP/6-31g*` で HOMO 取得 → HOMO-HOMO オーバーラップ積分
- BSSE補正: Counterpoise法
- 分散補正は **D3 (zero damping) = `gd3`** で**確定**（2026-07-04確認）。
  大野さんの実コード `legacy/ono_scripts/stepwise_optimization/make_step*.py` の
  ルート行が `#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2`。
  論文 METHOD の引用文献69（Grimme JCP 2010 = zero damping版）とも一致
- **要確認**: transfer integral の MO 計算レベル。論文 METHOD は「B3LYP/6-31G*」だが、
  大野さんの `legacy/ono_scripts/tcal_csv/tcal_csv.py` は `#pbepbe/6-311G**` を使用。
  どちらが論文の図に使われた設定か大野さんに確認する

**対象分子（プリセット）**:
- naphthalene (n=2)、anthracene (n=3)、tetracene (n=4)、pentacene (n=5)、hexacene (n=6)
- カスタム分子: ユーザーが XYZ ファイルをアップロード

---

## 計算パイプライン / Calculation Pipeline

> **注意**: 以下の Step 1a〜5 は論文本文（§2.1〜2.8）から起こした初期ドラフト。
> 大野さんの実コード受領後（2026-07-04）、実際のパイプラインは下の
> 「大野コード対応表」の3ステップ構成であることが判明した。実装・GUI は
> 対応表を正とする。

論文の §2.1〜2.8 に対応する 5ステップ構成。全ステップで同じ 11 変数（α, a, b, θ_incl, φ_incl, θ_twist, θ'_incl, φ'_incl, x, y, z）を段階的に導入する。

### Step 1a: 層内 vdW 粗探索（論文 §2.1 前半）

**目的**: 各αに対してvdW接触が成り立つ最小のユニットセル面積 S=a×b を求める。

**変数**: α (5°刻み, 5°〜85°), a, b

**アルゴリズム**:
1. グライド対称を仮定してモノマー対を配置
2. 隣接8分子（T型4つ、SP型4つ）の原子間距離を計算
3. 全原子ペアについて `d_ij > R_i + R_j` (vdW和) を満たす最小の a,b を数値最適化
4. S=a×b を記録

**出力**: α vs S プロット、各αの初期構造 XYZ

---

### Step 1b: 層内 DFT-D 精密化（論文 §2.1 後半）

**目的**: DFT-D計算でE_intra(8)を最小化し、R形を確定する。

**変数**: α（固定、各点ごと）, a, b（最適化変数）

**計算式**:
```
E_intra(8) = 4 × E(T型接触) + 4 × E(SP型接触)
E(接触) = E(ダイマー) - 2×E(モノマー)  [BSSE補正済み]
```

**Gaussian16 入力テンプレート**:
```
# B3LYP empiricaldispersion=gd3 6-311g** counterpoise=2 nosymm

{分子名} dimer step1 alpha={alpha:.1f} a={a:.3f} b={b:.3f}

0 1
{モノマー1座標}
--
0 1
{モノマー2座標}
```

**HPC ジョブ**: SGE (qsub)、auto_opt の `cluster.py` を流用

**出力**: α vs E_intra(8) プロット、最安定α（≈25°）での R形構造

---

### Step 2: 長軸傾斜マップ（論文 §2.2）

**目的**: R形から長軸を傾けたときの E_intra(6) 2Dマップを作成し、G形・N形を特定する。

**変数**: θ_incl (0〜40°, 1°刻み), φ_incl (0〜360°, 5°刻み)

**計算式**:
```
E_intra(6) = 2×E(T型接触) + 2×E(SP型接触) + 2×E(第2SP型接触)
             [= E_intra(8) から第2T型2本を除いた6分子]
```

**グライド対称の扱い**:
- φ_incl = 0°, 180° → グライド対称保持 → **G形**
- φ_incl ≈ 48° (±), 132° (±) → グライド対称破れ → **N形**

**分子配置の生成**:
```python
# θ_incl, φ_incl による長軸傾斜
dz_T = (a*b / sqrt(a²+b²)) * tan(θ_incl) * cos(φ_incl - φ_T)
# φ_T: T型接触方向の角度
```

**出力**: θ_incl×φ_incl ヒートマップ、G形・N形の θ_incl/φ_incl 値

---

### Step 3: 層間スタッキング（論文 §2.3）

**目的**: 上下層のズレ (x,y) に対してvdW体積 V(x,y) と DFT-D エネルギー E_inter(7) をマップ化する。

**変数**: x (-a/2〜a/2, 0.1Å刻み), y (-b/2〜b/2, 0.1Å刻み), z (vdW接触距離で自動計算)

**Einter(7) の計算**:
```
E_inter(7) = E(層間最近接1分子) + 6×E(層間→層内隣接6分子)
           = E_inter(1,0) + Σ E_inter(1,j)  [j=1..6]
```

**N形の分岐**:
- N1: (x,y) ≈ (0.0Å, +1.4Å) → **Type II** 候補（長ポリアセン: n≥5）
- N2: (x,y) ≈ (+1.4Å, -0.7Å) → **Type IV** 候補（中ポリアセン: n=4）

**出力**: V(x,y) ヒートマップ、E_inter(7) ヒートマップ（R形/G形/N形各3枚）

---

### Step 4a: ねじれ精密化・G形用（論文 §2.4）

**目的**: G形にねじれ θ_twist を導入し E_int(near) を最小化する（Type III: n=2,3 向け）。

**変数**: θ_twist (0〜20°, 1°刻み), a, b, ΔzT, x, y, z（全最適化）

**Eint(near) の計算式**:
```
E_int(near) = E_intra(6) + 2 × E_inter(7)
```

**分子配置**: グライド対称を保ちつつ隣接T型分子を逆方向にねじる

**期待結果**:
- naphthalene: θ_twist ≈ 13° で最安定
- anthracene: θ_twist ≈ 9° で最安定
- tetracene以上: θ_twist 導入でエネルギー低下なし

**出力**: θ_twist vs E_int(near) プロット

---

### Step 4b: 不均一傾斜精密化・N形用（論文 §2.5）

**目的**: N2型スタッキングにユニットセル内一方の分子だけ独立した傾きを導入し Type II / Type IV の安定性を比較する。

**変数**: θ'_incl (0〜5°, 0.5°刻み), φ'_incl (固定: Step4a の E_intra(6) 最小値の方向)

**計算式**: E_int(near) = E_intra(6) + 2 × E_inter(7)（θ'_incl, φ'_incl を変数として最適化）

**期待結果**:
- tetracene: θ'_incl ≈ 2.5° で Type IV が Type II より安定
- pentacene/hexacene: θ'_incl 導入でも Type II が最安定

**出力**: θ'_incl vs E_int(near) プロット（Type II vs Type IV 比較）

---

### Step 5: Transfer Integral（論文 §2.8）

**目的**: 各多形（Type I〜IV）の最安定構造について T型・SP型接触の Transfer integral J を計算する。

**計算手法**: 松井先生提供コードを使用
- Gaussian16: `B3LYP/6-31g*` で HOMO 計算 → fchk/chk 出力
- J = ⟨ψ_HOMO^A | Ĥ | ψ_HOMO^B⟩（HOMO-HOMO オーバーラップ積分）

**計算対象**:
1. 各多形の最安定構造での J（T型・SP型）
2. α依存性: R形で α を変化させたときの J(α)
3. θ_incl依存性: N形で θ_incl を変化させたときの J(θ_incl)

**出力**: 多形別 J 棒グラフ、J(α) プロット、J(θ_incl) プロット

---

## 実行方式 / Execution Policy（2026-07-04 ミーティング決定）

- **Step 1a（vdW 粗探索）**: ローカル PC で完結するので **GUI 内で直接実行**する
- **DFT 系ステップ（Step 1b, 2, 3, 4a/4b, 5）**: 計算機接続・計算時間の問題から **GUI からは実行しない**。
  GUI には「このコード（CLI）で実行できます」という紹介と実行コマンドを表示するのみ
- **結果表示**: 各タブに結果 CSV の **drag & drop 欄**を置く。
  未ドロップ時は `example/pentacene/` の大野さんによる事前計算サンプルを表示。
  ユーザーがファイルをドロップしたらそちらのグラフ・3D モデル表示に切り替わる
  ※ **サンプル（論文で使った図のデータ）の同梱は後回しにする（2026-07-05 決定）**。
  GUI はサンプル未配置でも動作し「Sample results are not bundled yet」表示になる。
  受領したら `example/pentacene/README.md` の表に従って配置するだけでよい
- 機能は最小限に：「基本の流れが分かりやすくパッケージ化されている」ことを優先（ミーティングでの先生方針）

以降のワイヤーフレーム中の `[Run DFT Jobs]` 等の DFT 実行ボタンはこの方針により
「CLI コマンド紹介 + 結果 drag & drop」に置き換えて実装する。

---

## 大野コード対応表 / Legacy Code Mapping（2026-07-04 受領・読解・統合）

大野さんの実パイプラインは3ステップ構成（step1 → step2 para/twist → step3
para/twist → tcal）で、論文§の5ステップドラフトとは粒度が異なる。
GUI タブは受領コードに合わせて再構成済み（旧ワイヤーフレーム Tab 4〜6 は置換）。

**タブ構成をさらに5タブへ再編成（2026-07-05, 三好さんの指摘を受けて）**:
7タブのままだとpara/twistの扱いが非対称だった（paraだけTab4/5で独立、twistだけ
Tab6にstep2/step3両方まとめられていた）。パイプラインの粒度（3ステップ）に
一致させ、**ステップ番号でタブを分け、para/twistはステップごとのサブタブ**にした：

| Tab | 内容 |
|---|---|
| 1. Molecule Setup | 変更なし |
| 2. Step 1 – Intralayer | vdW粗探索（GUI実行）→DFT精密化（CLI）の2セクション（順番に行う工程なのでサブタブでなく区切り） |
| 3. Step 2 – Long-Axis Shift | サブタブ: para / twist（どちらかを選ぶ工程なのでサブタブ） |
| 4. Step 3 – Interlayer Stacking | サブタブ: para / twist |
| 5. Transfer Integrals | 変更なし |

| legacy (ono_scripts/) | 内容 | csp 統合先 / GUI |
|---|---|---|
| stepwise_optimization/vdw.py + step1.py --init | vdW接触モデルで初期(a,b,θ)候補（T接触方向0〜90°掃引、Sの極小＋両端点） | `src/csp/vdw/contact.py step1a_scan` → **Tab 2 vdWセクション（GUI内実行・実装済み）**。CLI版: `csp/vdw/sweep.py` + `extract_init.py` + `csp/plot/vdw_scan.py`（下記） |
| make_step1.py + step1.py | (a,b) 0.1Å刻み3×3山登りDFT。E_intra(8)=4Et+2Ep1+2Ep2 | gjf生成: `src/csp/dft/make_gjf.py`、ログ解析: `parse_log.py` → **Tab 2 DFTセクション** |
| make_step2_para.py + step2_para.py | 長軸シフト z スキャン（0〜4Å×0.1、--Link1で41点/1ファイル）。出力 z,Et,Ep | **Tab 3 paraサブタブ** |
| make_step2_twist.py + step2_twist.py | Rt（T接触長軸シフト）と A2（ねじれ）導入後 (a,b) 再最適化 | **Tab 3 twistサブタブ** |
| make_step3_para.py + step3_para.py + step3_para_vdw.py | 層間cベクトル(cx,cy,cz) 3×3×3山登り（10ダイマー/点、2パターン平均）; vdW層間距離マップ z(x,y) | **Tab 4 paraサブタブ** |
| make_step3_twist.py + step3_twist.py | twist形（Type III）の層間最適化 | **Tab 4 twistサブタブ** |
| tcal_csv/ | transfer integral 一括計算（松井研 tcal ラッパー） | **Tab 5** |

**変数対応**:
- `theta` (= A3) = ヘリンボーン二面角の半分（論文のα、5〜45°）。z軸まわり回転。T型副格子は −A3
- `A2` = ねじれ（−x軸まわり回転。step1 では 0）
- `Rt`, `Rp` = T型 / SP型接触での長軸方向シフト
- `cx, cy, cz` = 層間 c ベクトル
- 分子座標系: **長軸 = z**、格子 a = x, b = y。T型隣接 (±a/2, ±b/2)、SP型 (±a,0),(0,±b)
- モノマー CSV: X,Y,Z,R 列（R = Bondi半径、元素は R2atom で逆引き）。**Tab 1 からダウンロード可**（`intralayer.monomer_csv`）

**counterpoise ログ解析**（`src/csp/dft/parse_log.py`）:
「SCF Done: E(R…」が 5 本 1 組（超分子 / ghost付き単量体×2 / 単体基底単量体×2）。
E_int = (E[0] − E[1] − E[2]) × 627.510 kcal/mol。--Link1 連結で複数ペア/1ログ。

**結果CSVフォーマット（GUI の drag & drop 対応列）**:
- `step1_init_params.csv`: a, b, theta, status
- `step1.csv`: a, b, theta, E, E_p1, E_p2, E_t, status, file_name
- `step2_para.csv`: z, Et, Ep
- `step2_twist.csv`: a, b, theta, Rt, A2, E, E_p1, E_t, status, file_name
- `step3_para.csv`: cx, cy, cz, a, b, theta, Rt, Rp, E, E_i01…E_ip4, status, file_name
- tcal: result.txt（スペース区切り J_t J_p）→ CSV 変換して Tab 7 へ

**汎函数の扱い（2026-07-05 三好さん決定: 論文に合わせる）**:
1. ~~tcal の MO 計算レベル不一致~~ → **解決**: `tcal_csv.py` のルート行を
   `pbepbe/6-311G**` から論文 METHOD の `b3lyp/6-31g*` に修正済み
2. ~~make_step2_twist.py 冒頭コメント「pbepbe+d3bjで計算」~~ → **解決**: 誤解を招く
   コメントを削除（ルート行は元から b3lyp+GD3 で論文と一致）
   ※ 大野さんが pbepbe を使っていた経緯は参考までに聞いておくとよい

**要確認3: モノマー CSV の面内絶対配向 → 2026-07-06 解決済み**:
大野さんから実モノマーCSV（`legacy/ono_scripts/monomer/*.csv`, 5分子）を受領し、
`intralayer.monomer_csv()`の出力と数値比較した。結果:
- X軸（面法線）: 両方0 → **完全一致**
- Y軸（短軸）: 符号含め完全一致
- Z軸（長軸）: **符号が反転**（例: anthracene で Ono側 Z=-1.221 vs csp側 Z=+1.221）
ポリアセンは長軸中心に対して対称（D2h）なので、Z符号反転は分子形状を変えない
（原子ラベルが入れ替わるだけ）→ **vdW/DFTエネルギー計算に一切影響なし**。対応不要。

**確定事項**: 大野さんのローマ字表記は **Ohno**（2026-07-05 確認、LICENSE 等に使用）。
サンプル結果 CSV（example/pentacene/ 用）は未受領。

**`legacy/ono_scripts/tot_energy/`（2026-07-05 大野さんアップロード、2026-07-06 解読）**:
`tot_energy.py` + `L1_norm.py` + `utils.py`。論文 **§2.6「Comparison with Experimental
Structures Considering Distant Intermolecular Interactions」の実装**と判明:
```
Eint(G1) = Eintra(8) + 2×Einter(9)
Eint(all) = Eint(G1) + Eint(G2) + Eint(G3) + Eint(G4) + Eint(G5) + ...（外挿）
```
`L1_norm.py`の`L1norm(N)`/`L1norm2(N)`が格子座標上でN番目の近接殻(shell)の分子リストを
生成し、`tot_energy.py`が殻を広げながら（I=1,2,3）intra/inter両方のクラスターxyzを作って
DFT計算、遠方分子まで含めた収束挙動を評価する。

**→ CLAUDE.mdの当初方針「使わないもの: 遠距離相互作用外挿（G1〜G5）」に該当する機能**。
2026-07-06 三好さんに統合するか確認したが応答待ち。**いったんスコープ外のまま
`legacy/`に保持し、csp本体（`src/csp`, GUI）には統合しない**（当初方針を踏襲、
いつでも変更可）。

**Step 1a の CLI 版パイプライン（2026-07-05 三好さん実装 + プロット追加）**:
```bash
PYTHONPATH=src python -m csp.vdw.sweep \
    --monomer data/molecules/pentacene.xyz --out-dir runs/pentacene/ \
    --alpha-min 5 --alpha-max 85 --alpha-step 1 --theta-step 1
# → runs/pentacene/vdW_r_contact_pentacene.csv（alpha, theta_c, R_clps, TorF）

PYTHONPATH=src python -m csp.vdw.extract_init \
    --vdw-csv runs/pentacene/vdW_r_contact_pentacene.csv \
    --out runs/pentacene/step1_init_params.csv --minima
# → 候補CSV（alpha, theta_c, a, b, S, structure_type=a-stack/b-stack/local_min, status）

PYTHONPATH=src python -m csp.plot.vdw_scan \
    --init-csv runs/pentacene/step1_init_params.csv \
    --out runs/pentacene/vdw_scan_pentacene.png
# → S=a×b vs alpha プロット（論文 Fig.2(b) の vdW 版）
```
- 前提: モノマー XYZ は大野座標系（長軸=z、面法線=x）— data/molecules/ のファイルはそのまま使える
- pentacene 実行結果（α 1°刻み）: 全体最小 α=21°（S=42.3 Å²、b-stack）、
  **local_min は α≥23° で出現**（論文の議論と整合）。α↔90−α で a/b 入替の鏡映対称を確認
  （a-stack 最小 α=69°=90−21）→ 面内配向規約の懸念（要確認3）はスキャン全域では実質無害
- `runs/` は .gitignore 済み（結果はリポジトリに含めない）
- 注意: GUI Tab 2 の `contact.step1a_scan` と CLI 版は同じアルゴリズムの独立実装。
  将来どちらかに一本化するのが望ましい（TODO 参照）

**Tab 2 の2グラフ構成（2026-07-05 追加、論文 Fig. 2(b) / Fig. S1(c) のvdW版）**:
- `step1a_scan` の戻り値 `df_curves` に `valid` 列（bool）を追加。β（=theta_ab）掃引の
  全点を保持し、SP隣接のvdW接触を満たす点=True／破綻する点=Falseとして記録
  （論文 Fig. S1(c) の実線／破線に対応）
- **Fig. 2(b)-style タブ**: x=α, y=S。feasible領域のmin envelope（グレー線）＋
  候補点（low-β endpoint 青／high-β endpoint 赤／interior local min 緑）。
  点をクリックすると該当構造の9分子クラスターをプレビュー
- **Fig. S1(c)-style タブ**: x=β, y=S。複数αを選択（デフォルトはスキャン範囲を
  5等分した値 — α範囲5〜45°では偶然 [5,15,25,35,45] と論文の例示値に一致）。
  各αを1色とし、有効区間=実線／無効区間=破線で描画。点クリックでプレビュー
- クリック選択はplotlyの`on_select="rerun"`＋`customdata`（α,β,a,bを埋め込み）で実装。
  候補リストのセレクトボックスも含め3種の入力（Fig.2b点／Fig.S1c点／リスト選択）を
  「直前に変化したものだけを反映」する差分検知方式で統合（`st.session_state["s1vdw_current"]`）
- 制約: Fig. S1(c)-styleタブは `df_curves`（GUI内スキャン直後のみ）が必要。
  アップロードした`step1_init_params.csv`だけでは全掃引点がないため非表示

**レイアウト調整（2026-07-05 三好さんフィードバック反映）**:
- Fig. 2(b)-style: 同じ候補カテゴリ（low-β/high-β endpoint、interior local min）内で
  α順にソートして線で接続（論文のHB/PS/CH別々の連続曲線に対応する見せ方）
- Fig. S1(c)-style: y軸（S）の表示範囲を40からに固定（下側をカット、データ最大値まで自動）
- 各グラフを2カラム化（グラフ:3Dプレビュー = 2:1）し、クリックで隣に9分子構造が
  即座に表示されるように変更。3DビューアはFig.2b用・Fig.S1c用で別インスタンス
  （スタイル切替は独立、内部状態`s1vdw_current`は共有）
- **未検証**: plotlyグラフの実クリック操作はAppTest（自動テスト基盤）でシミュレートできないため、
  ドロップダウン・複数選択等のロジックはテスト済みだが、実ブラウザでのクリック動作は
  三好さんご自身での確認をお願いしたい

**バグ修正・凡例改善（2026-07-05 三好さんフィードバック反映 その2）**:
- 凡例名を大野コードの用語に統一: `low-β endpoint`→**b-stack**、`high-β endpoint`→**a-stack**、
  `interior local min of S`→**local min**（extract_init.py の structure_type と同じ語彙）
- **不具合修正 v1（不十分だった）**: Y軸範囲をスキャン時に一度だけキャッシュし
  `uirevision`固定文字列で安定化を試みたが、三好さんに実ブラウザで確認いただいたところ
  **0と40を交互に行き来する現象は解消しなかった**
- **不具合修正 v2（2026-07-05 再修正）**: `figB.update_yaxes(..., fixedrange=True)` を追加。
  `autorange=False`とキャッシュ済み固定値の指定に加えて、Y軸自体をplotly側で
  ズーム・パン不可（インタラクション無効）にすることで、クリック起因の再レンダリングが
  どんな経路であれY軸を動かせないようにする、より強い対策。
  自動テスト（AppTest）ではリランを跨いで値が不変であることを確認したが、
  **plotlyグラフの実クリックをヘッドレスブラウザ（Playwright）で再現しようとしたところ、
  Streamlitのタブ内レイアウト特有の座標計算の問題でクリック自体をうまく合成できず、
  「クリックしてY軸が動くか」という肝心の部分は自動検証しきれていない**。
  → **2026-07-05 三好さん実ブラウザ確認: 解消済み**

**視認性改善（2026-07-05）**: Fig. S1(c)-style の実線／破線が見分けにくいとのフィードバックを受け、
破線側を `dash="longdash"`（長め）・線幅太め(4 vs 実線2.5)・不透明度0.55（薄く）に変更。
実線側は線幅2.5・不透明度1.0のまま。マーカーサイズも実線5px/破線3pxで差をつけた

---

## タブ別 詳細ガイド / Per-Tab Deep-Dive（三好さん向けメモ、質問があるたびに追記）

以下は「GUI タブ構成」の元ワイヤーフレーム（論文ドラフトベース、実コード受領前）とは別に、
**実際に実装されている挙動**を tab ごとに詳しく書いたメモ。まだ全タブ分は書いていない
（聞かれた分だけ随時追記）。

### CLI系タブ（Tab 3, 4, 5, 6a, 6b, 7）の統一フォーマット（2026-07-05）

「GUIは大野さんのパッケージの使い方を説明するもの」という先生方針に立ち返り、
全CLI系タブの「How to run」を**4段構成**に統一（`app.py` の `cli_howto()` ヘルパー）：

1. **What** — このステップが何をするか、大野さんのどのスクリプトが対応するか
2. **Prepare inputs** — 事前に用意する入力ファイル名・列・値の出どころ（前のタブのどの結果を使うか）
3. **One-time setup** — `CSP_MONOMER_DIR` 環境変数、pjsub→qsub置き換えの必要性
4. **Run / Output** — 実行コマンドと出力ファイルの列

以前はタブによってこの情報の粒度がバラバラ（Tab3だけ詳しく、Tab6はコマンド一行だけ、
`step3_para.py`のRt/Rpがどこから来るか説明なし、等）だったのを解消した。

### 大野さんのコードに見つけた不具合の修正（2026-07-05、報告事項）

上記の統一作業でコマンドを実際に動かして確認する過程で、`legacy/ono_scripts/` 内に
以下の不具合を見つけて**修正済み**（三好さんの許可のもと）。修正後、実際に
`step1.py --init` のフルフロー（vdWスキャン→ジョブキュー登録の3イテレーション）と
`tcal_csv.py --init`（gjfファイル生成）を実行し、正常動作を確認した。

1. **`step2_para.py`**: `para_list.append(z,E1,E2)` は `list.append()` が引数を1つしか
   取らないため実行時エラーになる。`para_list.append((z,E1,E2))`（タプル化）に修正
2. **`step1.py`, `step2_twist.py`, `step3_para.py`, `step3_twist.py`**: いずれも
   `df_E = df_E.append(df_newline, ignore_index=True)` を使用しているが、
   `DataFrame.append()` は **pandas 2.0 で廃止**済み（本パッケージは `pandas>=2.0` を要求）。
   `df_E = pd.concat([df_E, pd.DataFrame([df_newline])], ignore_index=True)` に置き換え
3. **モノマーCSVパスのハードコード**: `make_step1.py`, `make_step2_para.py`,
   `make_step2_twist.py`, `make_step3_para.py`, `make_step3_twist.py`, `tcal_csv.py` の
   6ファイルすべてで `~/path/to/monomer/{monomer_name}.csv`（または`/path/to/tcal_csv/monomer/`）
   がハードコードされていて、使うたびにソースコードの手編集が必要だった。
   環境変数 `CSP_MONOMER_DIR`（未設定時は従来通り`~/path/to/monomer/`にフォールバック）で
   指定できるように変更。GUI側は「Tab 1 でダウンロードしたCSVを置いたディレクトリを
   `export CSP_MONOMER_DIR=...` で指定してください」という案内に統一
4. **`tcal_csv.py` の `/path/to/tcal_csv/{work_dir}` という無意味な固定プレフィックス**
   （25箇所）: `--auto-dir` で渡した実際のパスの前に、常にこの文字列が付加されてしまい、
   `--auto-dir` オプションが事実上機能していなかった（Ono さんのローカル環境の絶対パスが
   そのまま残っていたと思われる）。プレフィックスを削除し `--auto-dir` の値をそのまま使うよう修正。
   `job.sh`/`tcal_1.py` のコピー元は `work_dir` ではなく **スクリプト自身のディレクトリ**
   （`_SCRIPT_DIR`）を参照するよう修正（元コードは誤って同じ固定パス配下からコピーしようとしていた）
5. **`tcal_csv.py` の `subprocess.run(['cd',path])`**（3箇所）: `cd`はシェル組み込みコマンドで
   独立した実行ファイルではないため、Linux環境では`FileNotFoundError`でクラッシュする
   （直後にある実際に効く`os.chdir(path)`より前の、無意味な死んだコードだった）。削除

**変更していないもの**: `pjsub`（Fugaku/PJM向けジョブ投入）は環境依存性が高く、
我々のSGEクラスタ向けに安全に自動書き換えできる確証がないため、GUI側で
「qsubに置き換えてください」と案内するに留め、コード自体は変更していない。

### Tab 3: Step 2 – Long-Axis Shift（para / twist サブタブ、2026-07-06 Tab番号更新）

Tab 3 は para・twist いずれも**出発点は Tab 2 の DFT 結果**（最適化済み `a, b, theta`）。
para は「長軸方向にずらす」、twist は「ねじる＋格子を再最適化する」という、
Step 1 の配置に対する**別々の変形の入れ方**を試している、という位置づけ。

#### paraサブタブ

**やっていることの概要**:
Step 1（Tab 2）で決めた層内配置（格子定数 a, b、ヘリンボーン半角 α=theta）を**固定**した状態で、
隣接分子を「分子自身の長軸方向（z軸）」にずらしたとき、エネルギーがどう変わるかを見るステップ。
`legacy/ono_scripts/stepwise_optimization/step2_para.py`（+ `make_step2_para.py`）に対応。

**入力**:
- `step2_init_params.csv`: **1行だけ**のCSVで列は `a, b, theta`（Step 1 の最適値をそのまま使う）。
  GUIには生成ボタンはなく、Tab 2 の結果から手動で作って作業ディレクトリに置く
- 分子のモノマーCSV（Tab 1 からダウンロードしたもの）

**変数の意味**:
- `a, b, theta`: Step 1 で固定した層内パラメータ（thetaはヘリンボーン半角α）。Step 2 では**動かさない**
- `z`: このステップの主役。隣接分子を、その分子自身の長軸方向に平行移動させる量（Å）。
  0〜4Åの範囲を0.1Å刻みで41点計算し、対称性（+zと−zは物理的に等価という前提）を使って
  −4〜+4Åの81点に**ミラーリング**してCSVに保存する
- 計算対象は2種類の隣接分子ペア: **T型接触**（隣の分子が90°近い角度で接している配置）と
  **SP型（slipped-parallel）接触**（長軸がほぼ平行に少しずれて重なる配置）

**出力（`step2_para.csv`）の列**:
| 列 | 意味 |
|---|---|
| `z` | 長軸方向のシフト量（Å、-4〜+4） |
| `Et` | T型接触ペアの相互作用エネルギー（kcal/mol、counterpoise補正済み） |
| `Ep` | SP型接触ペアの相互作用エネルギー（kcal/mol） |

**グラフの見方（paraサブタブの Plotly図）**:
- 横軸: `z`（長軸方向のシフト量、Å）
- 縦軸: エネルギー（kcal/mol、負の方向が安定）
- 3本の線:
  - **E_t（実線+マーカー）**: T型接触のエネルギー変化
  - **E_p（実線+マーカー）**: SP型接触のエネルギー変化
  - **4·E_t + 2·E_p（点線）**: Step 1 の `E_intra(8) = 4Et + 2Ep1 + 2Ep2` の考え方をこの文脈に
    当てはめた「参考の合計値」。※ CSV自体にはこの列はなく、GUI側で計算して重ねて表示しているだけ
- **読み方**: 曲線の谷（極小）の位置を見る。谷が `z=0` にあれば「Step1で決めた配置がこの方向にも
  すでに最適」ということ。`z=0`からずれた位置に谷があれば、「長軸方向にずらした方がさらに安定になる」
  ことを意味し、次のステップ（ねじれ・不均一傾斜の精密化）の出発点になる

**論文の図との対応（2026-07-06 論文§2.2本文で確認・解決）**:
- 論文§2.2「Optimization by Uniform Long-Axis Inclination」に明記あり:
  R形の長軸傾斜は「層表面を平らに保ったまま**分子を長軸方向に一斉にスライドさせる**」ことで
  表現され（Fig.5a）、その結果をθ_incl（傾斜角）・φ_incl（層内での傾斜方向）の2変数に
  変換して図示している。つまり**θ_inclは物理量ではなく、実際の物理量である「長軸シフト」を
  角度に変換した表現**。`step2_para.py`の`z`はこの「長軸シフト」そのものであり、変換式も
  SI(Fig.S6a)に明記されている:
  ```
  N形: ΔzT = (ab/√(a²+b²))・tan(θ_inclN) ≈ 2.4 Å
  G形: ΔzT = (b/2)・tan(θ_inclG) ≈ 2.4 Å
  ```
  → `step2_para.py`のzスキャン範囲(0〜4Å)はこのΔzT≈2.4Åを問題なくカバーしており、
  **Tab 3 paraは論文§2.2をほぼそのまま実装したもの**と確認できた
- **φ_inclの扱いに注意**: 論文ではφ_incl=0°/90°でグライド対称保持（→G形）、中間値で
  対称性が破れる（→N形）。`step2_para.py`はT型・SP型接触に**同じz**を一律適用しているため、
  これは実質**φ_inclを0°/90°に固定した特殊ケース（G形方向）のみ**をスキャンしていることになる。
  N形（中間的なφ_incl）を探すにはT型・SP型で別々のシフト量が必要
- **確認済み（2026-07-06、三好さんの提案をコードで検証）**: `step3_para.py`の`Rt`・`Rp`が
  まさにこの「T型接触ごとに別々のシフト」を実現している。`make_step3_para.py`の実座標:
  ```
  it1 = (a/2,  b/2, Rt)       it2 = (a/2, -b/2, Rt-Rp)
  it3 = (-a/2,-b/2, -Rt)      it4 = (-a/2, b/2, -Rt+Rp)   ##it1/it2の反転対
  ip1 = (0, b, Rp)            ip2 = (0,-b, -Rp)            ##SP型（a>bの場合）
  ```
  Rt1=Rt, Rt2=Rt-Rp と置くと it1=(a/2,b/2,Rt1)、it2=(a/2,-b/2,Rt2)、
  ip1=(0,b,Rt1-Rt2) となり、**Rp = Rt1-Rt2（2つのT型接触のシフト差）** そのもの。
  Rp=0（Rt1=Rt2）→ 対称・G形方向、Rp≠0（Rt1≠Rt2）→ 対称性が破れる・**N形方向**。
  ただし `fixed_param_keys=['a','b','theta','Rt','Rp']` なので Rt,Rp は山登り法で
  自動最適化されず、**init CSVに複数の(Rt,Rp)行を用意して手動でスキャンする**設計
  （自動最適化されるのは cx,cy,cz のみ）

**コードを読んでいて見つけた気になる点**（2026-07-05 修正済み。詳細は下の「大野さんのコードに
見つけた不具合の修正」参照）:
- `step2_para.py` の `end_process()` 内、`para_list.append(z,E1,E2)` は
  Python の `list.append()` が引数を1つしか取れないため実行時エラーになるバグだった → 修正済み

#### twistサブタブ

**やっていることの概要**:
paraと同じくTab 2のDFT結果（`theta, a, b`）を出発点にするが、こちらは長軸方向のシフトではなく
**T型接触に「ねじれ」を導入**する。G形（グライド対称）から Type III 多形へ向かう精密化ステップ。
`legacy/ono_scripts/stepwise_optimization/step2_twist.py`（+ `make_step2_twist.py`）に対応。

**入力**:
- `step2_twist_init_params.csv`: 列は `theta, Rt, A2, a, b, status`。
  `theta` と初期 `a, b` は Tab 2 の DFT 結果から。`Rt`（T型接触の長軸方向シフト）と
  `A2`（ねじれ角、-x軸まわり回転）は新しく導入する変数で、試したい値のグリッドを自分で用意する
  （例: A2 = 0〜20°を数点）

**変数の意味**:
- `Rt`: T型接触分子を長軸方向にどれだけずらすか（Å）。paraの`z`と似ているがT型接触専用の変数
- `A2`: ねじれ角（度）。導入後、各`(theta, Rt, A2)`の組み合わせごとに`(a, b)`を
  0.1Å刻みの山登り法で**再最適化**する（paraは`a,b`固定、twistは`a,b`も動かす点が違う）
- `E = 4·E_t + 2·E_p1`（SP型は1方向のみ数える点がStep1のE_intra(8)と異なる）

**出力（`step2_twist.csv`）の列**: `a, b, theta, Rt, A2, E, E_p1, E_t, status, file_name`

**グラフの見方**: 横軸`A2`（ねじれ角）、縦軸`E`。極小の位置を探す。
論文の予測では naphthalene ≈13°、anthracene ≈9° で最安定、tetracene以上ではねじれを
入れてもエネルギーが下がらない（Type IIIが安定にならない）はず（論文§2.4）。

**Rt/Rpの全体像（2026-07-06、三好さんの質問を受けて整理）**:
Step2とStep3でRt（・Rp）の役割はまったく違う。混同注意:

| | 対象 | 自由度 | 座標での現れ方 |
|---|---|---|---|
| Step2 twist の `Rt` | 層内（1層のみ） | 1個（グライド対称を保つ） | T型接触: z=Rt、SP型(b方向): z=2Rt、SP型(a方向): z=0 |
| Step3 para の `Rt`,`Rp` | 層間（cx,cy,cz最適化の前提として固定入力） | 2個独立（Rp≠0で対称性が破れる） | it1=(a/2,b/2,Rt)、it2=(a/2,-b/2,Rt-Rp) |

つまり Step2 の twist は「G形→Type III」精密化（グライド対称を保ったまま）で、
**N形探索ではない**（以前の説明の通り）。RpがT型接触2つを独立に動かす自由度として
現れるのは Step3 だけ。ただし Step3 自身も Rt,Rp を「新しく探索する」わけではなく
（`fixed_param_keys`に含まれ、cx,cy,czだけが自動最適化される）、**外から与えられた
Rt,Rpの組に対して層間最適化を行うだけ**。

**気になる点（新規、大野さんに確認したい）**: 受領したスクリプト一式には、
「Rt・Rpの組み合わせ自体をDFTで探索・最適化する層内ステップ」
（spec.md当初ドラフトの「Step 4b: 不均一傾斜精密化・N形用」に相当するもの）が**見当たらない**。
N形を探すには、ユーザーが手動でRt,Rpの値を推測してinit CSVに並べる運用になっていると思われる。

**軽微な指摘**: `step3_para.py` 冒頭のコメント「##tetracene層内計算」は、実際の中身
（層間=interlayer計算、`opt_param_keys=['cx','cy','cz']`とコメントに明記）と矛盾する。
他にも見つかった消し忘れコメント（`make_step2_twist.py`の"pbepbe+d3bj"）と同種の
ノイズと思われ、実害はないため未修正。

### Tab 4: Step 3 – Interlayer Stacking（para / twist サブタブ）

#### paraサブタブ — 計算手順（2026-07-06、コードで確認）

Step1（層内 a,b 最適化）と全く同じ「山登り法＋ジョブキュー管理」の構造を、
層間オフセット `(cx, cy, cz)` に対して適用したもの。`step3_para.py`（+ `make_step3_para.py`）。

1. **（前段・任意）vdW粗探索**: `step3_para_vdw.py`でvdW接触モデルにより層間距離
   z(x,y)の粗いマップを作り、cx,cy,czの初期値の当たりをつける（Step1aのvdW版と同じ発想）
2. **初期化**: `step3_para_init_params.csv`に固定パラメータ`(a,b,theta,Rt,Rp)`と
   初期の層間オフセット`(cx,cy,cz)`を1行以上用意（status='NotYet'）
3. **ジョブ投入**: 現在の`(cx,cy,cz)`から**10種類の層間ダイマー**を作り、1本の
   Gaussianジョブ（--Link1で10ペア連結）として投げる。status→'InProgress'
4. **結果回収**: ログから10個のBSSE補正済み相互作用エネルギーを読み取り合計:
   ```
   E = (E_i01+E_ip1+E_ip2+E_i02+E_ip3+E_ip4)/2 + E_it1+E_it2+E_it3+E_it4
   ```
   （2つの層パターン分は平均、T型4本は2パターンで共有="identical"のためそのまま加算）
5. **山登り（3×3×3=27点）**: cx,cy,czをそれぞれ±0.1Åずらした27通りでEを計算
   （既計算分は再利用）。最小のEを与える組を新しい中心にする
6. **収束判定**: 新しい最小点が動かなければ収束='Done'。動けば5に戻る（座標降下法）
7. **次の初期値へ**: 収束したら次の行（NotYet）へ。`--num-init`で並行探索数を指定
8. 全行'Done'で終了。`step3_para.csv`に全評価点の履歴＋各行の収束値が残る

**10種類のダイマーの内訳**:
- `i01, ip1, ip2`: パターン1（+A3の中心分子）から見た層間最近接1分子＋SP型2方向
- `i02, ip3, ip4`: パターン2（-A3、グライド対称のもう一方）の同じ組
- `it1〜it4`: T型接触4本（2パターンで共有）

#### step3_para_vdw.py — vdW粗探索段階（2026-07-06、コードで確認）

`step3_para.py`（DFT精密化）の**前段**。Step1aのvdW粗探索と同じ「粗い総当たり→精密な山登り」の
2段構成を、層間版で繰り返している。中心関数は`get_c_vec_vdw(monomer_name, R3, R4, a_, b_, theta)`
（コメント「R3:t-shaped R4:slipped-parallel」= 前述のRt,Rpと同じもの）:

- **固定（外から渡す）**: `a_,b_`（格子定数）、`theta`（=A3、ヘリンボーン角）、`R3`(=Rt)、`R4`(=Rp)
- **グリッドサーチ（0.1Å刻みの総当たり、山登りではない）**: `Ra, Rb`（層間オフセット、
  `±a/2`・`±b/2`の範囲）。これが`step3_para.py`の`cx,cy`の候補元になる
- 各`(Ra,Rb)`点での`z_max`（vdW接触が許す最小層間距離）は**解析的に計算**（サーチではない）。
  `cz`の候補元。`V = a_*b_*z_max`も同時に計算（スペック当初ドラフトの「V(x,y) vdW体積マップ」）
- `detect_peaks()`で`-z_max`の極大（=z_maxの極小、最も密に詰まる安定候補）を検出し、
  その`(Ra,Rb,z_max)`が`step3_para.py`の初期`(cx,cy,cz)`になると考えられる
- **グリッド点数**: `a_,b_`と固定の0.1Å刻みから自動的に決まる（a,b≈6〜8Åなら数千点程度）。
  ユーザーが選ぶものではない
- **`(a,b,theta,Rt,Rp)`の組み合わせを何パターン試すかはコード上未規定**。この関数は
  1回の呼び出しにつき1組の固定値のみを扱い、複数パターンを自動で振るループがない
  （Rt,Rp自動探索コードがないという前述の指摘と同じ構造）。何パターン試すかは
  大野さんの手作業判断と思われる — 明日確認したい

**バグ修正（2026-07-06）**: `get_monomer_xyzR`は5引数（Ta,Tb,Tc,A2,A3）を取るが、
`step3_para_vdw.py`の呼び出しは`get_monomer_xyzR(monomer_name,0.,0.,0.,0.,0.,theta)`と
6引数渡していて実行時エラーになる状態だった。余分な`0.`を1つ削除して修正し、
小さな格子定数でのスモークテストで正常動作を確認済み。

#### Tab 4 paraサブタブへのvdW事前スキャンGUI実装（2026-07-06、三好さんの提案で実装）

論文Fig.6(b-d)上段のV(x,y)マップ（vdW粗探索）→ 極小領域だけDFTで精密化、という
ワークフローをTab 2（Step1）と同じ「GUI内vdW実行→クリックでCLI用CSV生成」の
パターンでTab 4 paraサブタブに実装した。

- **`src/csp/vdw/interlayer.py`（新規）**: `step3_para_vdw.py`の`get_c_vec_vdw`を
  numpyでベクトル化して移植した`interlayer_vdw_scan()`。3つの異なる
  `(a,b,theta,Rt,Rp)`の組み合わせで元の（バグ修正済み）関数と数値誤差1e-15で
  完全一致することを検証済み。ペンタセンの実サイズ格子（a≈6.8,b≈6.6Å、
  4623グリッド点）で約5秒（元のPython素朴ループ版なら数分かかる想定）
- **`bilayer_preview()`**: 選択した点における2層構造（下層9分子クラスター＋
  上層1分子）を3Dプレビュー用に生成。論文Fig.6eの「overlayer(red) vs
  underlayer(gray)」に相当
- **GUI**: a,b,theta（Tab2のDFT結果があれば自動入力）、Rt,Rpを指定して
  「Run interlayer vdW scan」→ Ra-Rb平面のヒートマップ（z or V で色分け切替可）を
  scipy `minimum_filter`で検出した極小点（黒い四角マーカー、論文の記述通り）付きで表示。
  グラフをクリック（または黒四角をクリック）すると隣に2層構造の3Dプレビューが更新され、
  その点の`step3_para_init_params.csv`（1行）をダウンロードできる
  → そのままTab 4下部の「How to run (CLI)」のDFT精密化ステップの入力になる
- AppTestで実行・値切替・クリア動作を確認済み

**バグ修正（2026-07-06、三好さん実ブラウザ確認で発覚・修正）**:
1. ヒートマップ本体をクリックすると`customdata[1]`でIndexError、次に辞書型で
   `KeyError: 0`とクラッシュを繰り返した。原因は`px.imshow`のヒートマップと
   自作の散布図（極小点マーカー）で`customdata`の形がトレースごとに異なっていたこと。
   → **customdataを完全に廃止**。ヒートマップも散布図も`x=Ra, y=Rb`で描画しているので、
   クリックした点の`x, y`をそのまま使えば良いと気づき、シンプルに修正
2. 上記修正後も「黒い四角以外をクリックすると画面はリロードされるが3Dモデルが
   更新されない」不具合が発覚。原因調査として`px.imshow(pivot, ...)`（DataFrameを直接渡す形）
   にx,yを明示的に渡す修正を試したが、**それでも解決しなかった**（三好さん再確認で判明）
3. **根本原因（2026-07-06 最終解決）**: Plotly.jsにおいて、Heatmap/imshowトレースは
   単純な1クリックに対して「selection」イベント（Streamlitの`on_select="rerun"`が
   購読しているもの）を確実には発火しない（box/lasso選択向けの挙動で、マーカー系
   トレースとは扱いが違う）。そのためヒートマップ本体をクリックしても
   `event.selection.points`が空のままで、何も更新されなかった
   （黒四角の散布図オーバーレイだけは、Tab 2と同じ「マーカートレースのクリック」の
   仕組みなので正しく発火していた）。
   → **ヒートマップ全体をpx.imshowではなく、Tab 2と全く同じ「四角マーカーの散布図」
   （`go.Scatter(mode="markers", marker=dict(symbol="square", color=値, colorscale=...))`）
   で描画するよう変更**。全グリッド点（数千個）を1つの散布図トレースとして描画し、
   マーカーサイズはグラフ幅とグリッド点数から概算して隣接マーカーが概ね埋まるよう調整。
   これでクリック判定がTab 2と同じ実績のある経路になり解決
4. **縦の黒い線（2026-07-06、三好さんが原因を特定）**: Ra軸とRb軸の刻みは同じ0.1Åで
   データ上は正方形のはずだが、`scaleanchor`で縦横比を1:1に固定していなかったため、
   実際の描画ではRa軸とRb軸のピクセル密度が食い違い、固定サイズの正方形マーカーが
   一方の軸で隙間を作っていた（それが黒い線に見えていた）。
   → `fig3.update_yaxes(scaleanchor="x", scaleratio=1)`で両軸のpx/Å比を固定し、
   マーカーサイズの計算も単一の基準（`650 / max(グリッド点数)`）に統一して解決
5. **2層構造が重なって見える不具合（2026-07-06、三好さん報告・重大なバグ）**: 実際に
   原子間距離を検証したところ、`Rt=Rp=0`（既定値）では正しく接触（誤差1e-15）していたが、
   **Rt・Rpを0以外にすると実際に原子が重なっていた**（例: Ra=2.4,Rb=-1.1で
   距離が半径和より0.22Å短い＝物理的にありえない重なり）。
   原因は`interlayer_vdw_scan`が返す`z`列が「上下層をvdW接触させるのに**追加で**必要な
   すき間」であり、Rt,Rpに由来する層自体の内部シフト`z_shift = (2·Rt−Rp)·Ra/a + Rp·Rb/b`
   を含んでいなかったこと（Rt=Rp=0なら z_shift=0 なので気づかなかった）。
   Ono さんの元コード自身もこの`z_shift`抜きの値を返す設計で、csp側の数値検証はこれと
   ビット単位で一致させていたため、**元コードに忠実であるがゆえに見落としていたバグ**。
   → `interlayer_vdw_scan`の戻り値に**実際の層間c-vector z成分**`cz = z_shift + z`を
   新しい列として追加し、GUI側（3Dプレビュー・CSVエクスポート）は`z`ではなく`cz`を
   使うよう修正。Rt,Rp≠0の30点でランダム検証し、重なりが完全に解消（誤差1e-15）することを確認。
   **これは`step3_para_init_params.csv`としてダウンロードされる値にも影響していたため、
   Rt/Rpを使った探索（N形方向）を試していた場合はDFTの出発点が物理的に不正だった可能性がある**

### legacy/ono_scripts/XRD_pattern/（2026-07-06 大野さんアップロード・解読）

論文**§2.7「Comparison with Experimental Structures」（Fig.10、R形偏差）**に対応する
粉末X線回折（XRD）パターンのシミュレーション＋実験構造との比較プログラムと判明。

**構成**:
- `src/powder pattern model.ipynb` — 本体。クラス `structure_factor_model`
- `atom_scattering_factor_list/*.txt` — 元素ごとの原子散乱因子パラメータ
  （Cromer-Mann型、4ガウシアン+定数の9個の数値）
- `cryst/*.xyz` — **実験結晶構造**（CSD refcodeそのもの: ANTCEN09, PENCEN01/09/10,
  TETCEN01/04, NAPHTHA12, HEXCEN01）
- `opt_model/*.xyz` — **こちらの計算による最適化構造**（pentacene_ht/sc/tf,
  tetracene_sc/tf 等、多形ごと）

**アルゴリズム**:
1. `structure_factor_model.__init__(a1,z1,b2,z2,cx,cy,cz,gjf_name)`:
   **csp本体と全く同じ7変数**（a,b軸傾き＋層間オフセットcx,cy,cz）から格子ベクトルa,b,cを構築
   （サブクラス`structure_factor_cif`は代わりに標準結晶学パラメータa,b,c,α,β,γを直接受け取る
   → 実験構造の計算用）
2. `atom_factor(atom,t)`: 散乱因子ファイルを読み、温度因子（Debye-Waller的な等方性B因子、
   H: B=8π²×0.06、他: B=8π²×0.05、固定値）で補正
3. `structure_factor(h,k,l)`: 標準的な結晶構造因子 F(hkl) = Σ 原子散乱因子×exp(iK·R)
4. `make_csv`: (h,k,l)を-N〜N（N=4、コメントには本番はN=20000相当と示唆）で網羅し保存
5. `gen_theta`: Braggの式で(h,k,l)→2θに変換
6. `LP`: Lorentz偏光補正（3D/2D=薄膜モード切替可）
7. `laue`/`peak_sum`: 有限結晶子サイズによるピーク広がり（Laue関数）を全反射で足し上げ、
   連続的な粉末パターンI(2θ)を生成・プロット

**現状**: CLAUDE.mdの「使わないもの」リスト（Amber力場、遠距離相互作用外挿）には
**含まれておらず**、tot_energy/とは違いスコープ外の決定はまだされていない。
csp本体・GUIへの統合は未着手（三好さんの判断待ち）

### Tab 2 Fig.2(b)再現機能（2026-07-06、三好さんの提案で実装）

`step1.csv`だけでは「各αごとの2〜3個の種（a-stack/b-stack/local min、Tab2上部の
vdW事前スキャンと同じ分類）のうちどれに属する点か」も「各種がどこに収束したか」も
分からない。両方とも新規計算なしで復元可能:
- 種の分類は`csp.vdw.contact.step1a_scan`の`df_init`（`kind`列）をαごとに再計算するだけ
  （大野さんの`init_process()`と同じvdW粗探索、値が一致することを検証済み）
- 各種の収束先は、`step1.py`の`get_opt_params_dict`と全く同じ3×3近傍降下法を、
  新規Gaussianジョブの代わりに`step1.csv`の(theta,a,b)→Eルックアップに対して
  再生（リプレイ）するだけで求まる（本物の山登り法が訪れた点は全て`step1.csv`に
  既にあるので、リプレイは必ず同じ収束点に辿り着く）

**新規実装**: `src/csp/plot/step1_results.py`の`classify_and_fold_step1_results()`。
さらに、ポリアセンはα↔(90-α)でa↔bが厳密に入れ替わる対称性を持つことをvdW接触モデルで
検証済み（`vdw_R(θ,'a')==vdw_R(90-θ,'b')`）。これは同一の物理構造をa/bの呼び方を
入れ替えただけなので、α=5〜45の範囲のスキャンを折り返すだけで論文Fig.2(b)相当の
広い範囲（0〜90°）が新規計算なしで再現できる。Tab2の「DFT results」の下に追加し、
プロットの点をクリックすると9分子クラスターの3Dプレビューが表示される
（既存のFig.2b風vdW事前スキャンと同じUIパターン）

### Tab 3 paraの再設計（2026-07-06、三好さんの指摘・提案で実装）

三好さんの指摘: 従来のEt(z)/Ep(z)グラフに表示していた「4·Et+2·Ep」の破線は
**G形（グライド対称、全接触に同じzを使う場合）にしか成り立たない式**で、
一般的な表示として不適切だった。削除し、代わりに以下2機能を追加:

1. **Et/Ep曲線 + クリックで3Dダイマープレビュー**: グラフ右に3Dモデルを追加し、
   Et(z)曲線上の点をクリックすると`(0,0,0)-(a/2,b/2,z)`のT型ダイマー、
   Ep(z)曲線上の点をクリックすると`(0,0,0)-(0,b,z)`または`(0,0,0)-(a,0,z)`
   （`step2_para.py`本来の規則通り、a>bならb方向、そうでなければa方向）の
   SP型ダイマーを表示。`csp.structure.intralayer.dimer()`に`z`引数を追加して対応
   （デフォルト0、既存呼び出しへの影響なし）
2. **Fig.5(b)風θ_incl・φ_inclマップ**: 「N形2次元マップの再構成方法」の式を
   `src/csp/plot/step2_results.py::build_theta_phi_map()`として実装。
   ヒートマップはTab4 vdWマップと同じ理由（Plotlyのheatmap/imshowはクリック
   selectionを確実に発火しない）で、正方形マーカーの散布図方式を最初から採用。
   クリックすると`csp.structure.intralayer.cluster6_inclined()`
   （新規関数、6隣接分子+中心の傾斜クラスター、N形2次元マップの導出式そのもの）
   で6分子クラスターを表示。a,b,theta入力欄を追加し、Tab2の選択結果
   （`s1fig2b_current`/`s1vdw_current`）があれば自動で初期値に反映

**設計上の注意**: `build_theta_phi_map`はOnoさんの`plot2d()`と同じ規約
（SP接触は常に"b方向"の役割で式に組み込む、a>bによる切替なし）を採用。
一方、Et/Ep曲線のダイマープレビューは実際の`step2_para.csv`を生成した
`make_step2_para.py`が使うa>b切替規則に忠実に従う。目的が違う2つの表示で
異なる規約を使っているが、それぞれの参照元コードには忠実

### Tab 2 vdWグラフの不具合修正（2026-07-06、三好さん報告）

三好さんがペンタセンで実際に触って発見した3つの問題:

1. **Fig.S1(c)サブタブがスキャン未実行時に表示されない**: `step1_init_params.csv`
   だけでは候補点しか無く、β全域スイープ曲線（`df_curves`）が無いため
   「スキャンを実行した直後しか使えません」という案内文のままだった。
   → `example/pentacene/step1_vdw_curves.csv`（`step1a_scan()`の副産物）を追加で
   バンドルし、サンプル表示時（`step1_init_params.csv`をアップロード扱いでなく
   サンプルとして表示している時）はこちらも自動で読み込むよう修正
2. **Fig.2(b)のDFT版（新機能）でθ=23〜39あたりが「ジグザグ」に見える**:
   原因は`local_min`の山登り結果が実は`b_contact`（またはa_contact）と
   **全く同じ収束点**になっていたこと（例: θ=30で両方とも(7.7,5.7)に収束）。
   θ=39→40で急に別の点にジャンプする（本当に別構造になる）ため不連続に見えた。
3. **local minとa-stack/b-stackの重複表示**: 2と同じ原因。
   「最終的にstack構造に落ち着いて極小点が2つしかない時はa-stack, b-stackだけ
   表示して」という三好さんの要望通り、`classify_and_fold_step1_results()`で
   a-stack/b-stackを先に処理し、`local_min`の収束点がどちらかと厳密に一致する
   場合はその`local_min`行を捨てるよう修正（ペンタセンの実データでθ=23〜39が
   該当し、θ=40〜45だけが真に別構造として残ることを確認）。
   なお、vdW粗探索段階（DFT前）でも同種の重複が起き得るため
   `step1a_scan()`側にも同じ考え方の重複除去（局所S極小点がa/b端点の
   (a,b)と一致する場合は除外）を追加済み

**表示順の見直し（同日、三好さんの提案）**: Fig.2(b)のDFT再構成が使えるようになった今、
従来の単純な「全パラメータでの最小値」チェックボックス付きプロットは冗長と判断し撤去。
molecule未選択時のフォールバック表示としてのみ残す（`agg_min`オプションは削除、
常に生データをそのまま表示）。分子選択済みの通常時はFig.2(b)ブランチ分類プロットが
`step1.csv`アップロード直後から即座に主表示になるよう並び替えた

### 論文SI Table S2の値とFig.5(b)検証（2026-07-06、三好さんの提案で調査）

論文SIの`Table S2`（Type IV of tetracene; Type I, II, IV of pentacene; Type II of
hexacene の最適化パラメータ表）からペンタセンの値を抽出:

| 構造 | a (Å) | b (Å) | α (deg) | θ_incl (deg) | φ_incl (deg) | cx (Å) | cy (Å) | cz (Å) |
|---|---|---|---|---|---|---|---|---|
| Type I (R形, calc) | 7.2 | 6.0 | 25 | 0.0 | − | 0.0 | 1.9 | 15.4 |
| Type II (単結晶相, calc) | 7.2 | 5.9 | 25 | 27 | 48 | 0.2 | 1.3 | 16.0 |
| Type IV (calc) | 7.4 | 5.8 | 25 | 22 | 43 | 1.5 | 1.3 | 15.6 |

（cx,cy,czはx,y,z列に対応。θ'_inclは不均一傾斜のパラメータで今回は使わず）

**検証**: 受領済みの実`step2_para.csv`を使い、`build_theta_phi_map(df, a=7.2, b=6.0)`で
2次元マップを作り、`scipy.ndimage.minimum_filter`で局所極小点を探索したところ
**θ_incl=25.4°, φ_incl=47.6°**という極小点が見つかり、論文のType II値
（27°, 48°）と非常に近い（グリッド分解能0.1Åの制約を考えれば妥当な一致）。
→ **Fig.5(b)再構成のロジック自体は正しいことを確認**。三好さんが「論文と全然違う」と
感じたのは、この極小点がR形（平坦、θ_incl=0）のグローバル最小値（E=-73.4）より
浅い（E=-70.8）ためプロットの色スケールに埋もれて見えなかったのが原因。
Tab4のvdWマップと同じ`minimum_filter`による局所極小点マーキング（黒四角）を追加して解決。
なお、論文のFig.5(b)自体もR/G/N形はいずれも局所極小点であり、Rが必ずしも
グローバル最小というわけではない（本文Fig.5参照）ので、この挙動は物理的に妥当

**Tab3・Tab4のデフォルト値も修正**: 従来はTab2で何も選択していない場合、
a,b,thetaが根拠のない仮の値（7.2,6.0,22.0など）だったため、サンプル表示時の
3Dプレビューが意味をなさなかった。上記Table S2の値（Tab3=Type II、
Tab4 vdW事前スキャン=Type I）をフォールバックのデフォルト値として設定

**軸の取り違えバグ修正（2026-07-06、三好さん指摘）**: 三好さんが「論文の図は横軸-45〜45、
縦軸-30〜30なのに全然違う」と指摘。大野さんの`plot2d()`をよく見ると、実際にプロットしている
軸は生の`(phi_incl, theta_incl)`ではなく、**極座標→直交座標変換**だった:
```python
x = theta_incl * cos(phi_incl)   # θ_inclを半径、φ_inclを角度とする
y = theta_incl * sin(phi_incl)
```
平坦なR形（θ_incl=0）が原点に来て、そこからの半径方向の距離が傾きの大きさ、
という結晶学的に自然な「傾きマップ」になる。これを見落として生の(φ,θ)を軸に
していたのが「論文と全然違う」の原因。`build_theta_phi_map()`に`x`,`y`列を追加し、
プロットの軸をこちらに変更。ペンタセンの実データで検証したところ
x範囲は±48°（論文±45°）、y範囲は±34°（論文±30°）とほぼ一致することを確認

### vdW版Fig.2(b)も90°まで折り返し表示（2026-07-06、三好さんの提案）

Tab2のvdW事前スキャン（DFT前の粗い vdW モデルによる S vs α プロット）にも、
DFT版Fig.2(b)再構成と同じα↔(90-α)・a↔b入れ替えの折り返しをチェックボックスで追加
（デフォルトON）。これにより既存の5〜45°スキャンだけで論文同様0〜90°の全域を表示できる

### Tab3の3Dモデルのα表示バグ修正（2026-07-06、三好さん報告）

三好さんが「Tab3の3Dモデルがα=25°に見えない」と報告。まず幾何計算式自体を検証
（`dimer()`, `cluster9()`, `cluster6_inclined()`いずれもPCAで分子平面の法線同士の
なす角を直接計算）した結果、**θ=25°を渡すと厳密に50.0°になることを確認**
（式は正しい）。

**本当の原因**: Tab3の`theta`入力欄が実際には`5`になっていた。Tab2のFig.2(b)は
未クリック時に「グローバル最小値（R形の最適点）」をデフォルト表示しているが、
その値を`session_state`（`s1fig2b_current`）に**保存していなかった**ため、
Tab3は「クリック履歴」を探しても何も見つからず、データの最初の行
（θ=5、必ずしも最適点ではない）由来の値にフォールバックしていた。

**修正**:
1. Tab2のFig.2(b)・vdW版Fig.2(b)ともに、未クリック時のデフォルト表示（グローバル最小値
   / 最小S点）を計算した時点で即座に`session_state`（`s1fig2b_current`/`s1vdw_current`）に
   保存するよう変更。これによりクリックしなくても「今表示されている最適点」が
   Tab3・Tab4に正しく引き継がれる
2. グローバル最小値が、α≤45の元データとα>45の折り返し（対称性による鏡像、
   物理的には全く同じ構造）の両方に**同値で存在する場合、折り返し前
   （α≤45、論文で一般的に引用される範囲）を優先して表示**するよう修正
   （なお折り返し後の表現でも3D構造自体は完全に同一なので、以前の表示も
   間違ってはいなかった。表示上の分かりやすさの問題）
3. vdW版のデフォルト表示ロジック（`_default_current`）も、従来は「先頭行」を
   使っていたのを「S最小の行」を使うよう修正し、同様に`session_state`へ保存するよう統一

### Tab 3 twistをFig.7(c)風の2次元マップに（2026-07-06、三好さんの提案）

三好さんから「Tab3のtwist結果はどうにかできる？Fig.7のcみたいな図を作りたい」と依頼。
論文Fig.7(c)は「naphthalene・anthraceneについて、θ_twist・ΔzT(=Rt)を軸にした
Eintra(6)の2次元マップ」。

**現状**: 大野さんからは`step2_twist.csv`（naphthalene/anthracene）はまだ未受領。
三好さんが進めているアントラセンのHPC計算（Step1）も、今回SSH接続がタイムアウトして
状況を確認できなかった（一時的なネットワークの問題の可能性）。

**先に実装したこと**: データが届いた時にすぐ使えるよう、Tab3 twistの表示を
Fig.2(b)/Fig.5(b)と同じパターンの**2次元マップ+局所極小点マーキング+クリックで
3Dダイマープレビュー**に書き換えた。
- `make_step2_twist.py`の幾何規約（T型ダイマー: 中心`(0,0,0,A2,+theta)`、
  隣接分子`(a/2,b/2,Rt,A2,-theta)`）は既存の`dimer()`関数（`z`引数を今回のTab3 para
  改修で追加済み）がそのまま使えることを確認、新しい幾何関数は不要だった
- (A2,Rt)グリッドの各点で(a,b)も山登り最適化されるため、`groupby([A2,Rt])[E].idxmin()`
  で各グリッド点の収束後の最小Eを取り出してからマップ化
- 合成データ（既知の最小値を仕込んだ疑似`step2_twist.csv`）で、2次元マップの
  極小点検出・3Dダイマー構築のロジックが正しく動くことを検証済み
- **潜在バグ修正**: Tab3 twistの3Dプレビューが、隣の`sub_para`ブロックでのみ
  条件付きで定義される`s2_a`/`s2_b`/`s2_theta`変数を参照していたため、
  `step2_para.csv`が無い状態で先にtwistサブタブを開くと`NameError`になる
  潜在的なバグがあった。twist専用の独立したフォールバック変数
  （`_twist_a`/`_twist_b`/`_twist_theta`）に切り替えて解消

**残作業**: 実データ（大野さんから受領 or 自前のアントラセン計算完了）が届き次第、
このTab3 twistの表示で実際にFig.7(c)相当の図が再現できるか検証する

### Fig.5(b)・Tab3 twistマップのクリックが反応しない不具合修正（2026-07-07）

三好さんが「Fig.5(b)の2次元マップをクリックしても3Dモデルが変わらない」と報告。
原因は`customdata`ベースのクリック識別で、これはTab4のvdWマップで以前遭遇したのと
同種の信頼性問題（過去にheatmap/imshowトレースでは`customdata`が確実に機能しない
ことが分かっていたが、今回は素朴な散布図トレースでも同様の問題が起きた）。

**修正方針**: `customdata`に頼らず、クリックされた点の生の`x, y`をそのまま
`df_map`自身の`x, y`列と`np.isclose`で突き合わせて対応する行を探す方式に統一
（この方式はTab2のFig.2(b)・Tab4のvdWマップで既に実績があり、最も安定している）。
Fig.5(b)と、同じ`customdata`方式で実装したばかりのTab3 twistマップの両方を修正。

**教訓**: Plotlyの`customdata`によるクリック識別は、このアプリの環境では
繰り返し信頼性の問題を起こしている（ヒートマップだけでなく単純な散布図でも）。
今後クリックで点を識別する新機能を作る際は、**最初から`customdata`を使わず、
プロットしている実際の`x, y`値をデータフレーム自身の列と照合する方式**を
デフォルトの実装パターンとすること

### 3Dモデルが枠の中心にない不具合修正（2026-07-07、三好さん報告）

三好さんが「3Dモデルが枠の中心にない、全ての枠で直して」と報告。原因は
`src/csp/plot/viewer3d.py`の`render_3d_html()`が生成するpy3Dmolキャンバスが
**固定サイズ**（480×420px）なのに対し、呼び出し側の`st.iframe()`は明示的な幅を
渡さず**カラム幅いっぱいに伸びる**設計だったこと。カラム幅が480pxと異なる場合
（特に`st.columns([2,1])`の狭い方の列など）、固定サイズのキャンバスが
左寄せで表示され中央に来ていなかった。この不整合は最初のスキャフォールディング
（プロジェクト最初期）からずっと存在していた模様。

**修正**: `render_3d_html()`の戻り値を`display:flex; justify-content:center;
align-items:center;`のdivで包み、キャンバスサイズと外枠サイズが一致しなくても
常に中央に来るようにした。この関数は`render_molecule_3d()`経由で全8箇所の
3Dプレビューから共通で呼ばれているため、**1箇所の修正で全ての枠に反映される**

### Tab 5: Transfer Integrals — DFTで何を計算しているか（2026-07-06、コードで確認）

`tcal_csv.py`の`pre_process()`は、1配置（`a,b,theta,A2,z`の組）× 接触タイプ（T型/SP型）
ごとに**6本のGaussianジョブ**を作る。全6本とも同じルート行`#b3lyp/6-31g*`
（2026-07-05修正後の値。受領時は`#pbepbe/6-311G**`だった）を使用:

| ファイル | 対象 | 目的 |
|---|---|---|
| `test_t_m1.gjf` | 単量体1だけ | HOMO軌道を得る |
| `test_t_m2.gjf` | 単量体2だけ | HOMO軌道を得る |
| `test_t.gjf` | ダイマー | フォック行列・重なり行列を得る |
| `test_p_m1.gjf`, `test_p_m2.gjf`, `test_p.gjf` | 同上（SP型接触版） | 同上 |

各ファイルは`--Link1--`で2段階：
1. 通常のSCF計算（電荷・多重度`0 1`）。波動関数をチェックポイントに保存
2. `geom=allcheck guess=read pop=full iop(3/33=4,5/33=3)` —
   チェックポイントを読み直し、**重なり行列**（IOp 3/33=4）と**フォック行列**（IOp 5/33=3）を出力

`tcal_1.py`がこれらを読み、ダイマーのフォック行列を単量体HOMOに射影する標準的な
手法（松井研tcal本体）でtransfer integral Jを算出する。1配置あたり計6本×2接触タイプ
= 実質1配置で処理は上記の通り（T型3本＋SP型3本）。

---

## GUI タブ構成 / Streamlit UI Layout

> **注意（2026-07-05）**: 以下は論文ドラフトベースの初期ワイヤーフレーム（実コード受領前）。
> **現在の実装は5タブ構成**（上の「大野コード対応表」参照）で、Tab番号・構成とも
> ここに書かれた内容とは一致しない。実装の正としては使わないこと。歴史的記録として残す。

### Tab 1: Molecule Setup

```
[Preset] naphthalene / anthracene / tetracene / pentacene / hexacene
[Upload] Custom XYZ file
─────────────────────────────────────────────────
Left: vdW radius table (editable)    Right: 3D viewer (py3Dmol)
      atom | vdW radius (Å)
      C    | 1.70
      H    | 1.20
      S    | 1.80
      N    | 1.55
      O    | 1.52
```

### Tab 2: Step 1 – Intralayer vdW Scan

```
[Settings]
α range: [5] to [85] deg, step [5] deg
vdW tolerance: [0.0] Å  (overlap threshold)

[Run vdW Scan]

─────────────────────────────────────────────────
Left: α vs S=a×b plot (Plotly)      Right: layer structure preview (2D)
      click a point → show structure      T-shaped / SP contacts visualized
      [Export CSV]
```

### Tab 3: Step 1 – DFT-D Optimization

```
[HPC Settings]
Host: [miyoshi@133.11.68.31]   Work dir: [/home/miyoshi/Working/ono_paper_dir]

[Gaussian Keywords]
# B3LYP empiricaldispersion=gd3 6-311g** counterpoise=2 nosymm

[Generate & Submit Jobs]  [Check Status]

─────────────────────────────────────────────────
Left: α vs E_intra(8) plot (kcal/mol)   Right: 3D structure viewer
      click → show structure                   R-form highlighted (α≈25°)
      [Confirm R-form α]
```

### Tab 4: Step 2 – Inclination Map

```
[Settings]
θ_incl range: [0] to [40] deg, step [2] deg
φ_incl range: [0] to [360] deg, step [5] deg
Starting structure: [R-form from Step 1]  α = 25°

[Run DFT Jobs]  [Check Status]

─────────────────────────────────────────────────
Left: 2D heatmap E_intra(6)              Right: 3D structure viewer
      x-axis: θ_incl cos(φ_incl)               click point to display
      y-axis: θ_incl sin(φ_incl)
      markers: ● R-form  ▲ G-form  ■ N-form

[Confirm] θ_inclN=[  ] φ_inclN=[  ] θ_inclG=[  ]
```

### Tab 5: Step 3 – Interlayer Stacking

```
[Settings]
Form: [R-form] / [G-form] / [N-form]
x range: [-a/2] to [a/2], step [0.2] Å
y range: [-b/2] to [b/2], step [0.2] Å

[Run vdW Volume Scan]  [Run DFT Einter(7)]

─────────────────────────────────────────────────
Left: V(x,y) map (vdW volume)    Right: E_inter(7) map (DFT-D)
      ● N1 (Type II)                    ● N1  ■ N2
      ■ N2 (Type IV)

[Select stacking for Step 4]  N1 / N2
```

### Tab 6: Step 4 – Refinement

**Sub-tab 6a: Twist (G-form → Type III)**
```
[Settings]
θ_twist range: [0] to [20] deg, step [1] deg
Optimize: a, b, ΔzT, x, y, z

[Run DFT Jobs]

─────────────────────────────────────────────────
Plot: θ_twist vs E_int(near)    Right: 3D structure
      optimal θ_twist marked
```

**Sub-tab 6b: Non-uniform Inclination (N-form → Type IV)**
```
[Settings]
θ'_incl range: [0] to [5] deg, step [0.5] deg
φ'_incl: [auto from Step 4 minimum]

[Run DFT Jobs]

─────────────────────────────────────────────────
Plot: θ'_incl vs E_int(near)    Table: Type II vs Type IV comparison
      for tetracene/pentacene/hexacene
```

### Tab 7: Transfer Integrals

```
[Select structures]
  Type I (R-form)  ☑
  Type II (N1)     ☑
  Type III (G+twist) ☑
  Type IV (N2)     ☑

[Run Gaussian (B3LYP/6-31g*)]  [Check Status]

─────────────────────────────────────────────────
Bar chart: J (meV) for T-type and SP-type contacts per polymorph
Plot: J vs α (R-form, varying α)
Plot: J vs θ_incl (N-form, varying θ_incl)
```

---

## ファイル構成 / Directory Structure

```
ono_paper/
├── src/
│   └── csp/
│       ├── __init__.py
│       ├── structure/
│       │   ├── molecule.py      # XYZ読込, PCA軸合わせ, vdW半径割り当て
│       │   ├── intralayer.py    # α,a,b → グライド対称層構造生成
│       │   └── interlayer.py    # x,y,z → 上下層スタッキング構造生成
│       ├── vdw/
│       │   ├── contact.py       # vdW接触判定, a,b 最小化
│       │   └── volume.py        # V(x,y) マップ計算
│       ├── dft/
│       │   ├── make_gjf.py      # Gaussian16 入力ファイル (.gjf) 生成
│       │   ├── parse_log.py     # .log からエネルギー・MO 抽出
│       │   └── job_cluster.py   # SGE ジョブ投入 (auto_opt/cluster.py を流用)
│       ├── transfer/
│       │   └── calc_J.py        # 松井先生コード組み込み: HOMO→J計算
│       └── plot/
│           ├── map2d.py         # Plotly インタラクティブ 2D マップ
│           └── viewer3d.py      # py3Dmol 3D 分子表示
├── app.py                       # Streamlit メイン (Tab 1〜7)
├── data/
│   └── molecules/
│       ├── naphthalene.xyz      # B3LYP-D3/6-311G** 最適化済み (大野さんから取得)
│       ├── anthracene.xyz
│       ├── tetracene.xyz
│       ├── pentacene.xyz
│       └── hexacene.xyz
├── example/
│   └── pentacene/               # 論文再現用 計算済み結果 (初期公開時に同梱)
│       ├── step1_vdw.csv        # α vs S
│       ├── step1_dft.csv        # α vs E_intra(8)
│       ├── step2_map.csv        # θ_incl, φ_incl, E_intra(6)
│       ├── step3_N1.csv         # x, y, V, E_inter(7) for N1
│       ├── step3_N2.csv
│       └── transfer_integrals.csv
├── pyproject.toml
├── requirements.txt
└── README.md                    # 英語, インストール〜実行手順
```

**auto_opt から流用できるコード**:
- `auto_opt/src/auto_opt/utils.py` の `place_monomer()` → `structure/intralayer.py` に組み込み
- `auto_opt/src/auto_opt/cluster.py` → `dft/job_cluster.py` にコピー・整理

---

## 依存ライブラリ / Dependencies

```toml
[project.dependencies]
numpy = ">=1.24"
scipy = ">=1.10"       # a,b 最適化 (minimize)
pandas = ">=2.0"
plotly = ">=5.0"       # インタラクティブマップ
streamlit = ">=1.30"
stmol = ">=0.0.9"      # py3Dmol Streamlit 連携
```

---

## 実装優先順位 / Implementation Phases

| Phase | 内容 | 目安 |
|-------|------|------|
| 1 | Step 1a (vdW scan) + Step 2 (傾斜マップ) の計算エンジン + Tab 2,4 UI | 最優先 |
| 2 | Step 1b + Step 3 の DFT ジョブ連携 + Tab 3,5 UI | 次 |
| 3 | Step 4a, 4b 精密化 + Tab 6 UI | |
| 4 | Transfer integral (松井先生コード統合) + Tab 7 UI | |
| 5 | README 英語化 + example/ 計算結果同梱 + GitHub 公開 + Zenodo DOI | 論文投稿前 |

---

## 残タスク / TODO（2026-07-05 時点）

**後回し（保留中）**:
- [x] ~~論文図のサンプル結果 CSV を `example/pentacene/` に同梱~~
      → **2026-07-06 一部完了**。大野さんから受領したペンタセン実データ2本
      （`step1.csv` 2692行、`step2_para.csv` 81行）を`example/pentacene/`に配置し、
      GUI（Tab 2・Tab 3 para）で正しく読み込めることを確認済み。
      さらに`step1_init_params.csv`（vdW初期候補）は大野さんから受領していないが、
      `step1.py`の`init_process()`をペンタセンでローカル実行するだけで再現可能
      （DFT不要、数秒で完了）と分かったので生成して配置済み。受領済みstep1.csvの
      22候補と完全一致することも検証済み。
      残り4ファイル（`step2_twist.csv`, `step3_para.csv`, `step3_twist.csv`,
      `transfer_integrals.csv`）は大野さんからの受領待ち
      （`example/pentacene/README.md`に受領状況を一覧化済み）

**受領・回答待ち（大野さん）**:
- [x] ~~モノマー CSV の座標規約照合~~ → 2026-07-06 実CSV受領・数値比較して確認済み
      （X・Y完全一致、Z符号反転だが分子対称性により無害。「大野コード対応表」参照）
- [x] ~~tcal で pbepbe を使っていた経緯の確認~~ → **2026-07-06 大野さんに確認済み: B3LYPで正しい**
      （論文Fig.11はb3lyp/6-31g*で計算されていた。2026-07-05の修正はそのままでよい）
- [x] ~~`step2_para.py` の `para_list.append(z,E1,E2)`~~ → 2026-07-05 修正済み。
      ただし2026-07-06、大野さんが自身のローカル版を再アップロードした際にこの修正が
      巻き戻ってしまい、マージで気づいて再修正（`legacy/`はGitHub上で直接編集する方針を
      共有した方がよいかもしれない）
- [x] ~~Tab 3 para の z シフトが論文のどの図に対応するか~~ → 2026-07-06 論文§2.2本文で確認済み
      （「タブ別詳細ガイド」Tab 3 参照。z = θ_inclの物理的実体、SI式で変換可能）
- [x] ~~`step3_para.py` の `Rt`・`Rp` が論文§2.2のφ_incl中間値（N形）に対応するか~~
      → 2026-07-06 コードで確認済み（Rp = Rt1-Rt2、「タブ別詳細ガイド」Tab 3 参照）。
      ただし Rt,Rp は自動最適化されず手動スキャン設計である点は大野さんに確認したい
      （山登り対象に含めない設計意図か、単に未実装か）
- [x] ~~`legacy/ono_scripts/tot_energy/`をcsp本体に統合するか~~
      → **2026-07-06 三好さん決定: 統合しない**。csp本体には組み込まず、
      README等で存在を「紹介」するに留める（CLAUDE.mdの当初方針通りスコープ外）
- [x] ~~受領コード一式にN形を発見・最適化するロジックが見当たらない件~~
      → **2026-07-06 完全解決**。三好さんの導出（線形代数で検証）と、同日大野さんが
      アップロードした`step2_para.py`の新`plot2d()`関数が完全一致。詳細は
      「N形2次元マップの再構成方法」参照。専用の2次元スキャンスクリプトは存在せず、
      `step2_para.py`の1次元Et(z)・Ep(z)カーブを組み合わせて再構成する方式と確定

### N形2次元マップの再構成方法（2026-07-06、三好さんの導出・大野さんの実コードで確認）

**2026-07-06 16:33、大野さんが`step2_para.py`に`plot2d()`関数を追加してアップロード。
三好さんの以下の導出と完全に一致する実装だった**（`--plot`フラグで実行）:
```python
Et1 = df[df['z']==z1]['Et'].values...   # z1 = zt1
Et2 = df[df['z']==round(z1-z2,1)]['Et'].values...  # z1-z2 = zt1-zp（偶関数性でzp-zt1と同値）
Ep  = df[df['z']==z2]['Ep'].values...   # z2 = zp
E = 2*(Et1+Et2+Ep)
```
さらに`(zt,zp)`→`(θ_incl,φ_incl)`への変換式も判明:
```python
za=2*zt-zp; zb=zp
Z=1/sqrt(1+(za/a)**2+(zb/b)**2); theta_incl=degrees(acos(Z))
phi_incl=phase(za/a + zb/b*1j)
```
これで「専用の2次元スキャンスクリプトは存在せず、step2_para.pyの1次元スキャン結果を
組み合わせて2次元マップを再構成する」という仮説が**確定**した。TODOの「大野さんに確認」
項目はこれで解消（さらに確認するまでもなく実コードで証明された）。

**同時に`step3_para_vdw.py`も大幅更新され、戻り値が`z_2dlist`から
`DataFrame(columns=['cx','cy','cz','V'])`に変更され、**`cz`列の計算式が
`cz=z_max+(2*Rt-Rp)*cx/a+Rp*cy/b`**——これは今日csp側で`interlayer_vdw_scan`に
追加した`cz = z_shift + z`と全く同じ式。三好さんの重なりバグ修正が独立に正しかったことの
裏付けになった。ただし新しい`get_monomer_xyzR1/2`関数（層をZ>0/Z<0で分割する設計に変更）に
旧知のバグ（引数過多、モノマーパスのハードコード）が再度混入していたため2026-07-06に再修正。
`step2_para.py`側も`para_list.append(z,E1,E2)`と`len(...)==41 & len(...)==41`
（`&`と`and`の優先順位バグ）が再度混入しており、あわせて修正

`step2_para.py`は中心分子を`z=0`に固定し、T型接触1本・SP型接触1本それぞれの
相対z変位に対するダイマーエネルギー Et(z), Ep(z) を独立に計算している
（Step1の9分子クラスター座標系: T型は`(±a/2,±b/2)`, SP型は`(0,±b)`または`(±a,0)`）。

**ポイント**: ペアワイズなダイマーエネルギーは「2分子間の相対位置」だけで決まるので、
Eintra(6)を計算する際、全T型・SP型接触に**同じz**を一律に使う必要はない。
各接触ごとに**異なる**z値でEt/Epを評価し、組み合わせて良い。

**一様傾斜（θ_incl, φ_incl）は平面 `z(x,y) = kx·x + ky·y` の式そのもの**なので、
6分子クラスターの各位置での必要なzは以下のように決まる（kx, kyがθ_incl,φ_inclに対応）:

```
中心      (0,      0,    0)
T型 a/2,b/2   → z = kx·a/2 + ky·b/2  = zt1
T型 -a/2,b/2  → z = -kx·a/2 + ky·b/2 = zp − zt1   (zp = ky·b と定義)
T型 a/2,-b/2  → z = kx·a/2 − ky·b/2  = zt1 − zp
T型 -a/2,-b/2 → z = -kx·a/2 − ky·b/2 = −zt1
SP型 0,b      → z = ky·b             = zp
SP型 0,-b     → z = -ky·b            = −zp
```

三好さんが提案したこの6点の座標を上式に代入して**完全に一致することを確認済み**
（線形平面モデルとの整合性を検証）。分子の対称性からEt(z)=Et(−z), Ep(z)=Ep(−z)
（偶関数）と仮定できるため、4つのT型は`Et1=Et(zt1)`と`Et2=Et(zp−zt1)`の2値、
2つのSP型は`Ep=Ep(zp)`1値に集約され、

```
Eintra(6) = 2×(Et1 + Et2 + Ep)
```

で計算できる（`zp=0`のときEt1=Et2=Et(zt1)となり、通常の`step2_para.py`が出す
`E=4Et+2Ep`のグライド対称ケースに一致する）。`(zt1, zp)`を独立に振ることで
任意の`(kx,ky)`＝任意の`(θ_incl,φ_incl)`（N形の中間角度を含む）を、
`step2_para.py`の1次元スキャン結果だけから再構成できる。

**TODO（2026-07-06、三好さん指示）**: このN形2次元マップ（論文Fig.5b相当）の
GUI実装をやった方がよい。Tab 3 paraで得られるEt(z)・Ep(z)の1次元スキャン結果
（`step2_para.csv`）から、上記の組み合わせ計算で(θ_incl,φ_incl)平面上の
Eintra(6)ヒートマップを再構成して表示する機能。Tab 4 paraのvdW事前スキャン
機能と似た位置づけ（既存の計算結果を活用した可視化）。実装時期は未定。

**格子定数 a, b は全過程で固定**（`step2_para.py`の`fixed_param_keys`にa,b,thetaが
含まれ、Step1の最適値をそのまま使う設計であることからも確認済み）。
- [x] ~~`step3_para_vdw.py`の`get_monomer_xyzR`呼び出しの引数過多バグ~~ → 2026-07-06 修正・動作確認済み
- [ ] **新規（2026-07-06）**: `(a,b,theta,Rt,Rp)`の固定パラメータの組み合わせを実際何パターン
      試したのか（コード上は1回の呼び出しにつき1組で、複数パターンを自動で振るループがない）
- [ ] **新規（2026-07-06）**: Fig.7(c)は論文で「a, b, ΔzTの3変数を最適化」と書かれているが、
      `step2_twist.py`は`Rt`(=ΔzT)を固定値として渡し`a,b`のみ自動最適化する設計。
      Rt自体を手動グリッドで振って最小値を選ぶ運用だったのか確認したい
      （Tab4のRt/Rp未探索問題と同根）
- [ ] **新規（2026-07-06）**: 論文§2.5（Fig.8、θ'_incl・φ'_inclの不均一傾斜、N2→Type IV）に
      対応するコード・タブが未特定。Tab4 paraのRp≠0探索と同じ仕組みなのか確認したい
- [ ] **新規（2026-07-06）**: Tab3 twist（Fig.7(c)相当）の実データがまだ無い。
      大野さんに`step2_twist.csv`（naphthalene/anthracene）をもらえないか確認、
      または三好さんのアントラセンHPC計算がStep2 twistまで進み次第、表示を検証
- [x] ~~`legacy/ono_scripts/XRD_pattern/`をcsp本体に統合するか~~
      → **2026-07-06 三好さん暫定決定: おそらく統合しない**（tot_energy/と同様、
      csp本体には組み込まず紹介に留める方向。最終確定ではないので変更の可能性あり）

**公開前の作業**:
- [ ] Zenodo 連携 → GitHub Release → DOI 取得（三好さんが対応。著者表記は Ohno で確定済み）
- [ ] CITATION.cff / .zenodo.json の追加（DOI 取得のタイミングで。著者: Mao Miyoshi, Ryota Ohno）
- [ ] README の英語チュートリアル拡充（インストール〜Tab 2 スキャン〜CLI 実行の通し手順）
- [ ] 論文 Methods への URL + DOI 記載（論文改訂時）
- [ ] **`paper_summary.txt`をgit履歴からも完全に削除**（2026-07-06 決定。現状は今後の
      ツリーからのみ除外済み＝`git log`を遡れば初期コミットにまだ残っている。
      論文が正式出版されて本格的にリポジトリを公開するタイミングで、`git filter-repo`等で
      履歴ごと消してforce push。**事前に大野さんに連絡し、force push後は
      彼にも再クローンしてもらう必要がある**（コミットハッシュが全部変わるため）。
      DOI取得・CITATION.cff追加と同じタイミングでまとめてやるのが良い）

**やった方がよい（2026-07-06、三好さん指示）**:
- [ ] N形2次元マップ（θ_incl・φ_incl、論文Fig.5b相当）のGUI実装。
      詳細は「N形2次元マップの再構成方法」セクション参照

**任意（余裕があれば）**:
- [ ] legacy スクリプトのハードコードパス（`~/path/to/monomer/` など）を引数化
- [ ] Tab 2 プロットのクリックで 3D プレビュー選択（現状はセレクトボックス）
- [ ] tcal result.txt → CSV の変換スクリプト（Tab 7 にそのままドロップできるように）
- [ ] Step 1a の GUI 版（`contact.step1a_scan`）と CLI 版（`sweep.py`+`extract_init.py`）の一本化
      （現状は同アルゴリズムの独立実装が2つ。CLI 版を GUI から呼ぶ形が自然か）

---

## GitHub 公開手順 / Publication

1. `git init` → `gh repo create MaoMiyoshi/ono-csp --public`
2. `pyproject.toml` に DOI プレースホルダーを記述
3. Zenodo 連携: GitHub Release を作成 → Zenodo が自動で DOI 発行
4. 論文中: `"The source code is available at https://github.com/MaoMiyoshi/ono-csp (DOI: ...)"` と記載
5. README に引用情報・インストール方法・最低限のチュートリアルを記述

---

## 論文対応表 / Paper Section Mapping

> **2026-07-06 更新**: 当初ドラフト（Step 1a/1b/2/3/4a/4b/5 という論文構成そのままの粒度）から、
> 大野さんの実コード受領後に確認できた**実際のTab/スクリプト対応**に更新。信頼度も付記。

| 論文 §| 図 | 対応するTab・スクリプト | 主変数 | 確認状況 |
|------|----|----|------|------|
| 2.1 | Fig.2, S1, S2 | Tab 2（vdW: `step1a_scan`／DFT: `step1.py`） | α(theta), a, b | 確認済み |
| 2.2 | Fig.5, S5, S6 | Tab 3 paraサブタブ（`step2_para.py`） | z(=θ_incl), θ_twist系はTab3 twist | 確認済み。ただしφ_incl中間値(N形)を探すコードは未発見 |
| 2.3 | Fig.6, S7 | Tab 4 paraサブタブ（vdW: `interlayer_vdw_scan`／DFT: `step3_para.py`） | Ra,Rb(x,y), z | 確認済み。vdWマップGUI実装済み(2026-07-06) |
| 2.4 | Fig.7 | Tab 3 twist（`step2_twist.py`）→ Tab 4 twist（`step3_twist.py`） | θ_twist(A2), ΔzT(Rt) | 確認済み。ただし論文はa,b,ΔzT**3変数**最適化と書くが`step2_twist.py`はRt固定・a,b自動最適化のみ（要確認） |
| 2.5 | Fig.8 | 未特定（おそらくTab4 para の Rp≠0探索と同じ仕組み） | θ'_incl, φ'_incl | **未確認**（N形探索コード不在の問題と同根） |
| 2.6 | Fig.S9, Table1 | `legacy/ono_scripts/tot_energy/`（csp本体には未統合、方針確認中） | G1〜G5外挿 | 確認済み（スコープ判断は三好さんの回答待ち） |
| 2.7 | Fig.10 | `legacy/ono_scripts/XRD_pattern/`（粉末XRDシミュレーション） | R形偏差 | 特定済み（2026-07-06、下記参照） |
| 2.8 | Fig.11 | Tab 5（`tcal_csv/`） | J | 確認済み |
