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
- 分散補正は **D3 (zero damping) = `gd3`**。論文 METHOD の引用文献69が
  Grimme のオリジナル D3 論文 (JCP 2010, 132, 154104) であり、BJ damping 論文
  (JCC 2011) は引用されていないため。※大野さんの実際の入力ファイルで要最終確認

**対象分子（プリセット）**:
- naphthalene (n=2)、anthracene (n=3)、tetracene (n=4)、pentacene (n=5)、hexacene (n=6)
- カスタム分子: ユーザーが XYZ ファイルをアップロード

---

## 計算パイプライン / Calculation Pipeline

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
  未ドロップ時は `example/pentacene/` の大野さんによる事前計算サンプルを表示（サンプルファイルは後日受領）。
  ユーザーがファイルをドロップしたらそちらのグラフ・3D モデル表示に切り替わる
- 機能は最小限に：「基本の流れが分かりやすくパッケージ化されている」ことを優先（ミーティングでの先生方針）

以降のワイヤーフレーム中の `[Run DFT Jobs]` 等の DFT 実行ボタンはこの方針により
「CLI コマンド紹介 + 結果 drag & drop」に置き換えて実装する。

---

## GUI タブ構成 / Streamlit UI Layout

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

## GitHub 公開手順 / Publication

1. `git init` → `gh repo create MaoMiyoshi/ono-csp --public`
2. `pyproject.toml` に DOI プレースホルダーを記述
3. Zenodo 連携: GitHub Release を作成 → Zenodo が自動で DOI 発行
4. 論文中: `"The source code is available at https://github.com/MaoMiyoshi/ono-csp (DOI: ...)"` と記載
5. README に引用情報・インストール方法・最低限のチュートリアルを記述

---

## 論文対応表 / Paper Section Mapping

| 論文 §| 図 | Package Step | 主変数 |
|------|----|----|------|
| 2.1 | Fig.2, S1, S2 | Step 1a + 1b | α, a, b |
| 2.2 | Fig.5, S5, S6 | Step 2 | θ_incl, φ_incl |
| 2.3 | Fig.6, S7 | Step 3 | x, y, z |
| 2.4 | Fig.7 | Step 4a | θ_twist |
| 2.5 | Fig.8 | Step 4b | θ'_incl, φ'_incl |
| 2.6 | Fig.S9, Table1 | (除外) | G1〜G5外挿 |
| 2.7 | Fig.10 | (参考) | R形偏差 |
| 2.8 | Fig.11 | Step 5 | J |
