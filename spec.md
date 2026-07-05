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
- [ ] 論文図のサンプル結果 CSV を `example/pentacene/` に同梱（大野さんから受領後。
      配置方法は `example/pentacene/README.md` 参照。後回し決定 2026-07-05）

**受領・回答待ち（大野さん）**:
- [x] ~~モノマー CSV の座標規約照合~~ → 2026-07-06 実CSV受領・数値比較して確認済み
      （X・Y完全一致、Z符号反転だが分子対称性により無害。「大野コード対応表」参照）
- [ ] tcal で pbepbe を使っていた経緯の確認（論文の図が pbepbe 由来なら
      2026-07-05 の b3lyp/6-31g* 修正を戻す必要がある）
- [x] ~~`step2_para.py` の `para_list.append(z,E1,E2)`~~ → 2026-07-05 修正済み
- [x] ~~Tab 3 para の z シフトが論文のどの図に対応するか~~ → 2026-07-06 論文§2.2本文で確認済み
      （「タブ別詳細ガイド」Tab 3 参照。z = θ_inclの物理的実体、SI式で変換可能）
- [x] ~~`step3_para.py` の `Rt`・`Rp` が論文§2.2のφ_incl中間値（N形）に対応するか~~
      → 2026-07-06 コードで確認済み（Rp = Rt1-Rt2、「タブ別詳細ガイド」Tab 3 参照）。
      ただし Rt,Rp は自動最適化されず手動スキャン設計である点は大野さんに確認したい
      （山登り対象に含めない設計意図か、単に未実装か）
- [ ] **新規**: `legacy/ono_scripts/tot_energy/`（論文§2.6の遠距離相互作用外挿 G1〜G5）を
      csp本体に統合するか、CLAUDE.mdの当初方針通りスコープ外のままにするか、三好さんに確認
      （2026-07-06 質問したが応答待ち。いったんスコープ外のまま進めている）
- [ ] **重要・新規（2026-07-06）**: 受領コード一式を`phi_incl`,`theta_incl`,`N_form`,`Rp`等で
      横断検索したが、N形を計算で発見・最適化するロジックは見当たらなかった
      （`step2_para/twist.py`はグライド対称方向のみ、`step3_para.py`・`tot_energy.py`は
      Rt,Rpを外から与えられた固定値として使うだけ）。**大野さんに直接聞くべき質問**:
      「論文Fig.5bのθ_incl・φ_incl 2次元マップ（N形の発見）はどのコードで計算しましたか？」
      未受領の別スクリプトがあるのか、Rt,Rpを手動グリッドで振って`step3_para.py`相当のもので
      計算したのかを確認する

**公開前の作業**:
- [ ] Zenodo 連携 → GitHub Release → DOI 取得（三好さんが対応。著者表記は Ohno で確定済み）
- [ ] CITATION.cff / .zenodo.json の追加（DOI 取得のタイミングで。著者: Mao Miyoshi, Ryota Ohno）
- [ ] README の英語チュートリアル拡充（インストール〜Tab 2 スキャン〜CLI 実行の通し手順）
- [ ] 論文 Methods への URL + DOI 記載（論文改訂時）

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
