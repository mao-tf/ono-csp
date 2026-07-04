# legacy/

大野さんがこれまで使ってきたCUIベースのPythonスクリプト置き場です（2026-07-04 受領）。

- `ono_scripts/stepwise_optimization/` — Step 1〜3 の山登り法最適化一式
- `ono_scripts/tcal_csv/` — transfer integral 一括計算（松井研 tcal ラッパー）

読解結果と `src/csp/` への統合状況は `spec.md`「大野コード対応表」を参照してください。
幾何・vdWスキャン・gjf生成・counterpoiseログ解析は移植済みです。
DFT 実行制御（山登りループ・ジョブ投入）は当面ここのスクリプトを直接使う想定のため、
このディレクトリは参照実装として残します（パスやジョブスケジューラ設定は各自の環境に
合わせて編集が必要です — 各タブの "How to run" 参照）。
