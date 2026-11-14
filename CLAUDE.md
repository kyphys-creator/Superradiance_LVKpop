# Superradiance_LVKpop

> ⚠️ リモートは public(2026-11-10 のユーザー決定「public のままでいいです」)。未発表の研究なので、push するものは公開してよい内容に限る。メール・返信の下書きなどのやりとりは入れない。visibility を変えるのはユーザーの指示があるときだけ。

## 概要
LIGO–Virgo–KAGRA(LVK)の GWTC-5 の BH スピンの集団で、自己相互作用(四次結合、強さ $\propto 1/f^2$)のある超軽量スカラーの超放射を検定し、質量 $\mu$ と崩壊定数 $f$ の面 $(\mu, f)$ で排除または支持を決める研究プロジェクト。

- **理論**: 自己相互作用のある多準位の雲の運動論を解き、BH の最終スピンの写像 $\chi_f(M,\chi_i,\mu,f,t)$ を作る。
- **集団解析**: GWTC-5 の事後分布と感度の注入を使う。階層合体の 2 世代目の BH(生まれたときのスピン $\approx0.69$)を、自然スピンの仮定に頼らない較正に使う。
- **問いの軸**: 同じ GWTC-5 から、排除(arXiv:2607.01317、自由場)と証拠 $\ln B\approx7.8$(arXiv:2609.02678)という逆の結論が出ている。違いを自己相互作用と自然スピンの扱いで判定する。

テーマを選んだ経緯は `~/Documents/Claude/LVKtest/research/topic_selection_theory.md`(理論の強みから逆算した第 2 回の評価)。文献値の正本は LVKtest の文献レビュー(`~/Documents/Claude/LVKtest/Refs/lvk_physics_review/`、finding の ID は `T06-06` のように書く)。同じ GWTC-5 のデータを使う姉妹プロジェクトが `~/Documents/Claude/PBH_LVKpop/` にある。

## リポジトリ情報
- パス: `~/Documents/Claude/Superradiance_LVKpop/`
- ブランチ: `main`
- リモート: `kyphys-creator/Superradiance_LVKpop`(2026-11-07 に private で作成。2026-11-10 にユーザーの決定で public)
- コミットの署名: Mac Studio でのコミットはホスト名の署名。MacBook でのコミットは、E19(2026-11-14)から、その Mac の git の設定の署名(名前とメールアドレスが公開の履歴に載る。ユーザーの決定「この Mac の署名」)。
- 親の `~/Documents/Claude` のリポは NS-PBH のフォルダ以外を無視するので、ここは独立したリポジトリ。

## 構造
```
research/research_note.tex   研究ノート本体(日付つきエントリ entries/e0N_*.tex を \input)
research/entries/            エントリ(過去のエントリは書き換えず、訂正は新しいエントリで)
research/calc/               計算。エントリ N は calc/e0N_*.py、共通部品は calc/*.py、出力は calc/out/
research/figs/               pgfplots の図(データは figs/data/、calc が生成)
data/                        公開データ(本体は git に入れない。出どころ・取得手順・サイズは data/README.md)
Refs/notes/                  論文の本文からの抜き書き(式番号つき): population_papers、theory_rates、astro_inputs、window_consistency
Refs/arxiv/                  ダウンロードした論文のソース・PDF(git に入れない)
DESIGN.md                    決めた設計とその理由
paper/outline.md             論文の構成案と、データ解析の前に残る理論の作業
share/chi_sr/                集団解析用の χ_SR の単体版(共同研究者に渡す。表は research/calc/export_chisr.py が作る)
study/study_note.tex         勉強ノート(学部の物理から、この研究の物理を省略せずに導出。chapters/、figs/。ビルド: cd study && latexmk -lualatex study_note.tex)
```

## 実行環境
- 推論用の Python 環境は未作成(2026-11-03)。作るときはユーザーに確認する(候補: jax、numpyro、h5py、fsspec)。相対論的な超放射の率は自前の連分数法で計算するので SuperRad は要らない(DESIGN.md D5)。PBH_LVKpop と共用にするかもそのとき決める。
- 雲の運動論など numpy・scipy だけで済む計算は `~/miniforge3/bin/python3`(numpy 2.2、scipy 1.16)。
- ノート: `cd research && latexmk -lualatex research_note.tex`
- 計算機は Mac Studio(14 コア、36 GB)。重い並列計算の前に、1 本だけ走らせてメモリ(RSS)を測る。
- MacBook(Apple M3、8 コア、16 GB)では、base の `~/miniforge3/bin/python3` が numpy を読めない(libgfortran の rpath の重複)。`~/miniforge3/envs/python313/bin/python`(numpy 2.3、scipy 1.16、mpmath)を使う。E19 の計算はこれで行い、E18 の保存した出力を完全に再現する。
- E18 の基準の表 `research/calc/out/reftab.npz`(700 MB、git に入れない)は Mac Studio にだけある。ない計算機では `chisr.build()` が `share/chi_sr/chi_sr_tables.npz`(12 MB、git にある)から同じ数の表を読むので、χ_SR と `export_chisr.py` の確かめはどちらの計算機でも動く(E19)。値の表を作り直す(`e18_reftab_build.py`、約 5 時間)のは Mac Studio で。

## ルール
- 数値はすべて `research/calc/` のスクリプトから出す。文献の式・数値を使うときは、原文(arXiv の本文)で照合した式番号を書く。
- 先行研究との比較は、このノートの模型で計算した量で行う。文献の図や値を写して比べない。
- 論文(`paper/`、まだない)では、先行研究を人名で前に出して不足を書かない(`~/Documents/Claude/CONVENTIONS.md` §5.9)。コードやリポジトリにも触れない(§5.8)。
- データ・パッケージのダウンロードは、ファイル名・出どころ・サイズを示してユーザーの許可を得てから行う。

## How to Resume
1. SESSION.md を読む → 現在状態と残タスクを把握
2. DESIGN.md で決めたこと、研究ノートの最新エントリで結果を確認
3. 残タスクに従って作業継続。不明点はユーザーに確認

## 自動更新ルール(必須)
- タスク完了時 → SESSION.md 更新(完了タスクは push 時に除去)
- 重要な判断時 → SESSION.md と DESIGN.md に記録
- push 前 → SESSION.md 更新・棚卸し → 4軸レビュー → commit
- 詳細は `~/Documents/Claude/CONVENTIONS.md` §3 参照
