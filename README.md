Game Programming Workshop 2016 投稿論文
「カーリングAIに対するモンテカルロ木探索の適用」
実験コード

Katsuki Ohto 2016/10/6

使用方法：基本的にunix系の環境を想定しています

1. 物理演算ライブラリのbox2dをインストール
   macならhomebrewにて brew install box2d でインストールできる
   実験に用いたバージョンは 2.3.0 または 2.3.1

2. projects/curling 下へ移動

3. Makefile を書き換える、特にbox2d周り
   macのhomebrewにてインストールした場合には LIBRARIES に -lbox2d があればOK

4. make release

5. トップから research/curling/ 内に移動

6. research/curling 内の exp_full.sh を実行
   ただしシード設定等の手順が異なるので論文と異なる実験結果になると思われます
   こちらの実験で論文に使用した棋譜