デジタルカーリング用カーリング思考プログラム
歩 (あゆむ, Ayumu)
大渡勝己 (Katsuki Ohto)
最終更新 2016/10/6

--- C++ Part ---

// ビルド準備
汎用物理演算エンジンのBox2D(バージョン2.3系)を使える環境が必須

macならHomebrew にて 
brew install box2d
としリンク時に
-lbox2d
をつければ良い

windowsなら何処かに.dllが落ちてあった
リンク時に
-llibbox2D
をつける

ソースからビルドしてリンクするのはどんな環境でもできるはず

// ビルド方法
本ディレクトリにて
make release 試合用
make default デバッグ用情報つきビルド
make debug 高速化なし, デバッグ出力(DERR)あり
のどれか

// バイナリ説明
client, policy_client TCP/IP接続する思考プログラム

usage: client -h HOST -p PORT -g GAMES

engine, policy_engine パイプ接続を仮定した思考プログラム

server, server_gat TCP/IP接続を待つサーバープログラム

usage: server -p PORT -g GAMES -e ENDS -t TIME_LIMIT(ms) -l RECORD_DIRECTORY

policy_learner 方策関数の学習
eval_learner 静的評価関数の学習

dc_test 基本演算テスト
error_test 着手にかかる外乱のチェック
variance_test
eval_test
simulator_test 高速物理演算システムの性能チェック
shot_test ダブルテイクなどのショット生成の性能チェック
interpolation_test 近い着手による着手後局面の線形補完の精度チェック

dcl_converter .dcl形式の棋譜を .dat形式のデータに変換
dcl_counter 棋譜から勝ち負けを集計(python版の方が手軽でよい)

table_generator 思考に利用するパラメータのファイル生成
image_log_maker NNでの学習に利用する盤面表現行列ファイルの生成

// コード概観
/projects/common/ 下
共通のユーティリティコード

/projects/curling 下
/c

  setings.h 名前, ビルドと戦略等設定

  client.cc 通信と思考部の呼び出しメイン

  dc.hpp 基本クラスと演算の定義

  /ai 現在未使用
  /ayumu 「歩」思考部
  /log デジタルカーリング公式クライアントが生成する .dcl 形式ファイルを扱う
  /server 公式コードベースのTCP/IP接続の対戦プログラム
  /simulation カーリングの物理演算のコード
  /solver 実験用で現在未使用
  /structure 特に「歩」の思考に限定しない汎用クラス
  /test 動作テストや実験, パラメータ生成のメイン

  /ayumu

    ayumu_dc.hpp 「歩」の演算のための基本クラス等の定義

    ayumu.hpp 歩の思考のメインクラス

    evlearner.cc 静的評価関数の学習

    initialize.hpp 初期化をまとめてやる

    pglearner.cc 方策の学習

    /db 探索中に過去の類似局面を探してその手を打つなどしようと思ったが現在未使用
    /eval エンド途中や末端局面の勝率を予測
    /learn 学習を回す
    /logic 判定
    /mc モンテカルロ法
    /model ニューラルネット系の学習のためのモデル作成
    /move 着手の生成
    /policy 方策
    /search 現在未使用
    /shot ショットの生成
    /structure 「歩」の思考に用いるクラス

  /log

    dcl_converter.cc .dcl 形式の棋譜の情報をまとめて .dat 系のデータに保存
    dcl_counter.cc .dcl 形式の棋譜を読んで勝ち負け等をカウント(Pythonの方が楽なので使わない)
    dclMaker.hpp サーバーが公式と同じ .dcl 形式のファイルを作成する用

  /server

    CurlingSimulator.cpp ほぼ公式の物理シミュレータ
    CurlingSimulator.h 公式由来
    Message.h 公式由来
    server.cpp 公式のコードを書き換えて作ったTCP/IP接続サーバーメイン
    System.h 公式由来

  /simulation

    b2dSimulator.h 汎用物理エンジンBox2Dを利用する物理演算関数の定義(公式由来の関数あり)
    b2dSimulator.hpp 汎用物理エンジンBox2Dを利用する物理演算関数の実装(公式由来の関数あり)
    collision.hpp 衝突処理(公式に衝突の式が定義されていないのでBox2Dの実装に近づけている)
    error.hpp 着手にかかる外乱の式(公式由来の関数あり)
    fastSimulator.h 高速シミュレータ関数定義
    fasutSimulator.hpp 高速シミュレータメイン
    fastSimulatorFunc.hpp 高速シミュレータ用のテーブルと演算
    intSimulator.h 整数計算メインのシミュレータを作ろうとして断念中
    intSimulator.hpp 整数計算メインのシミュレータを作ろうとして断念中
    primaryShot.hpp ドローやヒットなど基本的なショットの生成
    simuFunc.hpp 公式ベースの物理演算関数いろいろ
    simuRecord.hpp 現在未使用

  /solver

  /structure

    field.hpp 公式の値に比較的近いルートでの盤面表現
    grid.hpp 格子点から積分計算
    log.hpp 現在未使用

  /test

    ayumu_test.cc
    dc_test.cc
    error_test.cc
    eval_test.cc
    image_log_maker.cc
    interpolation_test.cc
    logic_test.cc
    reinforcement_test.cc
    shot_test.cc
    simulator_test.cc
    solver_test.cc
    table_generator.cc
    test150812.cc 昔書いた
    test150813.cc 昔書いた
    variance_test.cc

  /ayumu/db

  /ayumu/eval

    estimator.hpp 静的局面評価関数(17次元の得点予測ベクトル)を計算
    eval.hpp 現在未使用
    evaluator.hpp エンド末端評価と補助評価関数
    stat.hpp エンド得点分布や延長戦勝率のデータ

  /ayumu/learn

    classifier.hpp 棋譜のショットにラベルを付ける
    PG.hpp 方策学習を回す
    scorePG.hpp 局面評価関数学習を回す

  /ayumu/logic

    logic.hpp 相対位置のラベルや判定関数
    mate.hpp 詰み判定

  /ayumu/mc

    leaf.hpp プレイアウトのルート以外の処理
    mcFunction.hpp 連続状態内でモンテカルロを行うための関数
    mcPlayer.hpp モンテカルロ着手決定の入り口
    monteCarlo.h モンテカルロ木探索の定数や置換表等のクラス定義
    root.hpp プレイアウトのルートでの処理

  /ayumu/model

    image.hpp 局面画像生成
    rbf.hpp RBFネットワークを書こうとした

  /ayumu/move

    estimator.hpp ショット生成器が出すショットの成功予測への橋渡し
    generator.hpp 離散着手ラベル生成
    realizer.hpp 離散着手ラベルからショット生成器への橋渡し

  /ayumu/policy

    heuristic.hpp 手で組んだ方策
    L1.hpp 手で組んだラストストーン方策
    policy.hpp 線形行動価値関数によるソフトマックス方策の計算

  /ayumu/search

  /ayumu/shot

    allShots.hpp 全てのショットを一度にインクルードする
    backTake.hpp
    double.hpp ダブルテイクアウト 2つの石を一緒にテイクアウト
    drawRaise.hpp ドローレイズ 目的の石をハウス中心に押し込む
    frontTake.hpp
    l1Draw.hpp ラストストーンを空いたところにドロー
    param.h 主にショット生成に使うパラメータテーブル
    peel.hpp ピール ヒットして自分の石も一緒にプレーエリア外に出す
    raiseTakeOut.hpp 前の石を叩いて後ろの石に当てて後ろの石をテイクアウト
    shotGenerator.hpp ベースクラス(不要?)
    takeOut.hpp ガードを避けてテイクアウト

  /ayumu/structure

    stoneInfo.hpp 1つの石の情報を表すクラス
    thinkBoard.h 思考用局面表現ヘッダ
    thinkBoard.hpp 思考用局面表現実装
    thinkField.hpp ThinkBoard に時間制限などの情報を加えてルートの盤面情報として使用可能にしたもの


--- Python Part ---

client.py TCP/IP接続してサーバーにつなぐベースプログラム

udsage: python client.py