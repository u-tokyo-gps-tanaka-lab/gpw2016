デジタルカーリング用カーリング思考プログラム「歩」
ソースコード補助ドキュメント
2016/8/21

1. 盤面表現について

先手 BLACK (1) 後手 WHITE (0)
 - 基本的に「先手」「後手」は各エンド内の先後を表す

エンドは +∞ ~ 0 の減算型
エンド内のターンは 15 ~ 0 の減算型

エンド得点は 先手の得点が + 後手の得点が -
累計相対得点も同じく 先手リードのとき + 後手リードのとき -


2. 物理演算について

物理演算は1次元の長い配列に前計算データを入れておくことで高速化している


3. 連続空間での行動について

探索内で「着手(MoveXY, vmove_t)」はヒットやドローなどのラベル + 引数 による離散型
「ショット(fMoveXY)」は実際の系での連続的な値の組 を表すことが多い

離散着手は「ヒット」「ドロー」などのラベル以外に以下の離散的な要素を持つ
ショット選択の自由度が 2 なので第1引数と第2引数によってショットを調整できることが求められる

 要素 回転方向(R, L, より適した方(未実装)), コンタクト性, 相対性,
      第1引数(0 ~ 15), 第2引数(0 ~ 15), 第3引数(0 ~ 255)

 各着手の引数の意味 (詳しくは realizer.hpp にて)
                 第1        第2
   - PASS       なし
   - DRAW       ドロー位置
   - PREGUARD   ドロー位置
   - L1DRAW     なし
   - FREEZE     狙う石       相対位置
   - HIT        狙う石       初速
   - COMEAROUND 狙う石       相対長さ
   - POSTGUARD  狙う石       相対長さ
   - DOUBLE     前方の石     後方の石
   - PEEL       狙う石       狙う石を優先 <-> 両方出すことを優先

コンタクト着手の場合は 最初にコンタクトを狙う石を第1引数, それに付随して動かす石を第2引数に設定する決まりとしている