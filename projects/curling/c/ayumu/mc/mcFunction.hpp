/*
 mcFunction.hpp
 Katsuki Ohto
 */

#ifndef DCURLING_AYUMU_MC_MCFUNCTION_HPP_
#define DCURLING_AYUMU_MC_MCFUNCTION_HPP_

// 連続空間モンテカルロのための関数定義

namespace DigitalCurling{
    namespace Ayumu{
    
        double depthWeight(int d, int partNum){
            // ノードの深さによって掛ける重み
            // 深ければ深いほど大きくする
            return pow(d + 1, partNum);
        }
        
        double modifySize(double n, int partNum){
            // 重みつけ加算したトライアル回数を調整する
            // 調整しないと増えすぎる
            //return pow(n, 1 / (double)(partNum - 0.5));
            //return n;
#ifdef USE_DIVING
            return pow(n, 1 / (double)(2.5/* - 0.5*/));
#else
            return n;
#endif
        }
        
        template<class node_t>
        bool isExpansionCondition(const node_t& nd, const int range){
#ifdef MODE_SKY
            return range > RANGE_MIN_LEAFNODE && nd.size >= pow(1.3, RANGE_MAX_LEAFNODE - range);
#else
            return range > RANGE_MIN_LEAFNODE
            && nd.size >= N_DIV_TRIALS_LEAFNODE * pow(2, RANGE_MAX_LEAFNODE - range);
#endif
        }
    }
}

#endif // DCURLING_AYUMU_MC_MCFUNCTION_HPP_