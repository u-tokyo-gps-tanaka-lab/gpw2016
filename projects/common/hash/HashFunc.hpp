/*
 hashfunc.hpp
 Katsuki Ohto
 */

#ifndef UTIL_HASHFUNC_HPP_
#define UTIL_HASHFUNC_HPP_

//ハッシュ値生成
//とりあえず生成関数を毎回呼べば動くだろうが、
//１からの計算は出来る限り避けるように

//新旧いずれのデータ構造も含む

//ビット演算として2つから5つの64ビット値をクロスさせる
//crossBits64を定義しておく必要あり

//関数名
//gen...情報より0から計算
//genproc...ハッシュ値に対して進行情報を与えて進行させる
//proc...ハッシュ値に対して進行ハッシュ値を与えて進行させる
//merge...同種の情報のハッシュ値を合成、進行
//knit...主に異種の情報のハッシュ値から求めるハッシュ値を計算

//基本関数

//線形加算
constexpr uint64 addHash(uint64 h0, uint64 h1)noexcept{ return h0 ^ h1; }

/**************************ビット**************************/

//1<<0,1<<1,...,1<<63の値を0~63に転置する完全ハッシュ関数
//1<<0,1<<1,...,1<<31の値を0~31に転置する完全ハッシュ関数
//としてどの関数を用いるか
//SSEでビットの位置を知るか、乗算を用いるかなど

//ただし、この関数を使わずにインデックスを得てもよいことになっており
//あくまで汎用ユーティリティという位置づけ

inline int genPerfectHash_Bit64(uint64 b)noexcept{
    //return (int)(((b)*0x03F566ED27179461ULL)>>58);//乗算
    return bsf64(b);
}

inline int genPerfectHash_Bit32(uint32 b)noexcept{
    return bsf32(b);
}

/**************************交叉ハッシュ**************************/

//交叉
constexpr uint64 crossHash(uint64 h0, uint64 h1)noexcept{
    return crossBits64(h0, h1);
}
constexpr uint64 crossHash(uint64 h0, uint64 h1, uint64 h2)noexcept{
    return crossBits64(h0, h1, h2);
}
constexpr uint64 crossHash(uint64 h0, uint64 h1, uint64 h2, uint64 h3)noexcept{
    return crossBits64(h0, h1, h2, h3);
}
constexpr uint64 crossHash(uint64 h0, uint64 h1, uint64 h2, uint64 h3, uint64 h4)noexcept{
    return crossBits64(h0, h1, h2, h3, h4);
}

template<int N>
uint64 crossHash(const uint64 h[]){
    return crossBits64<N>(h);
}

//部分生成
template<int N>
constexpr uint64 genPartCrossedHash(const int n, const uint64 hash)noexcept{
    return hash & genCrossNumber<N>(n);
}

template<int N>
constexpr uint64 genPart2CrossedHash(const int n, const uint64 hash0, const uint64 hash1)noexcept{
    return (hash0 ^ hash1) & genCrossNumber<N>(n);
}

//交叉更新
template<int N>
constexpr uint64 procCrossedHash(const uint64 hash, const int n, const uint64 hash_dist)noexcept{
    return hash ^ (hash_dist & genCrossNumber<N>(n));
}

#endif // UTIL_HASHFUNC_HPP_

