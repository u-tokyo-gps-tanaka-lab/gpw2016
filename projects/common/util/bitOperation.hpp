/*
 bitOperation.hpp
 Katsuki Ohto
 */

#ifndef UTIL_BIT_OPERATION_HPP_
#define UTIL_BIT_OPERATION_HPP_

#include <cassert>

#include "../defines.h"

// ビット演算ユーティリティ
// 環境（プロセッサ、コンパイラ）依存する命令を使うものと、使わないものをいずれも用意したいところだが
// 現在はSSE4.?以降かつgccのみ対応

/**************************基本演算**************************/

// 包括性
template<typename T>
static inline constexpr bool holdsBits(const T a, const T b)noexcept{ return (((~a) & b) == T(0)); }

// 排他性
template<typename T>
static inline constexpr bool isExclusive(const T a, const T b)noexcept{ return ((a & b) == T(0)); }

/**************************ビット数判定**************************/

template<class T>
static constexpr T any2Bits(const T a)noexcept{
    return a & (a - static_cast<T>(1));
}

/**************************ビット数計算**************************/

// 環境依存なし

//畳み込み
static inline uint32 countManyBits64(uint64 a)noexcept{
    a = (a & 0x5555555555555555) + ((a >> 1) & 0x5555555555555555);
    a = (a & 0x3333333333333333) + ((a >> 2) & 0x3333333333333333);
    a = (a & 0x0f0f0f0f0f0f0f0f) + ((a >> 4) & 0x0f0f0f0f0f0f0f0f);
    a = (a & 0x00ff00ff00ff00ff) + ((a >> 8) & 0x00ff00ff00ff00ff);
    a = (a & 0x0000ffff0000ffff) + ((a >> 16) & 0x0000ffff0000ffff);
    return (int)((a & 0x00000000ffffffff) + ((a >> 32) & 0x00000000ffffffff));
}

static inline uint32 countManyBits32(uint32 a)noexcept{
    a = (a & 0x55555555) + ((a >> 1) & 0x55555555);
    a = (a & 0x33333333) + ((a >> 2) & 0x33333333);
    a = (a & 0x0f0f0f0f) + ((a >> 4) & 0x0f0f0f0f);
    a = (a & 0x00ff00ff) + ((a >> 8) & 0x00ff00ff);
    return (int)((a & 0x0000ffff) + ((a >> 16) & 0x0000ffff));
}

// 1ビットずつ数える
static inline uint32 countFewBits32(uint32 a)noexcept{
    int count = 0;
    for (; a; a &= (a - 1)){ ++count; }
    return count;
}

static inline uint32 countFewBits64(uint64 a)noexcept{
    int count = 0;
    for (; a; a &= (a - 1)){ ++count; }
    return count;
}

// 環境依存あり

#if defined(_MSC_VER)

static inline uint32 countBits32(uint32 a)noexcept{
    return __popcnt(a);
}
static inline uint32 countBits64(uint64 a)noexcept{
    return __popcnt64(a);
}

#elif defined(__GNUC__) && ( defined(__i386__) || defined(__x86_64__) )

static inline uint32 countBits32(uint32 a)noexcept{
    return __builtin_popcount(a);
}
static inline uint32 countBits64(uint64 a)noexcept{
    return __builtin_popcountll(a);
}

#else

static inline uint32 countBits32(uint32 a)noexcept{
    return countManyBits32(a);
}
static inline uint32 countBits64(uint64 a)noexcept{
    return countManyBits64(a);
}

#endif

template<typename T>inline uint32 countBits(T a)noexcept;

template<>inline uint32 countBits<uint8>(uint8 a)noexcept{ return countBits32(a); }
template<>inline uint32 countBits<uint16>(uint16 a)noexcept{ return countBits32(a); }
template<>inline uint32 countBits<uint32>(uint32 a)noexcept{ return countBits32(a); }
template<>inline uint32 countBits<uint64>(uint64 a)noexcept{ return countBits64(a); }


/**************************4ビットごと演算**************************/

// 部分ごとのビット数
template<typename T>static T countAll4Bits(T a)noexcept; // 4ビットごとのビット数

template<>uint64 countAll4Bits(uint64 a)noexcept{
    a = (a & 0x5555555555555555) + ((a >> 1) & 0x5555555555555555);
    return (a & 0x3333333333333333) + ((a >> 2) & 0x3333333333333333);
}

// 部分ごとのビット存在
template<typename T>static T gatherAll4Bits(T a)noexcept; // 4ビットごとのビット存在

template<>uint64 gatherAll4Bits(uint64 a)noexcept{
    a = (a & 0x5555555555555555) | ((a >> 1) & 0x5555555555555555);
    return (a & 0x3333333333333333) | ((a >> 2) & 0x3333333333333333);
}

// 部分ごとの全ビット存在
template<typename T>static T andAll4Bits(T a)noexcept; // 4ビットごとのビット存在

template<>
uint64 andAll4Bits(uint64 a)noexcept{
    a = (a & 0x5555555555555555)&((a >> 1) & 0x5555555555555555);
    a = (a & 0x3333333333333333)&((a >> 2) & 0x3333333333333333);
    assert(!(a & 0xEEEEEEEEEEEEEEEE));
    return a;
}

// 部分ごとに存在するビットを埋め尽くす
template<typename T>static T fillAll4Bits(T a)noexcept; // 4ビットごとの埋め尽くし

template<>uint64 fillAll4Bits(uint64 a)noexcept{
    a = ((a & 0x5555555555555555) << 1) | ((a >> 1) & 0x5555555555555555);
    a = ((a & 0x3333333333333333) << 2) | ((a >> 2) & 0x3333333333333333);
    assert(((uint32)countBits(a)) % 4 == 0);
    return a;
}

// 部分ごとに存在するビット枚数に対応した位置にのみビットを立てる
template<typename T>static T rankAll4Bits(const T& arg)noexcept;

template<>uint64 rankAll4Bits(const uint64& arg)noexcept{
    // 2ビットごとの枚数を計算
    uint64 a = (arg & 0x5555555555555555) + ((arg >> 1) & 0x5555555555555555);
    // 4ビットあったところを3に配置
    uint64 r = (a & 0x8888888888888888) & ((a << 2) & 0x8888888888888888);
    // 3ビットあったところを2に配置
    uint64 r3 = ((a << 2) & 0x4444444444444444) & ((a >> 1) & 0x4444444444444444);
    r3 |= ((a << 1) & 0x4444444444444444) & (a & 0x4444444444444444);

    // 残りは足すだけ。ただし3,4ビットがすでにあったところにはビットを置かない。
    uint64 r12 = (((a & 0x3333333333333333) + ((a >> 2) & 0x3333333333333333))) & 0x3333333333333333;
    if (r3){
        r |= r3;
        r |= (~((r3 >> 1) | (r3 >> 2))) & r12;
    }
    else{
        r |= r12;
    }

    assert(countBits(r) == countBits(gatherAll4Bits(arg)));

    return r;
}

/**************************ビット位置**************************/

// 名前と処理の関係がややこしいので注意
// ntz = ctz = bsf 最も下位のビットが下から何番目か
// clz 最も上位のビットが上から何番目か
// bsr 最も上位のビットが下から何番目か

// 環境依存無し
constexpr int NTZ64Table[64] = {
    0, 1, 59, 2, 60, 40, 54, 3,
    61, 32, 49, 41, 55, 19, 35, 4,
    62, 52, 30, 33, 50, 12, 14, 42,
    56, 16, 27, 20, 36, 23, 44, 5,
    63, 58, 39, 53, 31, 48, 18, 34,
    51, 29, 11, 13, 15, 26, 22, 43,
    57, 38, 47, 17, 28, 10, 25, 21,
    37, 46, 9, 24, 45, 8, 7, 6,
};

static constexpr inline int ntz64(const uint64 a)noexcept{
    return NTZ64Table[(int)(((a & (-a)) * 0x03F566ED27179461ULL) >> 58)];
}


//環境依存あり

#if defined(_MSC_VER)

static inline int bsf32(uint32 a)noexcept{
    unsigned long i;
    _BitScanForward(&i, a);
    return i;
}

static inline int bsf64(uint64 a)noexcept{
    unsigned long i;
    _BitScanForward64(&i, a);
    return i;
}

static inline int bsr32(uint32 a)noexcept{
    unsigned long i;
    _BitScanReverse(&i, a);
    return i;
}

static inline int bsr64(uint64 a)noexcept{
    unsigned long i;
    _BitScanReverse64(&i, a);
    return i;
}

static inline int ctz32(uint32 a)noexcept{
    return bsf32(a);
}

static inline int ctz64(uint64 a)noexcept{
    return bsf64(a);
}

static inline int clz32(uint32 a)noexcept{
    return 31 - bsr32(a);
}

static inline int clz64(uint64 a)noexcept{
    return 63 - bsr64(a);
}

#elif defined(__GNUC__) && ( defined(__i386__) || defined(__x86_64__) )

static inline int bsf32(uint32 a)noexcept{
    return __builtin_ctz(a);
}

static inline int bsf64(uint64 a)noexcept{
    return __builtin_ctzll(a);
}

static inline int bsr32(uint32 a)noexcept{
    int r;
    __asm__("bsrl %1, %0;" :"=r"(r) : "r"(a));
    return r;
}

static inline int bsr64(uint64 a)noexcept{
    int64 r;
    __asm__("bsrq %1, %0;" :"=r"(r) : "r"(a));
    return (int)r;
}

static inline int ctz32(uint32 a)noexcept{
    return __builtin_ctzl(a);
}

static inline int ctz64(uint64 a)noexcept{
    return __builtin_ctzll(a);
}

static inline int clz32(uint32 a)noexcept{
    return __builtin_clzl(a);
}

static inline int clz64(uint64 a)noexcept{
    return __builtin_clzll(a);
}

#else

static_assert(0, "no bsf-bsr.");

#endif

template<typename T>inline uint32 bsf(T a)noexcept;
template<typename T>inline uint32 bsr(T a)noexcept;

template<>inline uint32 bsf<uint8>(uint8 a)noexcept{ return bsf32(a); }
template<>inline uint32 bsf<uint16>(uint16 a)noexcept{ return bsf32(a); }
template<>inline uint32 bsf<uint32>(uint32 a)noexcept{ return bsf32(a); }
template<>inline uint32 bsf<uint64>(uint64 a)noexcept{ return bsf64(a); }
template<>inline uint32 bsr<uint8>(uint8 a)noexcept{ return bsr32(a); }
template<>inline uint32 bsr<uint16>(uint16 a)noexcept{ return bsr32(a); }
template<>inline uint32 bsr<uint32>(uint32 a)noexcept{ return bsr32(a); }
template<>inline uint32 bsr<uint64>(uint64 a)noexcept{ return bsr64(a); }

/**************************ビット取り出し、性質ビット抽出**************************/

//上位ビットの取り出し

template<typename T>
static inline T highestBit(const T a)noexcept;

template<>inline int32 highestBit(const int32 a)noexcept{
    return 1 << bsr32(a);
}
template<>inline uint32 highestBit(const uint32 a)noexcept{
    return 1U << bsr32(a);
}
template<>inline int64 highestBit(const int64 a)noexcept{
    return (int64)(1ULL << bsr64(a)); // 論理シフト
}
template<>inline uint64 highestBit(const uint64 a)noexcept{
    return 1ULL << bsr64(a);
}

template<typename T>
static T highestNBits(T a, int n)noexcept{
    assert(n > 0);
    // n <= 0には未対応
    // 流石に0では呼ばれないと思うが...
    T ans = static_cast<T>(0);
    T h;
    while (1){
        h = highestBit<T>(a);
        ans |= h;
        if (n == 1){ break; }
        --n;
        a -= h;
    }
    return ans;
}

// 下位ビット取り出し
template<typename T>
constexpr inline T lowestBit(const T a)noexcept{
    return (a & (-a));
}

template<typename T>
static T lowestNBits(T a, int n)noexcept{
    // n <= 0には未対応
    assert(n > 0);
    T ans = static_cast<T>(0);
    T l;
    while (1){
        l = lowestBit<T>(a);
        ans |= l;
        if (n == 1){ break; }
        --n;
        a -= l;
    }
    return ans;
}

static inline uint64 NthLowestBit64(uint64 a, int n)noexcept{
    uint64 bit = 0ULL;
    while (n){
        if (n == 1){ bit = lowestBit<uint64>(a); return bit; }
        a = a & (a - 1ULL);
        n--;
    }
    return bit;
}

template<typename T>
constexpr static T allLowerBits(T a)noexcept{
    // 最下位ビットより下位のビット全て
    return ((~a) & (a - static_cast<T>(1)));
}

template<typename T>
constexpr T allLowerBitsThan1Bit64(T a)noexcept{
    return a - static_cast<T>(1);
}

template<typename T>
static inline T allHigherBits(T a)noexcept;

template<>
inline uint64 allHigherBits<uint64>(uint64 a)noexcept{
    return 0xFFFFFFFFFFFFFFFFULL ^ ((1ULL << bsr64(a)) - 1ULL);
}

template<typename T>
constexpr static T allHigherBitsThan1Bit64(T a)noexcept{
    return (~((a - static_cast<T>(1)) << 1));
}

static inline uint64 NthHighestBit64(uint64 a, int n)noexcept{
    uint64 bit = 0ULL;
    while (n){
        assert(a);
        if (n == 1){ bit = highestBit<uint64>(a); return bit; }
        a &= ((1ULL << bsr64(a)) - 1ULL); // 最上位ビットだけを外す
        --n;
    }
    return bit;
}

static inline uint32 msb32(const uint32 a)noexcept{ return 1U << bsr32(a); }
static inline uint64 msb64(const uint64 a)noexcept{ return 1ULL << bsr64(a); }
static inline uint32 lsb32(const uint32 a)noexcept{ return (a & (-a)); }
static inline uint64 lsb64(const uint64 a)noexcept{ return (a & (-a)); }

static inline uint32 msb(const uint32 a)noexcept{ return msb32(a); }
static inline uint64 msb(const uint64 a)noexcept{ return msb64(a); }
static inline uint32 lsb(const uint32 a)noexcept{ return lsb32(a); }
static inline uint64 lsb(const uint64 a)noexcept{ return lsb64(a); }

/**************************ビットクロス**************************/



template<int N>
static uint64 genCrossNumber(const int p)noexcept{

    assert(p >= 0 && p<N);

    uint64 r;
    switch (N){
        case 0: UNREACHABLE; break;
        case 1: r = 0xFFFFFFFFFFFFFFFF << p; break;
        case 2: r = 0x5555555555555555 << p; break;
        case 3: r = 0x9249249249249249 << p; break;
        case 4: r = 0x1111111111111111 << p; break;
        case 5: r = 0x1084210842108421 << p; break;
        case 6: r = 0x4141414141414141 << p; break;
        case 7: r = 0x8102040810204081 << p; break;
        case 8: r = 0x0101010101010101 << p; break;
        default: UNREACHABLE; break;
    }
    return r;
}

constexpr static inline uint64 crossBits64(const uint64 a, const uint64 b)noexcept{
    return (a & 0x5555555555555555) | (b & 0xAAAAAAAAAAAAAAAA);
}

constexpr static inline uint64 crossBits64(const uint64 a, const uint64 b, const uint64 c)noexcept{
    return (a & 0x9249249249249249) | (b & 0x2492492492492492) | (c & 0x4924924924924924);
}

constexpr static inline uint64 crossBits64(const uint64 a, const uint64 b, const uint64 c, const uint64 d)noexcept{
    return (a & 0x1111111111111111) | (b & 0x2222222222222222) | (c & 0x4444444444444444) | (d & 0x8888888888888888);
}

constexpr static inline  uint64 crossBits64(const uint64 a, const uint64 b, const uint64 c, const uint64 d, const uint64 e)noexcept{
    return (a & 0x1084210842108421) | (b & 0x2108421084210842) | (c & 0x4210842108421084) | (d & 0x8421084210842108) | (e & 0x0842108421084210);
}

template<int N>static uint64 crossBits64(const uint64 a[]);

template<>uint64 crossBits64<2>(const uint64 a[]){
    return crossBits64(a[0], a[1]);
}

template<>uint64 crossBits64<3>(const uint64 a[]){
    return crossBits64(a[0], a[1], a[2]);
}

template<>uint64 crossBits64<4>(const uint64 a[]){
    return crossBits64(a[0], a[1], a[2], a[3]);
}

template<>uint64 crossBits64<5>(const uint64 a[]){
    return crossBits64(a[0], a[1], a[2], a[3], a[4]);
}

/**************************ビットマスク**************************/

template<typename T, int W, int MAX_SIZE = sizeof(T) * 8>T fillBits(T a)noexcept;
template<typename T, int W, int MAX_SIZE = sizeof(T) * 8>T fillBits_sub(T a, int n)noexcept;

template<typename T, int W, int MAX_SIZE>
T fillBits(T a)noexcept{
    return
        (W <= 0) ? static_cast<T>(0) :
        ((W >= MAX_SIZE) ? a :
        (((fillBits_sub<T, W, MAX_SIZE>(a, (MAX_SIZE - 1) / W)) << W) | a));
}

template<typename T, int W, int MAX_SIZE>
T fillBits_sub(T a, int n)noexcept{
    return
        (n <= 0) ? a :
        ((W <= 0) ? static_cast<T>(0) :
        ((W >= MAX_SIZE) ? a :
        ((fillBits_sub<T, W, MAX_SIZE>(a, n - 1) << W) | a)));
}

template<int W>
uint32 fillBits32(uint32 a)noexcept{
    return fillBits<uint32, W>(a);
}
template<int W>
uint64 fillBits64(uint64 a)noexcept{
    return fillBits<uint64, W>(a);
}

/**************************ビットを利用した基本演算**************************/

// 対数
template<typename T>static uint32 log2i(T a)noexcept;

template<>uint32 log2i(uint32 a)noexcept{
    return (uint32)bsr32(a);
}
template<>uint32 log2i(uint64 a)noexcept{
    return (uint32)bsr64(a);
}

// 2の累乗に切り上げ
static inline uint32 roundUpPow2(uint32 v)noexcept{ return msb(v - 1) << 1; } // 2の累乗に切り上げ
static inline uint64 roundUpPow2(uint64 v)noexcept{ return msb(v - 1) << 1; } // 2の累乗に切り上げ

/**************************ビットによる組み合わせ表現**************************/

// bit permutation。nビット立っているものを列挙する
static inline uint32 getNextBitPermutation(uint32 x)noexcept{
    int t = x | (x - 1);
    return (t + 1) | (uint32)((~t & -~t) - 1) >> (bsf(x) + 1);
}

// bit permutation。nビット立っているものを列挙する
static inline uint64 getNextBitPermutation(uint64 x)noexcept{
    uint64 t = x | (x - 1);
    return (t + 1) | (uint64)((~t & -~t) - 1) >> (bsf(x) + 1);
}

/**************************ビット入れ替え**************************/

template<class T>
T swapBits(T x, std::size_t p1, std::size_t p2, std::size_t len)noexcept{
    // x のp1からlen bitとp2からlen bitを入れ替え
    T mask = (static_cast<T>(1) << len) - 1;
    T ope = (x >> p1 ^ x >> p2) & mask;
    return x ^ (ope << p1 | ope << p2);
}

/**************************ビット逆順化**************************/

static inline uint32 reverseBits32(uint32 v)noexcept{
    v = (v >> 16) | (v << 16);
    v = ((v >> 8) & 0x00ff00ff) | ((v & 0x00ff00ff) << 8);
    v = ((v >> 4) & 0x0f0f0f0f) | ((v & 0x0f0f0f0f) << 4);
    v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    return ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
}

static inline uint64 reverseBits64(uint64 v)noexcept{
    v = (v >> 32) | (v << 32);
    v = ((v >> 16) & 0x0000ffff0000ffff) | ((v & 0x0000ffff0000ffff) << 16);
    v = ((v >> 8) & 0x00ff00ff00ff00ff) | ((v & 0x00ff00ff00ff00ff) << 8);
    v = ((v >> 4) & 0x0f0f0f0f0f0f0f0f) | ((v & 0x0f0f0f0f0f0f0f0f) << 4);
    v = ((v >> 2) & 0x3333333333333333) | ((v & 0x3333333333333333) << 2);
    return  ((v >> 1) & 0x5555555555555555) | ((v & 0x5555555555555555) << 1);
}

static inline uint32 reverseBits(uint32 v)noexcept{ return reverseBits32(v); }
static inline uint64 reverseBits(uint64 v)noexcept{ return reverseBits64(v); }

/**************************PEXT,PDEP**************************/

#if defined (HAVE_BMI2)
static inline uint32 pext32(uint32 a,uint32 msk)noexcept{
    return _pext_u32(a, msk);
}

static inline uint64 pext64(uint64 a, uint64 msk)noexcept{
    return _pext_u64(a, msk);
}

static inline uint32 pext(uint32 a, uint32 msk)noexcept{ return pext32(a, msk); }
static inline uint64 pext(uint64 a, uint64 msk)noexcept{ return pext64(a, msk); }

static inline uint32 pdep32(uint32 a, uint32 msk)noexcept{
    return _pdep_u32(a, msk);
}

static inline uint64 pdep64(uint64 a, uint64 msk)noexcept{
    return _pdep_u64(a, msk);
}

static inline uint32 pdep(uint32 a, uint32 msk)noexcept{ return pdep32(a, msk); }
static inline uint64 pdep(uint64 a, uint64 msk)noexcept{ return pdep64(a, msk); }
    
#else

static inline uint64 pext(uint64 bits, uint64 mask)noexcept{
    uint64 ans;
    if(mask){
        int a = bsf(mask);
        ans = (bits >> a) & 1;
        uint64 tmask = mask & (mask - 1ULL);
        if(tmask){
            int cnt = 1;
            while(1){
                int a = bsf(tmask);
                ans |= ((bits >> a) & 1) << cnt;
                tmask &= (tmask - 1ULL);
                if(!tmask)break;
                ++cnt;
            }
        }
    }else{
        ans = 0;
    }
    return ans;
}

static inline uint64 pdep(uint64 bits, uint64 mask)noexcept{
    uint64 ans;
    if(mask){
        int a = bsf(mask);
        ans = ((bits & 1ULL) << a);
        uint64 tmask = mask & (mask - 1ULL);
        uint64 tbits = bits >> 1;
        while(tmask){
            int a = bsf(tmask);
            ans |= ((tbits & 1ULL) << a);
            tmask &= (tmask - 1ULL);
            tbits >>= 1;
        }
    }else{
        ans = 0;
    }
    return ans;
}
    
#endif

#endif // UTIL_BIT_OPERATION_HPP_
