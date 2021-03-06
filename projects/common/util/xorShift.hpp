/*
 xorShift.hpp
 Katsuki Ohto
 */


#ifndef UTIL_XORSHIFT_HPP_
#define UTIL_XORSHIFT_HPP_

#include "../defines.h"

// xorshiftによる乱数生成

// クラス
class XorShift64{
private:
    uint64 x, y, z, t;
public:

    uint64 rand()noexcept{
        uint64 tmp = x ^ (x << 11);
        x = y;
        y = z;
        z = t;
        t = (t ^ (t >> 19)) ^ (tmp ^ (tmp >> 8));
        return t;
    }
    
    double drand()noexcept{
        uint64 tmp = x ^ (x << 11);
        x = y;
        y = z;
        z = t;
        t = (t ^ (t >> 19)) ^ (tmp ^ (tmp >> 8));
        return t / static_cast<double>(0xFFFFFFFFFFFFFFFFULL);
    }
    void srand(const uint64 s)noexcept{
        if (!s){ // seedが0だとまずい
            x = 0x0123456789ABCDEFULL;
        }
        else{
            x = (s << 32) ^ s;
        }
        y = (x << 8) | ((x & 0xff00000000000000ULL) >> 56);
        z = (y << 8) | ((y & 0xff00000000000000ULL) >> 56);
        t = (z << 8) | ((z & 0xff00000000000000ULL) >> 56);
    }
    constexpr static uint64 min()noexcept{ return 0ULL; }
    constexpr static uint64 max()noexcept{ return 0xFFFFFFFFFFFFFFFFULL; }

    constexpr XorShift64()noexcept
        :x(), y(), z(), t(){}

    XorShift64(const uint64 s)
        : x(), y(), z(), t(){
        srand(s);
    }
};


#endif // UTIL_XORSHIFT_HPP_