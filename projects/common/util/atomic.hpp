/*
atomic.hpp
Katsuki Ohto
*/

#ifndef UTIL_ATOMIC_HPP_
#define UTIL_ATOMIC_HPP_

// データ構造のatomic版

#include <iostream>
#include <atomic>

#include "../defines.h"

#include "bitArray.hpp"
#include "bitSet.hpp"

template<int N = 2, int SIZE = 8 / N>
using AtomicBitArray8 = BitArrayInRegister<uint8, 8, N, SIZE, std::atomic<uint8>>;

template<int N = 4, int SIZE = 16 / N>
using AtomicBitArray16 = BitArrayInRegister<uint16, 16, N, SIZE, std::atomic<uint16>>;

template<int N = 4, int SIZE = 32 / N>
using AtomicBitArray32 = BitArrayInRegister<uint32, 32, N, SIZE, std::atomic<uint32>>;

template<int N = 4, int SIZE = 64 / N>
using AtomicBitArray64 = BitArrayInRegister<uint64, 64, N, SIZE, std::atomic<uint64>>;

using AtomicBitSet8 = BitSetInRegister<uint8, 8, std::atomic<uint8>>;
using AtomicBitSet16 = BitSetInRegister<uint16, 16, std::atomic<uint16>>;
using AtomicBitSet32 = BitSetInRegister<uint32, 32, std::atomic<uint32>>;
using AtomicBitSet64 = BitSetInRegister<uint64, 64, std::atomic<uint64>>;

#endif // UTIL_ATOMIC_HPP_