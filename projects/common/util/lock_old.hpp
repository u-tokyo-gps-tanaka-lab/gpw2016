/*
 lock.hpp
 Katsuki Ohto
 */


#ifndef UTIL_LOCK_HPP_
#define UTIL_LOCK_HPP_

#include "../defines.h"

//排他処理用スピンロック
//モンテカルロのコンピュータ囲碁の本を参考に

void lock(volatile int *ptr){
#if defined(_MSC_VER)
    __asm
    {
        mov  ecx, &l
    la: mov  eax, 1
        xchg eax, [ecx]
        test eax, eax
        jz   end
    ln:
        pause
        mov  exax, [ecx]
        test eax, eax
        jz   la
        jmp  lb
    end:
    }
#else
    int itmp;
    for(;;){
        asm( "1: movl  $1,  %1 \n\t"
            "   xchgl (%0),%1 \n\t"
            :"=g" (ptr), "=r" (itmp) : "0" (ptr) );
        if( !itmp ){return;}
        while(*ptr);
    }
#endif
    assert(0);
}
void unlock(volatile int *ptr){
    *ptr = 0;
}

class SpinLock{
public:
    void lock()noexcept{
        lockInner(&flag_);
    }
    void unlock()noexcept{
        flag_ = 0;
    }
private:
    std::atomic<int> flag_;
    
    void lockInner(std::atomic<int> *ptr)noexcept{
        ASSERT(ptr != nullptr,);
#if defined(_MSC_VER)
        __asm
        {
            mov  ecx, &l
        la: mov  eax, 1
            xchg eax, [ecx]
            test eax, eax
            jz   end
        ln:
            pause
            mov  exax, [ecx]
            test eax, eax
            jz   la
            jmp  lb
        end:
        }
#else
        int tmp;
        for(;;){
            asm("1: movl  $1,  %1 \n\t"
                "   xchgl (%0),%1 \n\t"
                :"=g" (ptr), "=r" (tmp) : "0" (ptr)  );
            if(!tmp){ return; }
            while(*ptr);
        }
#endif
        UNREACHABLE;
    }
};

#endif // UTIL_LOCK_HPP_