#ifndef UTIL_HASHBOOK_HPP_
#define UTIL_HASHBOOK_HPP_

// ハッシュによる置換表
// 用途により様々なタイプを用意しておくつもり

#include "../analyze/analyzer.hpp"

class TwoValuePage32{
public:
    uint32 any()const noexcept{ return data_; }
    void setHash(const uint64 hash)noexcept{ data_ = (uint32)(hash >> 32); }
    void setResult(const uint32 res)noexcept{
        assert(!(res & (~3U)));
        data_ = (data_ & (~3U)) | res;
    }
    uint32 cmpHash(const uint64 hash)const noexcept{ return ((data_ ^ ((uint32)(hash >> 32))) & (~3U)); }
    uint32 getResult()const noexcept{ return (data_ & (3U)); }
    uint32 data()const noexcept{ return data_; }
    
private:
    std::atomic<uint32> data_;
};

class TwoValuePage64{
public:
    uint64 any()const noexcept{ return data_; }
    void setHash(const uint64 hash)noexcept{ data_ = hash; }
    void setResult(const uint64 res)noexcept{
        assert(!(res & (~3ULL)));
        data_ = (data_ & (~3ULL)) | res;
    }
    uint64 cmpHash(const uint64 hash)const noexcept{return ((data_ ^ (hash)) & (~3ULL)); }
    uint64 getResult()const noexcept{ return (data_ & (3U)); }
    uint64 data()const noexcept{ return data_; }
    
private:
    std::atomic<uint64> data_;
};

template<int SIZE>
class TwoValueBook{
    // 2(+中間1)値を保存するハッシュ表
    // 値1...1
    // 値2...2
    // 中間値...0
    // 未登録...-1
public:
    using page_t = TwoValuePage32;
    
    void init(){
        memset(page_, 0, sizeof(page_));
    }
    
    TwoValueBook(){
        init();
    }
    
    TwoValueBook(const std::string& argName)
    :ana(argName, SIZE, sizeof(*this))
    {
        init();
    }
    
    ~TwoValueBook(){}
    
    int read(const uint64 hash)noexcept{
        const page_t& fpage = page_[convHash_Index(hash)];
        if(fpage.any()){ // 結果が記入されている
            if(fpage.cmpHash(hash)){ // Hash値が違う
                ana.addUnfounded();
                return -1;
            }else{
                ana.addHit();
                return fpage.getResult();
            }
        }
        ana.addWhite();
        return -1;
    }
    
    void regist(const int result, const uint64 hash)noexcept{
        page_t& fpage = page_[convHash_Index(hash)];
        if(!fpage.any()){
            ana.addFilled();
        }
        fpage.setHash(hash);
        fpage.setResult(result);
        ana.addRegistration();
    }
private:
    
    HashBookAnalyzer ana;
    
    page_t page_[SIZE];
    
    constexpr static int convHash_Index(const uint64 hash)noexcept{
        return hash % (unsigned int)SIZE;
    }
};

template<int N_BITS>
class BitPage{
    // Nビットを保存する形式
    // (64 - N)ビットをハッシュ値の確認として利用
public:
    constexpr static uint64 MASK = (1ULL << N_BITS) - 1ULL;
    
    constexpr uint64 any()const noexcept{ return data_; }
    void set(const uint64 arg)noexcept{ data_ = arg; }
    constexpr uint64 cmpHash(const uint64 hash)const noexcept{
        return (data_ ^ hash) & (~MASK);
    }
    constexpr uint64 data()const noexcept{ return data_; }
    
private:
    std::atomic<uint64> data_;
};

template<int N_BITS, int SIZE>
class BitBook{
    
private:
    
    HashBookAnalyzer ana;
    
    BitPage<N_BITS> page_[SIZE];
    
    constexpr static int convHash_Index(const uint64 hash)noexcept{
        return hash % (unsigned int)SIZE;
    }
public:
    
    void init(){
        memset(page_, 0, sizeof(BitPage<N_BITS>) * SIZE);
    }
    
    BitBook(){
        init();
    }
    
    BitBook(const std::string& argName)
    :ana(argName, SIZE, sizeof(BitBook<N_BITS, SIZE>))
    {
        init();
    }
    
    ~BitBook(){}
    
    uint64 read(const uint64 hash)noexcept{
        const BitPage<N_BITS>& fpage = page_[convHash_Index(hash)];
        if(fpage.any()){//結果が記入されている
            if(fpage.cmpHash(hash)){//Hash値が違う
                ana.addUnfounded();
                return 0ULL;
            }else{
                ana.addHit();
                return fpage.data() & (BitPage<N_BITS>::MASK);
            }
        }
        ana.addWhite();
        return 0ULL;
    }
    
    void regist(const uint64 bits,const uint64 hash)noexcept{
        BitPage<N_BITS>& fpage = page_[convHash_Index(hash)];
        if(!fpage.any()){
            ana.addFilled();
        }
        fpage.set((hash & (~BitPage<N_BITS>::MASK)) | (bits & BitPage<N_BITS>::MASK));
        ana.addRegistration();
    }
};

template<class _data_t>
struct HashData{
    
    using data_t = _data_t;
    
    uint64 hash_;
    data_t data_;
    
    uint64 hash()const noexcept{ return hash_; }
    uint64 any()const noexcept{ return hash_; }
    void setData(const data_t& arg)noexcept{ data_ = arg; }
    void setHash(const uint64 ahash)noexcept{ hash_ = ahash; }
    bool cmp(const HashData<data_t>& arg)const noexcept{
        return (hash_ != arg.hash());
    }
    const data_t& data()const noexcept{ return data_; }
};

template<class page_t, int SIZE>
class HashBook{
    
private:
    HashBookAnalyzer ana;
    
    page_t page_[SIZE];
    
    constexpr static int convHash_Index(const uint64 hash)noexcept{
        return hash % (unsigned int)SIZE;
    }
public:
    
    void init(){
        memset(page_, 0, sizeof(page_));
    }
    
    HashBook(const std::string& argName):
    ana(argName, SIZE, sizeof(*this))
    {
        init();
    }
    
    ~HashBook(){}
    
    typename page_t::data_t read(const page_t& arg)noexcept{
        const page_t& fpage = page_[convHash_Index(arg.hash())];
        if(fpage.any()){ // 結果が記入されている
            if(fpage.cmp(arg)){ // 違うデータ
                ana.addUnfounded();
                return page_t::data_t(0);
            }else{
                ana.addHit();
                return fpage.data();
            }
        }
        ana.addWhite();
        page_t::data_t(0);
    }
    
    const page_t& page(const uint64 ahash)noexcept{
        return page_[convHash_Index(ahash)];
    }
    
    void regist(const typename page_t::data_t& data, const uint64 hash)noexcept{
        page_t& fpage = page_[convHash_Index(hash)];
        if(!fpage.any()){
            ana.addFilled();
        }
        fpage.setHash(hash);
        fpage.setData(data);
        ana.addRegistration();
    }
};

#endif // UTIL_HASHBOOK_HPP_