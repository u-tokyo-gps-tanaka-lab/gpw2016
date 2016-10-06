/*
 analyzer.hpp
 Katsuki Ohto
 */

#ifndef UTIL_ANALYZER_HPP_
#define UTIL_ANALYZER_HPP_

#include <iostream>

#include "../defines.h"

// メソッドアナライザー
// メソッドのcpu時間や使用回数など諸々を調べる

// part...どれかが単体で動くことを想定
// phase...一度の処理の中で連続して存在することを想定

// 簡単に使うため、アナライザクラスには解析可能性のあるデータについて全て実装しておく
// 最適化で邪魔な部分が消えてくれない限り、相当なメモリの無駄だが
// 本番用では全体をオフにするので不問とする

// std:atomicでスレッドセーフにしましたが少し遅くなっていると思います

namespace Analysis{
    // アナライザーの表示タイプ
    
    enum{
        TYPE_PLANE = 0, // 回数、時間計測、失敗回数のみ
        TYPE_SEARCH, // 探索結果
        TYPE_PLAYOUT, // 未実装
    };
    
}

// アナライザと呼ぶほどもでない静的カウンタ、時間計測器

// デストラクタで回数を表示するカウンタ
struct Counter{
#ifdef USE_ANALYZER
    uint64 c_;
    std::string nm_; // 表示名
    Counter() : c_(0), nm_("COUNTER"){}
    Counter(uint64 ac) : c_(ac), nm_("COUNTER"){}
    Counter(const std::string& nm) : c_(0), nm_(nm){}
    Counter(uint64 ac, const std::string& nm) : c_(ac), nm_(nm){}
    ~Counter(){ std::cerr << nm_ << " : " << c_ << std::endl; }
    operator uint64()const noexcept{ return c_; }
    Counter& operator=(uint64 n)noexcept{ c_ = n; return *this; }
    Counter& operator++()noexcept{ c_ += 1; return *this; }
    Counter& operator+=(uint64 n)noexcept{ c_ += n; return *this; }
#else
    Counter(){}
    Counter(uint64 ac){}
    Counter(const std::string& nm){}
    Counter(uint64 ac, const std::string& nm){}
    ~Counter(){}
    operator uint64()const noexcept{ return 0ULL; }
    Counter& operator=(uint64 n)noexcept{ return *this; }
    Counter& operator++()noexcept{ return *this; }
    Counter& operator+=(uint64 n)noexcept{ return *this; }
#endif
};

#ifdef USE_ANALYZER
std::ostream& operator<<(std::ostream& ost, const Counter& c){
    ost << c.c_;
    return ost;
}
#else
std::ostream& operator<<(std::ostream& ost, const Counter& c){
    ost << "none";
    return ost;
}
#endif

// アトミック演算でより正確なカウント
struct AtomicCounter{
#ifdef USE_ANALYZER
    std::atomic<uint64> c_;
    std::string nm_; // 表示名
    AtomicCounter() : c_(0), nm_("COUNTER"){}
    AtomicCounter(uint64 ac) : c_(ac), nm_("COUNTER"){}
    AtomicCounter(const std::string& nm) : c_(0), nm_(nm){}
    AtomicCounter(uint64 ac, const std::string& nm) : c_(ac), nm_(nm){}
    ~AtomicCounter(){ std::cerr << nm_ << " : " << c_ << std::endl; }
    uint64 count()const noexcept{ return c_; }
    void add(uint64 n = 1)noexcept{ c_ += n; }
    AtomicCounter& operator=(uint64 n)noexcept{ c_ = n; return *this; }
    AtomicCounter& operator++()noexcept{ c_ += 1; return *this; }
    AtomicCounter& operator+=(uint64 n)noexcept{ c_ += n; return *this; }
#else
    AtomicCounter(){}
    AtomicCounter(uint64 ac){}
    AtomicCounter(const std::string& nm){}
    AtomicCounter(uint64 ac, const std::string& nm){}
    ~AtomicCounter(){}
    static uint64 count()noexcept{ return 0ULL; }
    void add(uint64 n = 1)const noexcept{}
    AtomicCounter& operator=(uint64 n)noexcept{ return *this; }
    AtomicCounter& operator++()noexcept{ return *this; }
    AtomicCounter& operator+=(uint64 n)noexcept{ return *this; }
#endif
};

struct StaticClock{
#ifdef USE_ANALYZER
    Clock cl;
    std::atomic<uint64> c_, t_;
    StaticClock() : c_(0ULL), t_(0ULL){}
    StaticClock(uint64 ac, uint64 at) : c_(ac), t_(at){}
    ~StaticClock(){ report(); }
    void start()noexcept{ cl.start(); }
    void stop()noexcept{
        uint64 tmp = cl.stop();
        ++c_;
        t_ += tmp;
    }
    void report()const{
        std::cerr << "CLOCK : " << t_ << " clock ( " << c_ << " times)";
        if (c_ > 0){
            uint64 avg = t_ / c_;
            std::cerr << " " << avg << " per trial.";
        }
        std::cerr << std::endl;
    }
#else
    StaticClock(){}
    StaticClock(uint64 ac, uint64 at){}
    ~StaticClock(){}
    void start()const noexcept{}
    void stop()const noexcept{}
    void report()const noexcept{}
#endif
};


template<int PART_NUM = 1, int PHASE_MAX = 1, int MODE = 0>
struct StaticAnalyzer{
    //静的なやつ 手軽
#ifdef USE_ANALYZER
    Clock clock;
    const std::string& name;
    
    std::atomic<uint64> trials[PART_NUM];
    //時間
    std::atomic<uint64> time[PART_NUM][PHASE_MAX];
    //失敗
    std::atomic<uint64> failures[PART_NUM];
    //探索解析
    std::atomic<uint64> nodes[PART_NUM];
    std::atomic<uint64> childs[PART_NUM];
    //プレイアウト解析
    std::atomic<uint64> turns[PART_NUM];
#endif
    
    void startClock()noexcept{
#ifdef USE_ANALYZER
        clock.start();
#endif
    }
    
    uint64 stopClock()noexcept{
#ifdef USE_ANALYZER
        return clock.stop();
#else
        return 0ULL;
#endif
    }
    
    uint64 restartClock(){
#ifdef USE_ANALYZER
        return clock.restart();
#else
        return 0ULL;
#endif
    }
    
    void addTrial(const int pa = 0){
#ifdef USE_ANALYZER
        ++trials[pa];
#endif
    }
    
    void addTime(const uint64 argTime, const int pa = 0, const int ph = 0){
#ifdef USE_ANALYZER
        time[pa][ph] += argTime;
#endif
    }
    
    void addFailure(const int pa = 0){
#ifdef USE_ANALYZER
        ++failures[pa];
#endif
    }
    
    void addNodes(const int argNodes, const int pa = 0){
#ifdef USE_ANALYZER
        nodes[pa] += argNodes;
#endif
    }
    
    void addChilds(const int argChilds, const int pa = 0){
#ifdef USE_ANALYZER
        childs[pa] += argChilds;
#endif
    }
    
    void addTurns(const int argTurns, const int pa = 0){
#ifdef USE_ANALYZER
        turns[pa] += argTurns;
#endif
    }
    
    void start()noexcept{
        startClock();
    }
    
    void restart(uint64 argTime, int pa, int ph){ // 外部クロックの場合
        addTime(argTime, pa, ph);
    }
    
    void restart(int pa = 0, int ph = 0){
        addTime(stopClock(), pa, ph);
        startClock();
    }
    
    void end(uint64 argTime, int pa, int ph){ // 外部クロックの場合
        addTime(argTime, pa, ph);
        addTrial(pa);
    }
    
    void end(int pa = 0, int ph = 0){
        end(stopClock(), pa, ph);
    }
    
    void init()noexcept{
#ifdef USE_ANALYZER
        for (int i = 0; i < PART_NUM; ++i){
            trials[i] = 0U;
            failures[i] = 0U;
            for (int j = 0; j < PHASE_MAX; ++j){
                time[i][j] = 0ULL;
            }
        }
#endif
    }
    
    void report()const{
#ifdef USE_ANALYZER
        uint64 time_all_sum = 0ULL;
        uint64 time_part_sum[PART_NUM] = { 0ULL };
        uint64 time_per_trial_all = 0ULL;
        uint64 time_per_trial_part[PART_NUM] = { 0ULL };
        uint64 time_per_trial_phase[PART_NUM][PHASE_MAX] = { 0ULL };
        
        uint64 trial_sum = 0U;
        uint64 failure_sum = 0U;
        
        //sum
        for (auto i = 0; i < PART_NUM; ++i){
            if (trials[i]){
                trial_sum += trials[i];
                failure_sum += failures[i];
                for (int j = 0; j < PHASE_MAX; ++j){
                    time_part_sum[i] += time[i][j];
                    time_per_trial_phase[i][j] = time[i][j] / trials[i];
                }
                time_all_sum += time_part_sum[i];
                time_per_trial_part[i] = time_part_sum[i] / trials[i];
            }
        }
        if (trial_sum){
            time_per_trial_all = time_all_sum / trial_sum;
        }
        
        std::cerr << "****** Analysis of " << name << " ******" << std::endl;
        
        std::cerr << " < Normal Analysis >" << std::endl;
        
        std::cerr << " " << trial_sum << " trials.  " << time_all_sum << " clock ( " << time_per_trial_all << " clock/trial).  " << failure_sum << " failures." << std::endl;
        if (PART_NUM > 1 || PHASE_MAX > 1){
            for (auto i = 0; i < PART_NUM; ++i){
                std::cerr << "   Part " << i << " : " << trials[i] << " trials.  " << time_part_sum[i] << " clock ( " << time_per_trial_part[i] << " clock/trial).  " << failures[i] << " failures." << std::endl;
                if (PHASE_MAX > 1){
                    for (auto j = 0; j < PHASE_MAX; ++j){
                        std::cerr << "      Phase " << j << " : " << time[i][j] << " clock ( " << time_per_trial_phase[i][j] << " clock/trial)." << std::endl;
                    }
                }
            }
        }
        
        if (MODE == Analysis::TYPE_SEARCH){
            uint64 node_sum = 0U;
            uint64 child_sum = 0U;
            
            uint64 node_per_trial_part[PART_NUM] = { 0U };
            uint64 node_per_trial_all = 0U;
            
            uint64 time_per_node_part[PART_NUM] = { 0ULL };
            uint64 time_per_node_all = 0ULL;
            
            uint64 child_per_trial_part[PART_NUM] = { 0U };
            uint64 child_per_trial_all = 0U;
            
            uint64 time_per_child_part[PART_NUM] = { 0ULL };
            uint64 time_per_child_all = 0ULL;
            
            for (auto i = 0; i < PART_NUM; ++i){
                if (trials[i]){
                    node_sum += nodes[i];
                    child_sum += childs[i];
                    
                    if (nodes[i]){
                        node_per_trial_part[i] = nodes[i] / trials[i];
                        time_per_node_part[i] = time_part_sum[i] / nodes[i];
                    }
                    
                    if (childs[i]){
                        child_per_trial_part[i] = childs[i] / trials[i];
                        time_per_child_part[i] = time_part_sum[i] / childs[i];
                    }
                }
            }
            
            if (trial_sum){
                if (node_sum){
                    node_per_trial_all = node_sum / trial_sum;
                    time_per_node_all = time_all_sum / node_sum;
                }
                if (child_sum){
                    child_per_trial_all = child_sum / trial_sum;
                    time_per_child_all = time_all_sum / child_sum;
                }
            }
            
            
            std::cerr << " < Search Analysis >" << std::endl;
            std::cerr << " " << node_sum << " nodes ( " << node_per_trial_all << " nodes/trial ; " << time_per_node_all << " clock/node).  ";
            std::cerr << " " << child_sum << " childs ( " << child_per_trial_all << " childs/trial ; " << time_per_child_all << " clock/child).  " << std::endl;
            if (PART_NUM > 1){
                for (auto i = 0; i < PART_NUM; ++i){
                    std::cerr << "   Part " << i << " : " << nodes[i] << " nodes ( " << node_per_trial_part[i] << " nodes/trial ; " << time_per_node_part[i] << " clock/node).  ";
                    std::cerr << childs[i] << " childs ( " << child_per_trial_part[i] << " childs/trial ; " << time_per_child_part[i] << " clock/child)." << std::endl;
                }
            }
        }
#endif
    }
    
    StaticAnalyzer(){ init(); }
    
    StaticAnalyzer(const std::string& argName)
#ifdef USE_ANALYZER
    :name(argName)
#endif
    {
        init();
    }
    
    ~StaticAnalyzer(){
        report();
    }
};

struct HashBookAnalyzer{
    //置換表の解析用
    
#ifdef USE_ANALYZER
    const std::string name;
    const uint64 entry;
    const uint64 memory;
    
    // state
    std::atomic<uint64> filled;
    std::atomic<uint64> deleted;
    
    // read
    std::atomic<uint64> hit;
    std::atomic<uint64> unfounded;
    std::atomic<uint64> white;
    
    // regist
    std::atomic<uint64> registration;
    std::atomic<uint64> registrationFailure;
    
    void init()noexcept{
        filled = 0; deleted = 0;
        hit = 0; white = 0; unfounded = 0;
        registration = 0; registrationFailure = 0;
    }
    
    void addFilled()noexcept{ ++filled; }
    void addDeleted()noexcept{ ++deleted; }
    
    void addHit()noexcept{ ++hit; }
    void addWhite()noexcept{ ++white; }
    void addUnfounded()noexcept{ ++unfounded; }
    
    void addRegistration()noexcept{ ++registration; }
    void addRegistrationFailure()noexcept{ ++registrationFailure; }
    
    void report()const{
        using std::cerr; using std::endl;
        double realRate = filled / (double)entry;
        double idealRate = 1.0 - pow((1.0 - 1.0 / (double)entry), registration);
        cerr << "****** Analysis(HashBook) of " << name;
        cerr << " (" << entry << " entries, " << memory << " bytes) ******" << endl;
        cerr << "State  : filled = " << realRate << " ( ideally... " << idealRate << " )" << endl;
        cerr << "Regist : success = " << registration << "  failure = " << registrationFailure << endl;
        cerr << "Read   : hit = " << hit << "  unfounded = " << unfounded << "  white = " << white << endl;
    }
    
#else
    
    void init()const noexcept{}
    
    void addFilled()const noexcept{}
    void addDeleted()const noexcept{}
    
    void addHit()const noexcept{}
    void addWhite()const noexcept{}
    void addUnfounded()const noexcept{}
    
    void addRegistration()const noexcept{}
    void addRegistrationFailure()const noexcept{}
    
    void report()const noexcept{}
    
#endif
    
    //HashBookAnalyzer(){init();}
    
    HashBookAnalyzer(const std::string& argName, const int aent, const int amem)
#ifdef USE_ANALYZER
    :name(argName), entry(aent), memory(amem)
#endif
    {
        init();
    }
    
    ~HashBookAnalyzer(){
        report();
    }
};


#endif // UTIL_ANALYZER_HPP_