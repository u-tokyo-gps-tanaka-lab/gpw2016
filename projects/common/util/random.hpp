/*
 random.hpp
 Katsuki Ohto
 */

#ifndef UTIL_RANDOM_HPP_
#define UTIL_RANDOM_HPP_

#include <random>
#include <cmath>
#include <array>
#include <iostream>
#include <sstream>

#include "math.hpp"

// 確率分布いろいろ

// mean() ... 平均値
// med() ... 中央値
// mod() ... 最頻値
// var() ... 分散
// sd() ...　標準偏差
// cov(i, j) ... 共分散
// skew() ... 歪度
// kur() ... 尖度

// dens(x) ... 確率密度
// relative_dens(x) ... 相対確率密度(計算量削減)

// dist(x) ... 累積分布関数

// ent ... エントロピー

// rand(pdice) ... 確率密度に従う乱数発生



// 指数分布
struct ExponentialDistribution{
    double lambda;
    
    template<class dice_t>
    double rand(dice_t *const pdice)const{
        return (-1.0 / lambda * (1.0 - pdice->drand()));
    }
    
    double relative_dens(double x)const noexcept{
        return exp(-lambda * x);
    }
    double dens(double x)const noexcept{
        return lambda * relative_dens(x);
    }
    double dist(double x)const noexcept{
        return 1 - exp(-lambda * x);
    }
    double mean()const{
        return 1 / lambda;
    }
    double med()const{
        return log(2) / lambda;
    }
    double var()const{
        return 1 / (lambda * lambda);
    }
    double ent()const{
        return 1 - log(lambda);
    }
    
    static constexpr double mod()noexcept{ return 0; }
    static constexpr double skew()noexcept{ return 2; }
    static constexpr double kur()noexcept{ return 6; }
    
    ExponentialDistribution& set(double l)noexcept{
        lambda = l;
        return *this;
    }
    
    std::string toString()const{
        std::ostringstream oss;
        oss << "Exp(" << lambda << ")";
        return oss.str();
    }
    
    constexpr ExponentialDistribution():
    lambda(){}
    explicit constexpr ExponentialDistribution(double l):
    lambda(l){}
};

std::ostream& operator<<(std::ostream& out, const ExponentialDistribution& e){
    out << e.toString();
    return out;
}

// ガンマ分布
struct GammaDistribution{
    double a;
    
    template<class dice_t>
    double rand(dice_t *const pdice)const{
        double x, y, z;
        double u, v, w, b, c, e;
        if (a < 1){
            ExponentialDistribution ex(1);
            e = ex.rand(pdice);
            do {
                x = pow(pdice->drand(), 1 / a);
                y = pow(pdice->drand(), 1 / (1 - a));
            } while (x + y > 1);
            return (e * x / (x + y));
        }else{
            b = a - 1;
            c = 3 * a - 0.75;
            while(true){
                u = pdice->drand();
                v = pdice->drand();
                w = u * (1 - u);
                y = sqrt(c / w) * (u - 0.5);
                x = b + y;
                if(x >= 0){
                    z = 64 * w * w * w * v * v;
                    if(z <= 1 - (2 * y * y) / x){
                        return x;
                    }else{
                        if(log(z) < 2 * (b * log(x / b) - y)){
                            return x;
                        }
                    }
                }
            }
            return x;
        }
    }
    void add(double arg)noexcept{
        a += arg;
    }
    void subtr(double arg)noexcept{
        a -= arg;
    }
    void set(double arg)noexcept{
        a = arg;
    }
    
    constexpr GammaDistribution():
    a(){}
    explicit constexpr GammaDistribution(double aa):
    a(aa){}
};

// ベータ分布
struct BetaDistribution{
    double a, b;
    
    template<class dice_t>
    double rand(dice_t *const pdice)const{
        double r1 = GammaDistribution(a).rand(pdice);
        double r2 = GammaDistribution(b).rand(pdice);
        return r1 / (r1 + r2);
    }
    constexpr double size()const noexcept{
        return a + b;
    }
    constexpr double mean()const{
        return a / (a + b);
    }
    constexpr double rate()const{
        return a / b;
    }
    double var()const{
        double sum = a + b;
        return (a * b) / (sum * sum * (sum + 1.0));
    }
    double sd()const{
        double sum  = a + b;
        return std::sqrt(a * b / (sum + 1.0)) / sum;
    }
    double relative_dens(double x)const{
        return pow(x, a - 1) * pow(1 - x, b - 1);
    }
    double dens(double x)const{
        return relative_dens(x) / beta(a, b);
    }
    double log_relative_dens(double x)const{
        return (a - 1) * log(x) + (b - 1) * log(1 - x);
    }
    double log_dens(double x)const{
        DERR << log_beta(a, b) << endl;
        return log_relative_dens(x) - log_beta(a, b);
    }
    
    void add(const BetaDistribution& arg)noexcept{
        a += arg.a;
        b += arg.b;
    }
    void subtr(const BetaDistribution& arg)noexcept{
        a -= arg.a;
        b -= arg.b;
    }
    
    void operator+=(const BetaDistribution& rhs)noexcept{
        a += rhs.a;
        b += rhs.b;
    }
    void operator-=(const BetaDistribution& rhs)noexcept{
        a -= rhs.a;
        b -= rhs.b;
    }
    void operator*=(const double m)noexcept{
        a *= m;
        b *= m;
    }
    void operator/=(const double d){
        (*this) *= 1 / d;
    }
    
    void mul(const double m)noexcept{
        a *= m;
        b *= m;
    }
    
    BetaDistribution& rev()noexcept{
        std::swap(a, b);
        return *this;
    }
    
    BetaDistribution& set(const double aa, const double ab)noexcept{
        a = aa;
        b = ab;
        return *this;
    }
    BetaDistribution& set_by_mean(double m, double size)noexcept{
        set(m * size, (1 - m) * size);
        return *this;
    }
    BetaDistribution& set_by_rate(double r, double size)noexcept{
        set_by_mean(r / (1 + r), size);
        return *this;
    }
    
    BetaDistribution& resize(double h){
        // サイズをhにする
        double s = size();
        
        assert(s);
        
        double h_s = h / s;
        a *= h_s;
        b *= h_s;
        return *this;
    }
    
    bool exam()const noexcept{
        if(a < 0 || b < 0){ return false; }
        if(a == 0 && b == 0){ return false; }
        return true;
    }
    
    std::string toString()const{
        std::ostringstream oss;
        oss << "Be(" << a << ", " << b << ")";
        return oss.str();
    }
    
    constexpr BetaDistribution():
    a(), b(){}
    explicit constexpr BetaDistribution(const double aa, const double ab):
    a(aa), b(ab){}
};

std::ostream& operator<<(std::ostream& out, const BetaDistribution& b){
    out << b.toString();
    return out;
}

// ディリクレ分布
// 各確率変数が報酬を持つとし、その報酬を返すことも可能とする
template<std::size_t N>
struct DirichletDistribution{
    
    std::array<double, N> a;
    double sum;
    
    static constexpr std::size_t dim()noexcept{ return N; }
    
    template<class dice_t, typename reward_t>
    double rand(reward_t reward[], dice_t *const pdice)const{
        double sum = 0.0;
        reward_t sum_rew = 0.0;
        for(auto i = 0; i < N ; ++i){
            double r = GammaDistribution(a[i]).rand(pdice);
            sum += r;
            sum_rew += reward[i] * r;
        }
        sum_rew /= sum;
        return sum_rew;
    }
    
    template<class dice_t>
    std::array<double, N> rand(dice_t *const pdice)const{
        std::array<double, N> ans;
        double sum = 0;
        for(auto i = 0; i < N ; ++i){
            double r = GammaDistribution(a[i]).rand(pdice);
            ans[i] = r;
            sum += r;
        }
        double _sum = 1 / sum;
        for(auto i = 0; i < N ; ++i){
            ans[i] *= _sum;
        }
        return ans;
    }
    
    double size()const noexcept{
        return sum;
    }
    double mean(int i)const{
        return a[i] / sum;
    }
    std::array<double, N> mean()const{
        std::array<double, N> ans;
        for(auto i = 0; i < N ; ++i){
            ans[i] = mean(i);
        }
        return ans;
    }
    
    double var(int i)const{
        return a[i] * (sum - a[i]) / (sum * sum * (sum + 1.0));
    }
    std::array<double, N> var()const{
        std::array<double, N> ans;
        for(auto i = 0; i < N ; ++i){
            ans[i] = var(i);
        }
        return ans;
    }
    
    double cov(int i, int j)const{
        return - a[i] * a[j] / (sum * sum * (sum + 1.0));
    }
    std::array<std::array<double, N>, N> cov()const{
        std::array<std::array<double, N>, N> ans;
        for(auto i = 0; i < N ; ++i){
            for(auto j = 0; j < N; ++j){
                ans[i][j] = cov(i, j);
            }
        }
        return ans;
    }
    
    double relative_dens(const std::array<double, N>& x)const{
        double m = 1;
        for(auto i = 0; i < N ; ++i){
            m *= pow(x[i], a[i] - 1);
        }
        return m;
    }
    double dens(const std::array<double, N>& x)const{
        return relative_dens(x) / multivariate_beta(a, sum);
    }
    double log_relative_dens(const std::array<double, N>& x)const{
        double r = 0;
        for(auto i = 0; i < N ; ++i){
            r += (a[i] - 1) * log(x[i]);
        }
        return r;
    }
    double log_dens(const std::array<double, N>& x)const{
        return log_relative_dens(x) - log_multivariate_beta(a, sum);
    }
    
    DirichletDistribution& set(std::size_t i, double v){
        ASSERT(i < N, cerr << i << " in " << N << endl;);
        sum -= a[i];
        a[i] = v;
        sum += v;
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
        return *this;
    }
    DirichletDistribution& add(std::size_t i, double v){
        ASSERT(i < N, cerr << i << " in " << N << endl;);
        a[i] += v;
        sum += v;
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
        return *this;
    }
    DirichletDistribution& subtr(std::size_t i, double v){
        ASSERT(i < N, cerr << i << " in " << N << endl;);
        a[i] -= v;
        sum -= v;
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
        return *this;
    }
    
    void operator+=(const DirichletDistribution& rhs)noexcept{
        for(auto i = 0; i < N; ++i){
            a[i] += rhs.a[i];
            sum += rhs.a[i];
        }
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
    }
    void operator-=(const DirichletDistribution& rhs)noexcept{
        for(auto i = 0; i < N; ++i){
            a[i] -= rhs.a[i];
            sum -= rhs.a[i];
        }
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
    }
    void operator*=(const double m)noexcept{
        for(auto i = 0; i < N; ++i){
            a[i] *= m;
        }
        sum *= m;
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
    }
    void operator/=(const double d){
        (*this) *= 1 / d;
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
    }
    
    DirichletDistribution& set(const double aa[]){
        sum = 0;
        for(auto i = 0; i < N; ++i){
            double t = aa[i];
            a[i] = t;
            sum += t;
        }
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
        return *this;
    }
    DirichletDistribution& set(const std::array<double, N>& aa){
        sum = 0;
        for(auto i = 0; i < N; ++i){
            double t = aa[i];
            a[i] = t;
            sum += t;
        }
        ASSERT(exam(), std::cerr << toDebugString() << std::endl;);
        return *this;
    }
    
    std::string toString()const{
        std::ostringstream oss;
        oss << "Di(";
        for(auto i = 0; i < N - 1; ++i){
            oss << a[i] << ", ";
        }
        if(N > 0){
            oss << a[N - 1];
        }
        oss << ")";
        return oss.str();
    }
    std::string toDebugString()const{
        std::ostringstream oss;
        oss << toString() << " sum = " << sum;
        return oss.str();
    }
    
    DirichletDistribution(){
        std::array<double, N> tmp;
        tmp.fill(0);
        set(tmp);
    }
    DirichletDistribution(const double aa[]){
        set(aa);
    }
    DirichletDistribution(const std::array<double, N>& aa){
        set(aa);
    }
    
    bool exam()const noexcept{
        double tsum = 0;
        for(double v : a){
            tsum += v;
        }
        if(fabs(sum - tsum) > 0.0000001 * sum){ return false; }
        return true;
    }
};

template<std::size_t N>
std::ostream& operator<<(std::ostream& out, const DirichletDistribution<N>& d){
    out << d.toString();
    return out;
}

// 正規分布
template<typename f_t>
struct NormalDistribution{
    f_t mu, sigma;
    
    template<class dice_t>
    f_t rand(dice_t *const pdice)const{
        // Box-Muller
        f_t r1 = pdice->drand();
        f_t r2 = pdice->drand();
        f_t z1 = std::sqrt(-2.0 * std::log(r1)) * std::cos(2.0 * M_PI * r2);
        return z1 * sigma + mu;
    }
    template<class dice_t>
    void rand(f_t *const pa, f_t *const pb, dice_t *const pdice)const{
        // 2つ同時に発生させる
        f_t r1 = pdice->drand();
        f_t r2 = pdice->drand();
        
        //CERR<<r1<<","<<r2<<" : ";
        f_t z1 = std::sqrt(-2.0 * std::log(r1)) * std::cos(2.0 * M_PI * r2);
        f_t z2 = std::sqrt(-2.0 * std::log(r1)) * std::sin(2.0 * M_PI * r2);
        *pa = z1 * sigma + mu;
        *pb = z2 * sigma + mu;
        return;
    }
    
    f_t relative_dens(f_t x)const{
        return exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
    }
    f_t dens(f_t x)const{
        return relative_dens(x) / sigma * (1 / sqrt(2 * M_PI));
    }
    f_t dist(f_t x)const{
        return (1 + erf((x - mu) / sigma * (1 / sqrt(2)))) / 2;
    }
    
    f_t ent()const{
        return sigma * sqrt(2 * M_PI * exp(1));
    }
    
    std::string toString()const{
        std::ostringstream oss;
        oss << "N(" << mu << ", " << sigma << ")";
        return oss.str();
    }
    
    constexpr NormalDistribution():mu(),sigma(){}
    constexpr NormalDistribution(f_t argMu, f_t argSigma)
    :mu(argMu), sigma(argSigma){}
    
    NormalDistribution& set(f_t argMu, f_t argSigma)noexcept{
        mu = argMu;
        sigma = argSigma;
        return *this;
    }
};

template<class float_t>
std::ostream& operator<<(std::ostream& out, const NormalDistribution<float_t>& n){
    out << n.toString();
    return out;
}

#endif // UTIL_RANDOM_HPP_