/*
 fastSimulator.h
 Katsuki Ohto
 */

// Digital Curling
// Fast Simulation using array ( floating point number calculation )

#include "../dc.hpp"
#include "b2dSimulator.hpp"
#include "simuFunc.hpp"

#ifndef DCURLING_SIMULATION_FASTSIMULATOR_H_
#define DCURLING_SIMULATION_FASTSIMULATOR_H_

namespace DigitalCurling{
    namespace FastSimulator{

        /**************************物理モデルの性質にあまり依存しない定数、関数*************************/

        constexpr fpn_t FR_COLLISION = 0.00004; // 衝突閾値距離

        // 全区間のV_ALL変換
        struct FSimuData{
            // 原点(0,0)で,(0,1)方向に静止する石の、
            // それぞれの地点での状態
            fpn_t x, y; // 位置
            fpn_t th; // 運動の向き
            fpn_t r;

            void set(fpn_t ax, fpn_t ay, fpn_t ath, fpn_t ar)noexcept{
                x = ax; y = ay; th = ath; r = ar;
            }
        };

        // V -> R Conversion
        constexpr int FV_FR_TABLE_RESOLUTION = 10;
        constexpr fpn_t FV_SIMU_MAX = FVY_LEGAL_MAX + 2.0;
        constexpr fpn_t FV_TABLE_MAX = FV_SIMU_MAX + 0.1;

        constexpr pfpn_t PFV_FR_TABLE_STEP = pfpn_t(FV_TABLE_MAX) / pfpn_t(1 << FV_FR_TABLE_RESOLUTION);
        constexpr pfpn_t PFV_FR_TABLE_INVSTEP = pfpn_t(1 << FV_FR_TABLE_RESOLUTION) / pfpn_t(FV_TABLE_MAX);

        constexpr fpn_t FV_FR_TABLE_STEP = (fpn_t)PFV_FR_TABLE_STEP;
        constexpr fpn_t FV_FR_TABLE_INVSTEP = (fpn_t)PFV_FR_TABLE_INVSTEP;

        fpn_t FVtoFR(fpn_t);

        // R -> V Conversion (Long Range)
        constexpr int FR_FV_TABLE_RESOLUTION = 12;
        constexpr fpn_t FR_TABLE_MAX = 49.0;

        constexpr pfpn_t PFR_FV_TABLE_STEP = pfpn_t(FR_TABLE_MAX) / pfpn_t(1 << FR_FV_TABLE_RESOLUTION);
        constexpr pfpn_t PFR_FV_TABLE_INVSTEP = pfpn_t(1 << FR_FV_TABLE_RESOLUTION) / pfpn_t(FR_TABLE_MAX);

        constexpr fpn_t FR_FV_TABLE_STEP = (fpn_t)PFR_FV_TABLE_STEP;
        constexpr fpn_t FR_FV_TABLE_INVSTEP = (fpn_t)PFR_FV_TABLE_INVSTEP;

        fpn_t FRtoFV(fpn_t);

        // R -> V Conversion (Short Range)
        constexpr int FR_FV_SR_TABLE_RESOLUTION = 12;
        constexpr fpn_t FR_SR_TABLE_MAX = FR_TABLE_MAX / 256;

        constexpr pfpn_t PFR_FV_SR_TABLE_STEP = pfpn_t(FR_SR_TABLE_MAX) / pfpn_t(1 << FR_FV_SR_TABLE_RESOLUTION);
        constexpr pfpn_t PFR_FV_SR_TABLE_INVSTEP = pfpn_t(1 << FR_FV_SR_TABLE_RESOLUTION) / pfpn_t(FR_SR_TABLE_MAX);

        constexpr fpn_t FR_FV_SR_TABLE_STEP = (fpn_t)PFR_FV_SR_TABLE_STEP;
        constexpr fpn_t FR_FV_SR_TABLE_INVSTEP = (fpn_t)PFR_FV_SR_TABLE_INVSTEP;

        fpn_t FRtoFVinSR(fpn_t);

        // V -> All Conversion
        constexpr int FV_TABLE_RESOLUTION = 12;

        constexpr pfpn_t PFV_TABLE_STEP = pfpn_t(FV_TABLE_MAX) / pfpn_t(1 << FV_TABLE_RESOLUTION);
        constexpr pfpn_t PFV_TABLE_INVSTEP = pfpn_t(1 << FV_TABLE_RESOLUTION) / pfpn_t(FV_TABLE_MAX);

        constexpr fpn_t FV_TABLE_STEP = (fpn_t)PFV_TABLE_STEP;
        constexpr fpn_t FV_TABLE_INVSTEP = (fpn_t)PFV_TABLE_INVSTEP;

        uint32 FVtoID(fpn_t fv)noexcept{ assert(fv >= 0); return (uint32)(fv * FV_TABLE_INVSTEP); }
        fpn_t FVtoFRAC(fpn_t fv)noexcept{ assert(fv >= 0); return fv - (int(FVtoID(fv))) * FV_TABLE_STEP; }

        template<typename float_t = fpn_t>
        constexpr fpn_t IDtoFV(uint32 idx)noexcept{
            return ((float_t)idx) * FV_TABLE_STEP;
        }

        fpn_t FVtoFX(fpn_t); fpn_t FVtoFY(fpn_t); fpn_t FVtoFTh(fpn_t);
        fpn_t FVtoFRbyAllTable(fpn_t);

        // R_ALL変換
        uint32 FRtoID(fpn_t fr){ return FVtoID(FRtoFV(fr)); }
        fpn_t IDtoFR(uint32 idx){ return FVtoFR(IDtoFV(idx)); }
        
        // R_CZ変換
        constexpr int FR_CZ_TABLE_RESOLUTION = 12;
        constexpr fpn_t FR_MAX_GTOS = 24.0; // ゴールからの距離の最大値
        
        constexpr fpn_t PFR_CZ_TABLE_STEP = pfpn_t(FR_MAX_GTOS) / pfpn_t(1 << FR_CZ_TABLE_RESOLUTION);
        constexpr fpn_t FR_CZ_TABLE_STEP = (fpn_t)PFR_CZ_TABLE_STEP;

        /**************************物理モデルの性質に強く依存する定数、関数*************************/

        // V <-> T Conversion (比例関係)
        constexpr fpn_t FTtoFV(fpn_t ft)noexcept{ return F_FRIC_RINK_STONE * ft; }
        constexpr fpn_t FVtoFT(fpn_t fv)noexcept{ return fv / F_FRIC_RINK_STONE; }
        constexpr fpn_t FDTtoFDV(fpn_t ft)noexcept{ return F_FRIC_RINK_STONE * ft; }
        constexpr fpn_t FDVtoFDT(fpn_t fv)noexcept{ return fv / F_FRIC_RINK_STONE; }

        // R <-> T Conversion
        fpn_t FTtoFR(fpn_t ft){ return FVtoFR(FTtoFV(ft)); }
        fpn_t FRtoFT(fpn_t fr){ return FVtoFT(FRtoFV(fr)); }

        // T -> ALL Conversion
        uint32 FTtoID(fpn_t ft)noexcept{ return FVtoID(FTtoFV(ft)); }
        fpn_t IDtoFT(uint32 idx)noexcept{ return FVtoFT(IDtoFV(idx)); }

        // V,THETA -> R,THETA Convertion
        
        // 第1回UEC杯~の物理モデルにおいて、
        // 石の速度ベクトルと、現在位置->停止位置の向きベクトルの間のなす角度が一定らしい
        // TODO: 運動方程式からの数理的な証明
        // ここで定義されている関数はその性質を利用した演算なので
        // 物理モデルが変わったら即座に使用を取りやめる必要があるので注意
        
        // 対数螺旋に関係する定数
        constexpr fpn_t FTHETA_CURL_MAX = 0.03333565; // 停止位置への回転最大角度
        constexpr fpn_t FTHETA_CURL = 0.03333564; // 停止位置への回転角度
        constexpr fpn_t FA_CURVE = 3.45628911574e-19;
        constexpr fpn_t FB_CURVE = 29.986811440344486;
        
        // 停止位置までの定数
        constexpr fpn_t FDX_FV_FR = 0.001386965639075;
        constexpr fpn_t FDY_FV_FR = 0.041588442394742;
        constexpr fpn_t FDR_FV_FR = XYtoR(FDX_FV_FR, FDY_FV_FR);
        
        void FVThtoFRTh(fpn_t *const, fpn_t *const, fpn_t, fpn_t, int);
        void FVThtoFXYRTh(fpn_t *const, fpn_t *const, fpn_t *const, fpn_t *const,
            fpn_t, fpn_t, Spin);
        void FXYVThtoFXYRTh(fpn_t *const, fpn_t *const, fpn_t *const, fpn_t *const,
            fpn_t, fpn_t, fpn_t, fpn_t, Spin);

        /**************************物理演算用の石クラス*************************/
        
        struct FMoment{
            fpn_t x, y;
            fpn_t v_, th_;
            fpn_t w;
            fpn_t t;

            void set(fpn_t ax, fpn_t ay, fpn_t av, fpn_t ath, fpn_t aw, fpn_t at)noexcept{
                x = ax; y = ay; v_ = av; th_ = ath; w = aw; t = at;
            }

            constexpr fpn_t getX()const noexcept{ return x; }
            constexpr fpn_t getY()const noexcept{ return y; }
            constexpr fpn_t v()const noexcept{ return v_; }
            constexpr fpn_t v2()const noexcept{ return v_ * v_; }
            constexpr fpn_t th()const noexcept{ return th_; }
            constexpr fpn_t vx()const noexcept{ return v_ * sin(th_); }
            constexpr fpn_t vy()const noexcept{ return v_ * cos(th_); }
            constexpr fpn_t getW()const noexcept{ return w; }
            constexpr Spin spin()const noexcept{ return FWtoSpin(w); }
            constexpr bool isLeftSpin()const noexcept{ return isLeftFW(w); }
            
            void setV(fpn_t av)noexcept{ v_ = av; }
            void setTh(fpn_t ath)noexcept{ th_ = ath; }

            std::string toString()const{
                std::ostringstream oss;
                oss << "x = " << x << " y = " << y
                    << " vx = " << vx() << " vy = " << vy() << " w = " << w << " t = " << t;
                return oss.str();
            }

            void setGoalInfo()const noexcept{}
            void stepGoalInfoByNextR(fpn_t nr)const noexcept{}
            void stepGoalInfo()const noexcept{}
            
            constexpr bool exam()const noexcept{ return true; }
        };

        struct FSMovingStone : public FMoment{
            BitSet32 asleep; // 自分と衝突可能性がある静止石のデータ
            BitSet32 awake; // 自分と衝突可能性がある動石のデータ
            fpn_t gx, gy; // position stopped without collision
            fpn_t gr; // 停止位置までの距離
            fpn_t gth; // theta to goal
            fpn_t gt; // goal time
            fpn_t dth; // theta distance between table

            StoneState state;
            int contactIndex; // 衝突した際に割り振られるインデックス(衝突木生成用)

            void stepGoalInfoByNextR(fpn_t nr)noexcept{
                gr = nr;
            }
            void stepGoalInfo()noexcept{
                gr = XYtoR(gx - x, gy - y);
            }

            void setGoalInfo()noexcept{
                gt = t + FDVtoFDT(v_);
                FXYVThtoFXYRTh(&gx, &gy, &gr, &gth, x, y, v_, th_, FWtoSpin(w));
                dth = (isRight(FWtoSpin(w))) ? (th_ - FVtoFTh(v_)) : (th_ + FVtoFTh(v_));
                //cerr<<FVtoFTh(v)<<endl;
            }

            void set(fpn_t ax, fpn_t ay, fpn_t av, fpn_t ath, fpn_t aw, fpn_t at)noexcept{
                FMoment::set(ax, ay, av, ath, aw, at);
                setGoalInfo();
                ASSERT(exam(),);
            }
            
            void init(StoneState ast)noexcept{
                contactIndex = -1;
                state = ast;
            }

            std::string toString()const{
                std::ostringstream oss;
                oss << FMoment::toString() << endl;
                oss << "gx = " << gx << " gy = " << gy
                    << " gr = " << gr << " gth = " << gth << " gt = " << gt << " dth = " << dth;
                return oss.str();
            }

            bool exam()const{
                bool ret = true;
                fpn_t tgx, tgy, tgr, tgth, tgt;
                FXYVThtoFXYRTh(&tgx, &tgy, &tgr, &tgth, x, y, v_, th_, FWtoSpin(w));
                tgt = t + FDVtoFDT(v_);
                if (fabs(tgx - gx) > 0.01){ cerr << "FSMovingStone::exam() : illegal gx = " << gx << " expected = " << tgx << endl; ret = false; }
                if (fabs(tgy - gy) > 0.01){ cerr << "FSMovingStone::exam() : illegal gy = " << gy << " expected = " << tgy << endl; ret = false; }
                if (fabs(tgr - gr) > 0.01){ cerr << "FSMovingStone::exam() : illegal gr = " << gr << " expected = " << tgr << endl; ret = false; }
                //if( fabs(tgth-gth)>0.01 ){cerr<<"FSMovingStone::exam() : illegal gth"<<gth<<" expected = "<<tgth<<endl;endl;ret=false;}
                if (fabs(tgt - gt) > 0.01){ cerr << "FSMovingStone::exam() : illegal gt = " << gt << " expected = " << tgt << endl; ret = false; }
                if (!ret){ cerr << toString() << endl; }
                return ret;
            }
        };

        /**************************物理演算用の局面クラス*************************/
        
        struct FastSimulatorBoard{
            
            // 盤面に静止石が2個以上あったような局面において
            // 初めての衝突が起こったらこちらのクラスに石の位置データをセットする
            // TODO: 思考時の局面情報クラスでそのままやるように書き換えた方がいいと思う
            
            using stone_t = FSMovingStone;
            
            stone_t stone_[N_STONES]; // 石の状態
            //std::pair<int,int> pList[16*(16-1)/2]; // 衝突可能性が否定出来ないawake_awakeの組み合わせ
            //int awakeList[16]; // awakeな石を時間の進んでいない順に並べたもの
            
            BitSet32 awake;
            int NAwake;
            BitArray64<4, 16> awakeArray; // awakeな石を衝突可能性がある最短時刻順に並べたもの
            
            BitSet32 asleep;
            int NAsleep;
            int turn; // 投げられた石番に等しいことを前提とする

            //int NStones; // 全石数

            int getTurnNum()const noexcept{ return turn; }
            
            const FSMovingStone& stone(int idx)const{ return stone_[idx]; }
            FSMovingStone& stone(int idx){ return stone_[idx]; }

            template<class board_t>
            void setAfterFirstColl(const board_t& bd, uint32, uint32); // 最初の衝突後にこのクラスに情報をセット
            
            void growContactTreeFirst(ContactTree *const pct, int iaw, int ias){
                // 最初の衝突時の ContactTree の生成
                pct->setContact(0, stone(ias).state.getIndex());
                stone(iaw).contactIndex = 1;
                stone(ias).contactIndex = 2;
            }
            void growContactTree(ContactTree *const pct, int i0, int i1){
                // 衝突が起きた場合に ContactTree に情報を追加する
                int ci0 = stone(i0).contactIndex, ci1 = stone(i1).contactIndex;
                if(!pct->isComplex()){
                    // 汎用物理エンジンに投げた場合contactIndexは正しくないが、それ以外の場合のチェック
                    ASSERT(ci0 != ci1, cerr << "ci0 = ci1 = " << ci0 << endl;);
                }
                if(ci0 > ci1){
                    pct->setContact(ci0, stone(i1).state.getIndex());
                    int nci = ContactTree::ItoNextI(ci0);
                    stone(i0).contactIndex = nci; stone(i1).contactIndex = nci + 1;
                }else{
                    pct->setContact(ci1, stone(i0).state.getIndex());
                    int nci = ContactTree::ItoNextI(ci1);
                    stone(i1).contactIndex = nci; stone(i0).contactIndex = nci + 1;
                }
            }
            
            bool exam()const{
                /*if( NStones != NAwake + asleep.count() ){
                cerr<<" FastSimulatorBoard::exam : invalid num of stones "<<NStones<<" != "<<NAwake<<" + "<<asleep.count()<<endl;
                return false;
                }*/

                // awake, asleep list
                if (!isExclusive((uint32)awake, (uint32)asleep)){
                    cerr << "FastSimulatorBoard::exam() : awakes and asleeps inexclusive." << endl;
                    return false;
                }
                for (BitSet32 tmp = awake; tmp.any(); tmp.pop_lsb()){
                    int iaw = tmp.bsf();
                    if (!awake.holds(stone(iaw).awake)){
                        cerr << "FastSimulatorBoard::exam() : all awakes doesn't hold aw_stone(" << iaw << ")'s awakes." << endl;
                        return false;
                    }
                    if (!asleep.holds(stone(iaw).asleep)){
                        cerr << "FastSimulatorBoard::exam() : all asleeps doesn't hold aw_stone(" << iaw << ")'s asleeps." << endl;
                        return false;
                    }
                    if (stone(iaw).awake.test(iaw)){
                        cerr << "FastSimulatorBoard::exam() : aw_stone(" << iaw << ")'s awakes holds myself." << endl;
                        return false;
                    }
                    if (stone(iaw).asleep.test(iaw)){
                        cerr << "FastSimulatorBoard::exam() : aw_stone(" << iaw << ")'s asleeps holds myself." << endl;
                        return false;
                    }
                    if (!awake.test(iaw)){
                        cerr << "FastSimulatorBoard::exam() : awakes doesn't hold aw_stone(" << iaw << ")" << endl;
                        return false;
                    }
                    if (asleep.test(iaw)){
                        cerr << "FastSimulatorBoard::exam() : asleeps holds aw_stone(" << iaw << ")" << endl;
                        return false;
                    }
                }
                // the number of stones
                int nst[2] = {0};
                for (BitSet32 tmp = awake; tmp.any(); tmp.pop_lsb()){
                    int iaw = tmp.bsf();
                    nst[stone(iaw).state.getColor()]++;
                }
                for (BitSet32 tmp = asleep; tmp.any(); tmp.pop_lsb()){
                    int ias = tmp.bsf();
                    nst[stone(ias).state.getColor()]++;
                }
                for (int c = 0; c < 2; ++c){
                    if (nst[c] > N_COLOR_STONES){
                        cerr << "FastSimulatorBoard : illegal number of stones( "
                        << toColorString(c) << nst[c] << " )." << endl;
                        return false;
                    }
                }
                return true;
            }
        };
        
        template<typename callback_t>
        void iterateAsleepStoneWithColor(const FastSimulatorBoard& bd, const callback_t& callback){
            iterate(bd.asleep,[&bd,callback](int idx)->void{
                callback(bd.stone(idx).state.getColor(), bd.stone(idx));
            });
        }
        template<typename callback_t>
        void iterateAsleepStoneWithIndex(const FastSimulatorBoard& bd, const callback_t& callback){
            iterate(bd.asleep,[&bd,callback](int idx)->void{
                callback(idx, bd.stone(idx));
            });
        }
        
        /**************************関数リスト*************************/

        // initializing
        void initV_R();
        void initR_V();
        void initR_VinSR();
        void initALL();
        void initR_CZ();
        void init();
        
        struct FastSimulatorInitializer{
            FastSimulatorInitializer(){
                init();
            }
        };
        FastSimulatorInitializer fsInitializer;

        // collision
        template<class mst_t, class st_t>
        void collisionAw_As(mst_t *const aw, st_t *const as);
        template<class mst0_t, class mst1_t>
        void collision2Aw(mst0_t *const aw0, mst1_t *const aw1);
        
        // collision possibility
        template<class pos0_t,class pos1_t>
        inline bool hasCollisionChanceAwAs(const pos0_t& aw, const pos1_t& as);

        // move <-> stop position
        template<class pos_t, class mv_t>
        void FXYtoFMV(const pos_t& pos, mv_t *const mv);

        template<class pos_t, class move_t>
        void FMVtoFXY(const move_t& mv, pos_t *const dst);
        
        template<class pos_t>
        void FV_FRtoFXY(fpn_t v, fpn_t r, pos_t *const pos, fpn_t *const theta, Spin spin);


        // fast step
        template<class mst_t>
        void step(mst_t *const dst, fpn_t fr, fpn_t fnr);

        template<class mst_t>
        void stepToStop(mst_t *const dst, fpn_t fr);

        template<class float_t>
        void stepToStop(float_t *const px, float_t *const py, float_t *const pt,
            const uint32 id, float_t theta, float_t w);

        template<int RINK_ONLY = 1, class board_t, class board2_t>
        int simulateF_b2d_FSB_inContact(board_t *const pbd, const fpn_t step);

        template<int RINK_ONLY = 1, class board_t, class move_t>
        ContactTree simulateSoloF(board_t *const pbd, const move_t& mv);

        template<int RINK_ONLY = 1, class board_t, class fsbd_t>
        ContactTree loopAfterFirstCollision(board_t *const pabd, fsbd_t *const pbd);

        template<int RINK_ONLY = 1, class board_t, class mst_t>
        ContactTree loop1_1(board_t *const, uint32 awNum, mst_t&);
        
        template<int RINK_ONLY = 1, class board_t, class mst_t>
        ContactTree loop1_N(board_t *const, uint32 awNum, mst_t&);

        template<int RINK_ONLY = 1, class board_t, class move_t>
        ContactTree simulateF(board_t *const pbd, uint32 stNum, const move_t& mv);

        template<class pos_t, class move_t>
        void rotateToPassPointF(move_t *const pmv, const pos_t& pos, fpn_t v){
            // 指定初速で、指定された位置を通る速度ベクトルを生成
            // ヒット系の着手の生成で呼ばれる

            const fPosXY<> tar(pos.x - FX_THROW, pos.y - FY_THROW);
            const fpn_t r = XYtoR(tar.getX(), tar.getY());

            fPosXY<> tp;
            fpn_t theta;
            
            FV_FRtoFXY(v, r, &tp, &theta, pmv->getSpin());

            fpn_t dtheta = XYtoT(tar.y, tar.x) - XYtoT(tp.y, tp.x) + theta;
            rotate(v, 0, dtheta, &pmv->y, &pmv->x);

            //CERR << tar << *mv << endl;
        }

        template<class pos_t, class move_t>
        void rotateToPassPointGoF(move_t *const pmv, const pos_t& pos, fpn_t r){
            // 指定された位置を通り、さらに指定距離分(投擲点から見て)進む速度ベクトルを生成
            // ヒット系の着手の生成で呼ばれる
            
            const fPosXY<> tar(pos.x - FX_THROW, pos.y - FY_THROW);
            const fpn_t pr = XYtoR(tar.getX(), tar.getY());
            const fpn_t v = FRtoFV(pr + r);

            fPosXY<> tp;
            fpn_t theta;
            
            FV_FRtoFXY(v, pr, &tp, &theta, pmv->getSpin());

            fpn_t dtheta = XYtoT(tar.y, tar.x) - XYtoT(tp.y, tp.x) + theta;

            rotate(v, 0, dtheta, &pmv->y, &pmv->x);
        }

        template<class pos_t, class move_t>
        void rotateToPassPointCircumferenceF(move_t *const pmv, const pos_t& pos, fpn_t av, fpn_t atheta, fpn_t ar = 2 * FR_STONE_RAD){
            // 指定初速で、指定された位置の石からの相対位置 ar, atheta の位置を通る速度ベクトルを生成
            // 基本は ar = 2 × 石半径、つまり石のどの点をたたくかという使い方
            rotateToPassPointF(pmv, fPosXY<>(pos.x + ar*sin(atheta), pos.y + ar*cos(atheta)), av);
        }

        // 軌道距離計算
        // awの軌道とasの位置の距離を近似計算
        //template<class pos_t>
        //fpn_t calcOrbitDist2AW_AS(const pos_t& org, const pos_t& stop, const int spin, const pos_t& as);
        //template<class pos_t>
        // fpn_t calcOrbitDistAW_AS(fpn_t v, fpn_t vtheta, const pos_t& as);
        
    }
    
    template<int RINK_ONLY, class board_t, class move_t>
    ContactTree makeMoveNoRand(board_t *const, uint32, const move_t&);
    
    template<int RINK_ONLY, class board_t, class move_t, class dice_t>
    ContactTree makeMove(board_t *const, uint32, const move_t&, dice_t *const);
}

#endif // DCURLING_SIMULATION_FASTSIMULATOR_H_