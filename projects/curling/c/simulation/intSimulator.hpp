//Digital Curling
//Fast Simulation using array ( integar calculation )

#include "../dc.hpp"
#include "simufunc.hpp"
#include "fastsimulator.h"
#include "intsimulator.h"

#ifndef DCURLING_SIMULATION_INTSIMULATOR_HPP_
#define DCURLING_SIMULATION_INTSIMULATOR_HPP_

namespace DigitalCurling{
    namespace IntSimulator{
        
        std::array<sim2int_t, (1 << IV_IR_TABLE_RESOLUTION) + 4> IV_IRTable;
        std::array<simint_t, (1 << IR_IV_TABLE_RESOLUTION) + 4> IR_IVTable;
        std::array<ISimuData, (1 << IR_TABLE_RESOLUTION) + 4> ISimuTable;
        
        sim2int_t ID_FRACtoX(uint32 id, uint32 frac){ return ISimuTable[id].x*((1U << IF_RESOLUTION) - frac) + ISimuTable[id + 1].x*frac; }
        sim2int_t ID_FRACtoY(uint32 id, uint32 frac){ return ISimuTable[id].y*((1U << IF_RESOLUTION) - frac) + ISimuTable[id + 1].y*frac; }
        int ID_FRACtoCOS(uint32 id, uint32 frac){ return ISimuTable[id].cos*((1U << IF_RESOLUTION) - frac) + ISimuTable[id + 1].cos*frac; }
        int ID_FRACtoSIN(uint32 id, uint32 frac){ return ISimuTable[id].sin*((1U << IF_RESOLUTION) - frac) + ISimuTable[id + 1].sin*frac; }
        
        //V_R変換
        sim2int_t IVtoIR(simint_t iv){
            uint32 vid = uint32(iv >> IV_RESOLUTION_OF_IV_IR_TABLE);
            sim2int_t hasu= sim2int_t(iv % IV_IR_TABLE_STEP);
            
            //CERR<<vid<<" - "<<hasu;
            ASSERT(0 <= vid && vid<(1 << IV_IR_TABLE_RESOLUTION) + 4, cerr << (I32VtoFV(iv) / (double)(FV_MAX+2.0)) << endl;);
            ASSERT(0 <= hasu && hasu < IV_IR_TABLE_STEP, cerr << hasu << endl;);
            
            return (IV_IRTable[vid] * ((sim2int_t)IV_IR_TABLE_STEP - hasu) + IV_IRTable[vid + 1] * hasu)>>IV_RESOLUTION_OF_IV_IR_TABLE;//16ビット余り
        }
        
        //R_V変換
        simint_t IRtoIV(sim2int_t ir){
            
            uint32 rid = uint32(ir >> IR_RESOLUTION_OF_IR_IV_TABLE);
            uint32 hasu = uint32((ir % IR_IV_TABLE_STEP)>>24);
            
            //CERR<<rid<<" - "<<hasu;
            ASSERT(0 <= rid && rid<(1 << IR_IV_TABLE_RESOLUTION) + 4, cerr << (I64RtoFR(ir) / (double)49.0) << endl;);
            ASSERT(0 <= hasu && hasu < I32R_IV_TABLE_STEP, cerr << hasu << endl;);
            
            return (IR_IVTable[rid] * (I32R_IV_TABLE_STEP - hasu) + IR_IVTable[rid + 1] * hasu)>>(IR_RESOLUTION_OF_IR_IV_TABLE-24);//7ビット余り
        }
        
        //T_V変換
        constexpr simint_t ITtoIV(simint_t it){ return I_FRIC_STONE * it; }
        //V_T変換
        constexpr simint_t IVtoIT(simint_t iv){ return iv / I_FRIC_STONE; }
        
        //T_R変換
        simint_t ITtoIR(simint_t it){ return IVtoIR(ITtoIV(it)); }
        //R_T変換
        simint_t IRtoIT(simint_t ir){ return IVtoIT(IRtoIV(ir)); }
        
        //T_ALL変換
        constexpr uint32 ITtoID(simint_t it){ return IVtoID(ITtoIV(it)); }
        constexpr simint_t IDtoIT(uint32 id){ return IVtoIT(IDtoIV(id)); }
        const ISimuData& ITtoFALL(simint_t it){ return V_AllIntTable[ITtoID(it)]; }
        
        //R_ALL変換
        uint32 IRtoID(simint_t ir){ return IVtoID(IRtoIV(ir)); }
        simint_t IDtoIR(uint32 id){ return IVtoIR(IDtoIV(id)); }
        
        
        //Base型
        //ALLテーブル参照のための基本型
        //constexpr int IBtoIV(int ib){ return ib; }
        //constexpr int IBtoIT(int ib){ return IVtoIT(IBtoIV(ib)); }
        //int IBtoIR(int ib){ return IVtoIR(IBtoIV(ib)); }
        
        struct IMoment{
            sim2int_t x, y;//current position
            sim2uint_t v;//current velocity
            simint_t cos,sin;//current direction of movement ( against table-data )
            sim2int_t w;//current angular velocity
            sim2uint_t t;//current time
            
            void set(sim2int_t ax, sim2int_t ay,
                     simint_t acos, simint_t asin,
                     sim2int_t av, sim2int_t aw
                     sim2int_t at,
                     ){
                x = ax; y = ay;
                cos = acos; sin = asin;
                v = av; w = aw;
                t = at;
            }
        };
        
        struct IMovingStone : public IMoment{
            sim2int_t gx,gy;//position stopped with no collision
            sim2uint_t gr;//distance from current position with no collision
            int gcos, gsin;
            sim2uint_t gt;//time stopped with no collision
            
            BitSet32 asleep;//asleep stones which have possibility to collide me
            BitSet32 awake;//awake stones which have possibility to collide me
            StoneState state;//information of stone
        };
        
        struct ISBoard{
            ISMovingStone stone[16];//石の状態
            //std::pair<int,int> pList[16*(16-1)/2];//衝突可能性が否定出来ないawake_awakeの組み合わせ
            //int NPList;
            
            //int awakeList[16];//awakeな石を時間の進んでいない順に並べたもの
            
            BitSet32 awake;
            int NAwake;
            BitSet32 asleep;
            int NAsleep;
            
            int col;
            
            //int NStones;//全石数
            
            int getTurnColor()const{ return col; }
            
            template<class board_t, class mst_t>
            void setAiterFirstColl(const board_t& bd, const int as_coll, const mst_t& aw0, const mst_t& aw1){
                //最初の衝突後にデータをセット
                stone[0].set(aw0.x, aw0.y, aw0.theta, aw0.t, aw0.v, aw0.w);
                stone[1].set(aw1.x, aw1.y, aw1.theta, aw1.t, aw1.v, aw1.w);
                
                awake.reset();
                awake.set(0); awake.set(1);
                NAwake = 2;
                //NPList=0;
                
                col = bd.getTurnColor();
                
                asleep = 0U;
                int cnt = NAwake;
                for (int c = 0; c<2; ++c){
                    for (int i = 0; i<bd.NActive[c]; ++i){
                        if (as_coll != c * 8 + i){
                            asleep.set(cnt);
                            stone[cnt].x = bd.stone[c][i].x;
                            stone[cnt].y = bd.stone[c][i].y;
                            stone[cnt].state = bd.state[c][i];
                            
                            DERR << bd.stone[c][i] << endl;
                            ++cnt;
                        }
                    }
                }
                
                NAsleep = cnt - NAwake;
                
                //NStones=cnt;
                
                //停止位置までの距離を入れる
                //stone[0].r_goal=IVtoIR(stone[0].v);
                //stone[1].r_goal=IVtoIR(stone[1].v);
                //動いている石のasleepListを作る
                stone[0].asleep = asleep;
                stone[1].asleep = asleep;
                
                stone[0].awake.reset();
                stone[1].awake.reset();
                
                //石の付加情報
                stone[0].state.init();
                stone[0].state.setColor(col);
                
                stone[1].state = bd.state[as_coll / 8][as_coll % 8];
            }
            bool exam()const{
                
                //awake,asleep list
                if (!isExclusive((uint32)awake, (uint32)asleep)){
                    cerr << "ISBoard::exam() : awakes and asleeps inexclusive." << endl;
                    return false;
                }
                for (BitSet32 tmp = awake; tmp.any(); tmp.pop_lsb()){
                    int aw = tmp.bsf();
                    if (!awake.holds(stone[aw].awake)){
                        cerr << "ISBoard::exam() : all awakes doesn't hold aw_stone(" << aw << ")'s awakes." << endl;
                        return false;
                    }
                    if (!asleep.holds(stone[aw].asleep)){
                        cerr << "ISBoard::exam() : all asleeps doesn't hold aw_stone(" << aw << ")'s asleeps." << endl;
                        return false;
                    }
                    if (stone[aw].awake.test(aw)){
                        cerr << "ISBoard::exam() : aw_stone(" << aw << ")'s awakes holds myself." << endl;
                        return false;
                    }
                    if (stone[aw].asleep.test(aw)){
                        cerr << "ISBoard::exam() : aw_stone(" << aw << ")'s asleeps holds myself." << endl;
                        return false;
                    }
                    if (!awake.test(aw)){
                        cerr << "ISBoard::exam() : awakes doesn't hold aw_stone(" << aw << ")" << endl;
                        return false;
                    }
                    if (asleep.test(aw)){
                        cerr << "ISBoard::exam() : asleeps holds aw_stone(" << aw << ")" << endl;
                        return false;
                    }
                }
                //the number of stones
                int nst[2] = { 0 };
                for (BitSet32 tmp = awake; tmp.any(); tmp.pop_lsb()){
                    int aw = tmp.bsf();
                    nst[stone[aw].state.getColor()]++;
                }
                for (BitSet32 tmp = asleep; tmp.any(); tmp.pop_lsb()){
                    int as = tmp.bsf();
                    nst[stone[as].state.getColor()]++;
                }
                for (int c = 0; c<2; ++c){
                    if (nst[c]>8){
                        cerr << "ISBoard : illegal number of stones( " << OutColor(c) << nst[c] << " )." << endl;
                        return false;
                    }
                }
                return true;
            }
        };
        //初期化
        
        //初期化スレッド
        void initV_R(){//V_R変換初期化
            constexpr pfpn_t step = static_cast<pfpn_t>(1) / static_cast<pfpn_t>(1000000);
            
            fMStone<pfpn_t> mst(0, 0, 0, 0, 0), lmst, tmst;
            FV_FRTable[0] = 0;
            pfpn_t tm = 0;//逆算時刻
            for (int i = 1; i<(1 << IV_IR_TABLE_RESOLUTION) + 4; ++i){
                //V_Rテーブルに現時点での値を入れていく
                pfpn_t line = I64VtoFV(IV_IR_TABLE_STEP * i);
                pfpn_t fv;
                while (1){
                    fv = XYtoR(mst.vx, mst.vy);
                    if (fv >= line){ break; }
                    Simulator::invStepSolo<pfpn_t>(&mst, step);
                    tm += step;
                }
                tmst = mst;
                lmst = mst;
                Simulator::stepSolo<pfpn_t>(&lmst, step);
                pfpn_t lfv = XYtoR(lmst.vx, lmst.vy);
                pfpn_t hasu;
                if (fv != lfv){
                    hasu = (line - lfv) / (fv - lfv);// lfv < line <= fv
                }
                else{
                    hasu = static_cast<pfpn_t>(0.5);
                }
                ASSERT(0 <= hasu && hasu <= 1.0, cerr << double(lfv) << "->" << double(fv) << "(" << double(line) << ")" << endl;);
                pfpn_t fr = XYtoR(lmst.x, lmst.y)*((pfpn_t)1 - hasu) + XYtoR(tmst.x, tmst.y)*hasu;
                
                sim2int_t ir=FRtoI32R(fr);
                IV_IRTable[i] = ir;
            }
        }
        void initR_V(){//R_V変換初期化
            constexpr pfpn_t step = static_cast<pfpn_t>(1 / 700000.0);
            
            fMStone<> mst(0, 0, 0, 0, 0), lmst, tmst;
            IR_IVTable[0] = 0;
            pfpn_t tm = 0;//逆算時刻
            for (int i = 0; i<(1 << IR_IV_TABLE_RESOLUTION) + 4; ++i){
                //V_Rテーブルに現時点での値を入れていく
                pfpn_t line = I64RtoFR(IR_IV_TABLE_STEP * i);
                pfpn_t ir;
                while (1){
                    ir = XYtoR(mst.x, mst.y);
                    if (ir >= line){ break; }
                    Simulator::invStepSolo<pfpn_t>(&mst, step);
                    tm += step;
                }
                tmst = mst;
                lmst = mst;
                Simulator::stepSolo(&lmst, step);
                pfpn_t lir = XYtoR(lmst.x, lmst.y);
                pfpn_t hasu;
                if (ir != lir){
                    hasu = (line - lir) / (ir - lir);
                }
                else{
                    hasu = static_cast<pfpn_t>(0.5);
                }
                ASSERT(0 <= hasu && hasu <= 1.0, cerr << double(lir) << "->" << double(ir) << "(" << double(line) << ")" << endl;);
                pfpn_t fv=XYtoR(lmst.vx, lmst.vy)*(1 - hasu) + XYtoR(tmst.vx, tmst.vy)*hasu;
                sim2int_t iv=FVtoI32V(fv);
                IR_IVTable[i] = iv;
                //CERR<<line<<" -> "<<FRtoI64R(line)<<"  "<<fv<<" -> "<<iv<<endl;
            }
        }
        
        //V_ALL変換初期化
        void initALL(){
            constexpr pfpn_t step = static_cast<pfpn_t>(1 / 700000.0);
            
            fMStone<pfpn_t> mst(0, 0, 0, 0, 0), lmst, tmst;
            V_AllIntTable[0].set(0, 0, 0, 0);
            int tm = 0;//逆算時刻
            for (int i = 1; i<(1 << IV_TABLE_RESOLUTION) + 4; ++i){
                pfpn_t line = I32VtoFV(IDtoIV(i));
                pfpn_t iv;
                while (1){
                    iv = XYtoR(mst.vx, mst.vy);
                    if (iv >= line){ break; }
                    Simulator::invStepSolo<pfpn_t>(&mst, step);
                    tm += step;
                }
                if (i % 64 == 0){
                    //cerr<<XYtoT(-mst.y,-mst.x)-XYtoT(mst.vy,mst.vx)<<endl;
                }
                tmst = mst;
                lmst = mst;
                Simulator::stepSolo<pfpn_t>(&lmst, step);
                pfpn_t liv = XYtoR(lmst.vx, lmst.vy);
                pfpn_t hasu;
                if (iv != liv){
                    hasu = (line - liv) / (iv - liv);
                }
                else{
                    hasu = static_cast<pfpn_t>(0.5);
                }
                ASSERT(0 <= hasu && hasu <= 1.0, cerr << double(liv) << "->" << double(iv) << "(" << double(line) << ")" << endl;);
                //進行方向基準に回転
                //CERR<<double(tmst.x)<<","<<double(tmst.y)<<endl;
                //rotate( tmst.y, tmst.x, -XYtoT(tmst.vy,tmst.vx), &tmst.y, &tmst.x );
                //CERR<<" -> "<<tmst.x<<","<<tmst.y<<endl;
                
                //rotate( lmst.y, lmst.x, -XYtoT(lmst.vy,lmst.vx), &lmst.y, &lmst.x );
                V_AllIntTable[i].set(
                                     FRtoI64R(-lmst.x*(1 - hasu) + -tmst.x*hasu),
                                     FRtoI64R(-lmst.y*(1 - hasu) + -tmst.y*hasu),
                                     simint_t(cosl(-XYtoT(lmst.vy, lmst.vx)*(1 - hasu) + -XYtoT(tmst.vy, tmst.vx)*hasu) *(1<<16)),
                                     simint_t(sinl(-XYtoT(lmst.vy, lmst.vx)*(1 - hasu) + -XYtoT(tmst.vy, tmst.vx)*hasu) *(1<<16))
                                     );
            }
        }
        
        void init(){
            
            
            thread th0(initV_R);
            thread th1(initR_V);
            thread th2(initALL);
            
            th0.join();
            th1.join();
            th2.join();
        }
        /*
         template<class mst_t, class st_t>
         void collisionAw_As(mst_t *const aw, st_t *const as){
         //衝突処理(動-静)
         
         //衝突面がx軸と平行になるように回転
         simint_t dy=as->y - aw->y;
         simint_t dx=as->x - aw->x;
         simint_t dr=sqrt(XYtoR2(dy,dx));
         int colcos=(dy<<16)/dr;
         int colsin=(dx<<16)/dr;
         
         //CERR<<coltheta<<endl;
         
         int vn(aw->v * cos(aw->theta - coltheta)), vt(aw->v * sin(aw->theta - coltheta));
         
         //力積を質量で割ったもの
         int n_m(vn);
         int f_m((1.0 - I_SLID_STONES) / 6.0 * (IR_STONE_RAD * aw->w + vt));
         
         //角速度変化量
         int dw(f_m * 2.0 / IR_STONE_RAD *1.056);
         
         as->v = XYtoR(n_m, f_m);
         as->theta = XYtoT(n_m, f_m) + coltheta;
         as->w = dw;
         
         //CERR<<aw->w<<" ,"<<dw<<endl;
         
         as->t = aw->t;
         
         aw->v = XYtoR(vn - n_m, vt - f_m);
         aw->theta = XYtoT(vn - n_m, vt - f_m) + coltheta;
         aw->w += dw;
         }
         template<class mst_t>
         void collision2Aw(mst_t *const aw0, mst_t *const aw1){
         //衝突処理(動-動)
         
         //衝突面がx軸と平行になるように回転
         int coltheta(XYtoT(aw1->y - aw0->y, aw1->x - aw0->x));
         //CERR<<coltheta<<endl;
         
         int v0n(aw0->v * cos(aw0->theta - coltheta));
         int v0t(aw0->v * sin(aw0->theta - coltheta));
         
         int v1n(aw1->v * cos(aw1->theta - coltheta));
         int v1t(aw1->v * sin(aw1->theta - coltheta));
         
         //力積を質量で割ったもの
         int n_m(v0n - v1n);
         int f_m((1.0 - I_SLID_STONES) / 6.0 * (IR_STONE_RAD * (aw0->w + aw1->w) + (v0t - v1t)));
         
         //角速度変化量
         int dw(f_m * 2.0 / IR_STONE_RAD * 1.04);
         
         aw1->v = min(IV_MAX, XYtoR(v1n + n_m, v1t + f_m));
         
         ASSERT(aw1->v <= IV_MAX, cerr << aw1->v << endl;);
         
         aw1->theta = XYtoT(v1n + n_m, v1t + f_m) + coltheta;
         aw1->w += dw;
         
         //CERR<<dw<<endl;
         
         aw0->v = min(IV_MAX, XYtoR(v0n - n_m, v0t - f_m));
         
         ASSERT(aw0->v <= IV_MAX, cerr << aw0->v << endl;);
         
         aw0->theta = XYtoT(v0n - n_m, v0t - f_m) + coltheta;
         aw0->w += dw;
         }
         */
        /*
         template<class pos_t, class mv_t>
         void IXYtoFMV(const pos_t& pos, mv_t *const mv){
         //posで静止するVX,VYを計算
         int dx(pos.getX() - IX_THROW);
         int dy(pos.getY() - IY_THROW);
         int r(XYtoR(dx, dy));
         int dtheta(XYtoT(dy, dx));
         int v(IRtoIV(r));
         
         int dtx, vtheta;//テーブルでの変化量
         int id = IVtoID(v);
         if (!mv->isLeitSpin()){
         dtx = +VID_AllIntTable[id].x;
         vtheta = +VID_AllIntTable[id].theta;
         }
         else{
         dtx = -VID_AllIntTable[id].x;
         vtheta = -VID_AllIntTable[id].theta;
         }
         
         int dttheta(XYtoT(VID_AllIntTable[id].y, dtx));
         
         rotate(v, 0,
         -dttheta + dtheta - vtheta,
         &mv->y, &mv->x
         );
         }
         */
        template<class mst_t>
        void step(mst_t *const dst,
                  const uint32 nid, const uint32 nfrac ){
            //前計算データから石の位置をステップさせる
            //id...現在状態のインデックス
            //nid...次状態のインデックス
            /*int
             dx(VID_AllIntTable[id].x),
             dy(VID_AllIntTable[id].y),
             dt(VID_AllIntTable[id].t - VID_AllIntTable[nid].t),
             dvtheta(VID_AllIntTable[id].theta - VID_AllIntTable[nid].theta),
             dtheta(VID_AllIntTable[id].theta);
             
             dx -= VID_AllIntTable[nid].x;
             dy -= VID_AllIntTable[nid].y;
             
             if (!(dst->w > 0)){
             dx = -dx;
             dvtheta = -dvtheta;
             dtheta = -dtheta;
             }
             rotateToAdd(
             dy, dx,
             dst->theta + dtheta,
             &dst->y, &dst->x
             );//位置を進める
             //dst->theta += dvtheta;//向き設定
             dst->t += dt;//時刻設定
             dst->v -= dt * I_FRIC_STONE;
             */
            
            //before rotation
            sim2int_t x = ID_FRACtoX(nid, nfrac);
            sim2int_t y = ID_FRACtoY(nid, nfrac);
            
            dst->x = dst->gx; dst->y = dst->gy;
            
            rotateToAdd(
                        y, x,
                        dst->cos, dst->sin,
                        dst->x, dst->y);
            
            dst->t = dst->gt-ID_FRACtoT(nid,nfrac);
            dst->r = ID_FRACtoR(nid, nfrac);
        }
        /*
         template<class mst_t>
         void stepToStop(mst_t *const dst,
         const uint32 id, const uint32 hasu){
         //前計算データから石の位置をステップさせる
         //衝突せずに静止するまで動かす
         //id...現在状態のインデックス
         int dx, dtheta;
         if (dst->w > 0){
         dx = VID_AllIntTable[id].x;
         dtheta = VID_AllIntTable[id].theta;
         }
         else{
         dx = -VID_AllIntTable[id].x;
         dtheta = -VID_AllIntTable[id].theta;
         }
         rotateToAdd(
         VID_AllIntTable[id].y, dx,
         dst->theta + dtheta,
         &dst->y, &dst->x
         );//座標設定
         dst->t += VID_AllIntTable[id].t;//時刻設定
         }
         */
        /*
         
         template<int RINK_ONLY = 1, class board_t, class board2_t>
         int simulateI_b2d_ISB(board_t *const pbd, board2_t *const tobd, const int step)
         {
         //シミュレータ、box2dへの入り口
         //Awakeな石が複数あっていい場合
         
         b2Board b2bd;
         
         //ゲーム情報の初期化
         BitSet32 tmp;
         //awakeな石をセット
         tmp = pbd->awake;
         for (; tmp.any(); tmp.pop_lsb()){
         int id = tmp.bsf();
         b2bd.setAwakeStone(pbd->stone[id], pbd->stone[id].state.getColor());
         }
         //asleepな石をセット
         tmp = pbd->asleep;
         for (; tmp.any(); tmp.pop_lsb()){
         int id = tmp.bsf();
         b2bd.setAsleepStone<RINK_ONLY>(pbd->stone[id], pbd->stone[id].state);
         }
         
         //物理演算
         int ret = Simulator::MainLoop<RINK_ONLY>(&b2bd, step, -1);
         
         assert(ret != -1);
         
         // シミュレーション結果をまとめる
         tobd->initStonePos();
         const int size = b2bd.NAwake + b2bd.NAsleep;
         
         //assert( size == pbd->NStones );
         
         DERR << size;
         
         for (int i = 0; i<size; ++i){
         assert(b2bd.body[i] != nullptr);
         const int col = b2bd.state[i].getColor();
         const b2Vec2 vec = b2bd.body[i]->GetPosition();
         if ((!RINK_ONLY) || isInPlayArea(vec.x, vec.y)){
         tobd->pushStone(col, vec.x, vec.y);
         //DERR<<"("<<vec.x<<","<<vec.y<<")"<<endl;
         }
         }
         return ret;
         }
         template<int RINK_ONLY = 1, class board_t, class move_t>
         int simulateSoloF(board_t *const pbd, const move_t& mv){
         //ただ停止位置を計算し、石を置くだけ
         //元々石がある場合にはおかしなことになりうるので、注意
         const int col = pbd->getTurnColor();
         //fPosXY<> dst;
         //calcDestination( mv, &dst );
         int v(XYtoR(mv.getVX(), mv.getVY()));
         int vtheta(XYtoT(mv.getVY(), mv.getVX()));
         int x(IX_THROW), y(IY_THROW), dx, dtheta;
         
         int id = IVtoID(v);
         if (!mv.isLeitSpin()){
         dx = +VID_AllIntTable[id].x;
         dtheta = +VID_AllIntTable[id].theta;
         }
         else{
         dx = -VID_AllIntTable[id].x;
         dtheta = -VID_AllIntTable[id].theta;
         }
         rotateToAdd(VID_AllIntTable[id].y, dx,
         vtheta + dtheta,
         &y, &x
         );
         
         if ((!RINK_ONLY) || isInPlayArea(x, y)){
         pbd->pushStoneWithState(col, x, y);
         return 1;
         }
         else{
         return 0;
         }
         }
         
         template<int RINK_ONLY = 1, class board_t>
         int simulateSoloF(board_t *const pbd, const int v, const int vtheta, const int spin){
         //ただ停止位置を計算し、石を置くだけ
         //元々石がある場合にはおかしなことになりうるので、注意
         const int col = pbd->getTurnColor();
         int x(IX_THROW), y(IY_THROW), dx, dtheta;
         
         int id = IVtoID(v);
         if (spin == Spin::RIGHT){
         dx = +VID_AllIntTable[id].x;
         dtheta = +VID_AllIntTable[id].theta;
         }
         else{
         dx = -VID_AllIntTable[id].x;
         dtheta = -VID_AllIntTable[id].theta;
         }
         rotateToAdd(VID_AllIntTable[id].y, dx,
         vtheta + dtheta,
         &y, &x
         );
         
         if ((!RINK_ONLY) || isInPlayArea(x, y)){
         pbd->pushStoneWithState(col, x, y);
         return 1;
         }
         else{
         return 0;
         }
         }
         */
        /*
         template<int RINK_ONLY = 1, class board_t, class pos_t>
         int loopHitZone1_1(board_t *const pbd,sim2int_t gx,sim2int_t gy){
         
         //ヒットゾーンでのシミュレーション
         //awake x 1、asleep x 1の場合のみ
         ASSERT(pbd->NActive[0] + pbd->NActive[1] == 1, cerr << pbd->NActive[0] << "," << pbd->NActive[1] << endl;);
         
         const int aw_col = pbd->getTurnColor();
         const int as_col = pbd->NActive[BLACK] ? BLACK : WHITE;
         
         sim2int_t asx(FRtoI64R(accessStone(*pbd, as_col, 0).x)), asy(FRtoI64R(accessStone(*pbd, as_col, 0).y));
         
         IMovingStone aw;
         
         aw.x = IX_THROW; aw.y = IX_THROW;
         aw.gx = gx; aw.gy = gy;
         aw.r=FRtoI64R(XYtoR(I64RtoFR(aw.gx-aw.x),I64RtoFR(aw.gy-aw.y)));
         
         uint32 id = IRtoID(aw.r);
         uint32 frac = IRtoFRAC(aw.r);
         
         int r2_hantei, r2_asleep;
         
         r2_asleep = XYtoR2(aw.x - as.x, aw.y - as.y);
         
         for (;;){
         //r2_hantei=( r_goal + 2*IR_STONE_RAD )*( r_goal + 2*IR_STONE_RAD );
         
         //まず非衝突判定
         if (r2_asleep < pow(IR_COLLISION + 2 * IR_STONE_RAD, 2)){//衝突
         //衝突処理
         //CERR<<"collision!"<<endl;
         
         collisionAw_As(&aw, &as);
         
         //CERR<<" AW "<<aw.x<<","<<aw.y<<","<<aw.v*sin(aw.theta)<<","<<aw.v*cos(aw.theta)<<","<<aw.w<<endl;
         //CERR<<" AS "<<as.x<<","<<as.y<<","<<as.v*sin(as.theta)<<","<<as.v*cos(as.theta)<<","<<as.w<<endl;
         
         //1_1では1度衝突したら2度と衝突することはない
         stepToStop(IVtoID(aw.v), &aw);
         stepToStop(IVtoID(as.v), &as);
         //フリーガード違反判定
         if (pbd->state[as_col][0].isFreeGuard()){
         if (RINK_ONLY && !isInPlayArea(as.x, as.y)){//違反
         DERR << "ireeguard foul 1-1" << endl;
         return 0;
         }
         pbd->initStonePos();
         pbd->pushStone(as_col, as.x, as.y);
         }
         else{
         pbd->initStonePos();
         if ((!RINK_ONLY) || isInPlayArea(as.x, as.y)){
         pbd->pushStone(as_col, as.x, as.y);
         }
         }
         if ((!RINK_ONLY) || isInPlayArea(aw.x, aw.y)){
         pbd->pushStone(aw_col, aw.x, aw.y);
         }
         return 2;
         }
         else{//1ステップ進める
         int r_asleep = sqrt(r2_asleep) - 2 * IR_STONE_RAD;
         
         if (r_goal < r_asleep){
         goto UNCOL;
         }
         
         int nr = sqrt(r_goal*r_goal + r_asleep*r_asleep - 2.0*r_goal*r_asleep*cos(ITHETA_CURL_MAX));
         
         uint32 nid = IRtoID(nr);
         
         //CERR<<aw.x<<","<<aw.y<<","<<aw.v*sin(aw.theta)<<","<<aw.v*cos(aw.theta)<<","<<aw.t;
         step(id, nid, &aw);
         //CERR<<" -> "<<aw.x<<","<<aw.y<<","<<aw.v*sin(aw.theta)<<","<<aw.v*cos(aw.theta)<<","<<aw.t<<endl;
         
         id = nid;
         r_goal = nr;
         }
         
         int r2_asleep_last = r2_asleep;
         r2_asleep = XYtoR2(aw.x - as.x, aw.y - as.y);
         if (r2_asleep > r2_asleep_last + 0.00001){//到達しない
         goto UNCOL;
         }
         }
         UNCOL:
         stepToStop(id, &aw);
         if ((!RINK_ONLY) || isInPlayArea(aw.x, aw.y)){
         pbd->pushStoneWithState(aw_col, aw.x, aw.y);
         return 1;
         }
         else{
         return 0;
         }
         assert(0);
         }
         
         template<int RINK_ONLY = 1, class board_t, class mst_t>
         int loopHitZone1_1(board_t *const pbd, const mst_t& arg){
         int v_goal = min(IV_MAX, XYtoR(arg.vx, arg.vy));
         int r_goal = IVtoIR(v_goal);
         int v_theta = XYtoT(arg.getVY(), arg.getVX());
         
         //CERR<<v_goal<<" "<<v_theta<<" "<<r_goal<<endl;
         
         return loopHitZone1_1<RINK_ONLY>(pbd, fPosXY<>(arg.getX(), arg.getY()), arg.getW(), v_goal, v_theta, r_goal);
         }
         
         template<int RINK_ONLY = 1, class board_t, class fsbd_t>
         int loopHitZone(board_t *const pabd, fsbd_t *const pbd){
         //ヒットゾーンでのシミュレーション(ISBoardによる)
         
         int& NAwake = pbd->NAwake;
         int& NAsleep = pbd->NAsleep;
         BitSet32& awake = pbd->awake;
         BitSet32& asleep = pbd->asleep;
         FSMovingStone *const stone = pbd->stone;
         
         int loop;
         for (loop = 0;; ++loop){
         
         LSTART:
         assert(pbd->exam());
         
         if (loop > 100){
         //ハマったのでここでBox2Dに投げる
         DERR << "to box2D" << endl;
         static StaticCounter c;
         c.add();
         simulateI_b2d_ISB<RINK_ONLY>(pbd, pabd, g_timeStep);
         return 2;
         }
         
         bool coll = false;//衝突が起きうるかどうか
         
         int minCollT = 9999;//衝突しうる再短時間
         
         //awakeなそれぞれの石に対して、asleepな石に対しての移動可能距離を求める
         for (BitSet32 tmp = awake; tmp.any(); tmp.pop_lsb()){
         int aw = tmp.bsf();
         
         //asleep_stones
         const int r_goal = IVtoIR(stone[aw].v);
         const int r2_hantei = (r_goal + 2 * IR_STONE_RAD)*(r_goal + 2 * IR_STONE_RAD);
         int r2_asleep = 9999;
         //int as_near;
         for (BitSet32 tas = stone[aw].asleep; tas.any(); tas.pop_lsb()){
         int last = tas.bsf();
         int r2 = XYtoR2(stone[last].x - stone[aw].x, stone[last].y - stone[aw].y);
         if (r2 < r2_hantei){
         if (r2 < pow(IR_COLLISION + 2 * IR_STONE_RAD, 2)){
         DERR << "collision AW_AS " << endl;
         collisionAw_As(&stone[aw], &stone[last]);//衝突
         
         asleep.reset(last);
         
         stone[aw].asleep = asleep;
         stone[last].asleep = asleep;
         
         ++NAwake; --NAsleep;
         
         //この2つ以外の石のリストを更新
         BitSet32 taw = awake;
         taw.reset(aw);
         
         stone[aw].awake = taw;
         stone[last].awake = taw;
         
         for (; taw.any(); taw.pop_lsb()){
         int an_aw = taw.bsf();
         stone[an_aw].awake.set(aw);
         stone[an_aw].awake.set(last);
         stone[an_aw].asleep.reset(last);
         }
         
         awake.set(last);
         ++loop;
         
         assert(pbd->exam());
         
         goto LSTART;
         }
         if (r2 < r2_asleep){
         r2_asleep = r2;
         //as_near=last;
         }
         }
         else{
         //非衝突が証明された
         stone[aw].asleep.reset(last);
         }
         }
         
         stone[aw].t_goal = IVtoIT(stone[aw].v);
         //CERR<<stone[aw].t_goal<<endl;
         stone[aw].r = r_goal;
         
         if (r2_asleep < r2_hantei){//静止まで動かせない
         
         int r_asleep = sqrt(r2_asleep) - 2 * IR_STONE_RAD;
         
         int nr = sqrt(r_goal*r_goal + r_asleep*r_asleep - 2.0*r_goal*r_asleep*cos(ITHETA_CURL_MAX));
         int mt = stone[aw].t_goal - IRtoIT(nr);
         
         minCollT = min(minCollT, mt);
         
         coll = true;
         }
         
         //awake_stones
         if (stone[aw].awake){
         //自分以下の番号の石のみ抽出
         BitSet32 taw = stone[aw].awake & ((1U << aw) - 1U);
         
         //CERR<<taw.count()<<endl;
         
         for (; taw.any(); taw.pop_lsb()){
         int last = taw.bsf();
         
         int r2 = XYtoR2(stone[aw].x - stone[last].x, stone[aw].y - stone[last].y);
         
         if (r2 < pow(IR_COLLISION + 2 * IR_STONE_RAD, 2)){
         //動-動衝突
         DERR << "collision 2AW" << endl;
         collision2Aw(&stone[aw], &stone[last]);//衝突
         //動-動関係を整理
         if (NAwake > 2){
         BitSet32 ttaw = awake;
         ttaw.reset(aw); ttaw.reset(last);
         
         stone[aw].awake = ttaw;
         stone[last].awake = ttaw;
         
         for (; ttaw.any(); ttaw.pop_lsb()){
         int an_aw = ttaw.bsf();
         stone[an_aw].awake.set(aw);
         stone[an_aw].awake.set(last);
         }
         }
         else{
         stone[aw].awake.reset();
         stone[last].awake.reset();
         }
         ++loop;
         
         goto LSTART;
         }
         
         int r_awake = sqrt(r2) - 2 * IR_STONE_RAD;
         
         if (r_awake < r_goal + stone[last].r){//静止まで動かせない
         
         //減速せずに互いに一直線に近づくときの衝突までの時間を計算
         int t = r_awake / (stone[aw].v + stone[last].v);
         
         minCollT = min(minCollT, t);
         coll = true;
         }
         else{//非衝突が証明された
         stone[aw].awake.reset(last);
         stone[last].awake.reset(aw);
         }
         }
         }
         }
         
         //全ての石が衝突せずに静止する事が証明されたので、全ての石を置いて終わり
         if (!coll){
         DERR << "no collision" << endl;
         for (BitSet32 taw = awake; taw.any(); taw.pop_lsb()){
         int aw = taw.bsf();
         stepToStop(IVtoID(stone[aw].v), &stone[aw]);
         if ((!RINK_ONLY) || isInPlayArea(stone[aw].x, stone[aw].y)){//エリア内
         asleep.set(aw);
         }
         else{//エリア外
         if (stone[aw].state.isFreeGuard()){//フリーガード違反
         DERR << "ireeguard foul" << endl;
         return 0;
         }
         }
         }
         goto END;
         }
         
         DERR << "Min Time : " << minCollT << endl;
         
         //それぞれの石を移動可能な分動かす
         {
         for (BitSet32 taw = awake; taw.any(); taw.pop_lsb()){
         int aw = taw.bsf();
         
         DERR << "(" << aw << ") " << stone[aw].x << "," << stone[aw].y << ","
         << stone[aw].v*sin(stone[aw].theta) << ","
         << stone[aw].v*cos(stone[aw].theta) << "," << stone[aw].t;
         
         if (stone[aw].t_goal < minCollT){
         
         //最短衝突時刻より停止時刻が早い石は静止させる
         
         stepToStop(IVtoID(stone[aw].v), &stone[aw]);
         
         DERR << " -> " << stone[aw].x << "," << stone[aw].y << ","
         << stone[aw].v*sin(stone[aw].theta) << ","
         << stone[aw].v*cos(stone[aw].theta) << "," << stone[aw].t << " stop" << endl;
         
         awake.reset(aw);
         --NAwake;
         
         if ((!RINK_ONLY) || isInPlayArea(stone[aw].x, stone[aw].y)){//エリア内
         asleep.set(aw);
         ++NAsleep;
         for (BitSet32 ttaw = awake; ttaw.any(); ttaw.pop_lsb()){
         int an_aw = ttaw.bsf();
         stone[an_aw].awake.reset(aw);
         stone[an_aw].asleep.set(aw);
         }
         }
         else{//エリア外
         if (stone[aw].state.isFreeGuard()){//フリーガード違反
         DERR << "ireeguard foul" << endl;
         return 0;
         }
         for (BitSet32 ttaw = awake; ttaw.any(); ttaw.pop_lsb()){
         int an_aw = ttaw.bsf();
         stone[an_aw].awake.reset(aw);
         }
         }
         }
         else{
         int nid = IVtoID(stone[aw].v - ITtoIV(minCollT));
         if (nid>0){
         //静止しない
         step(IVtoID(stone[aw].v), nid, &stone[aw]);
         
         DERR << " -> " << stone[aw].x << "," << stone[aw].y << ","
         << stone[aw].v*sin(stone[aw].theta) << ","
         << stone[aw].v*cos(stone[aw].theta) << "," << stone[aw].t << endl;
         
         if (RINK_ONLY && !isInPlayArea(stone[aw].x, stone[aw].y)){//エリア外
         if (stone[aw].state.isFreeGuard()){//フリーガード違反
         DERR << "ireeguard foul" << endl;
         return 0;
         }
         awake.reset(aw);
         --NAwake;
         for (BitSet32 ttaw = awake; ttaw.any(); ttaw.pop_lsb()){
         int an_aw = ttaw.bsf();
         stone[an_aw].awake.reset(aw);
         }
         }
         }
         else{//丸め込みにより静止
         stepToStop(IVtoID(stone[aw].v), &stone[aw]);
         
         DERR << " -> " << stone[aw].x << "," << stone[aw].y << ","
         << stone[aw].v*sin(stone[aw].theta) << ","
         << stone[aw].v*cos(stone[aw].theta) << "," << stone[aw].t << " stop" << endl;
         
         awake.reset(aw);
         --NAwake;
         
         if ((!RINK_ONLY) || isInPlayArea(stone[aw].x, stone[aw].y)){//エリア内
         asleep.set(aw);
         ++NAsleep;
         for (BitSet32 ttaw = awake; ttaw.any(); ttaw.pop_lsb()){
         int an_aw = ttaw.bsf();
         stone[an_aw].awake.reset(aw);
         stone[an_aw].asleep.set(aw);
         }
         }
         else{//エリア外
         if (stone[aw].state.isFreeGuard()){//フリーガード違反
         DERR << "ireeguard foul" << endl;
         return 0;
         }
         for (BitSet32 ttaw = awake; ttaw.any(); ttaw.pop_lsb()){
         int an_aw = ttaw.bsf();
         stone[an_aw].awake.reset(aw);
         }
         }
         }
         }
         }
         }
         
         if (NAwake == 0){
         break;
         }
         }
         END:
         //プレイアウト用の局面構造体に石を置く
         pabd->initStonePos();
         for (BitSet32 tmp = asleep; tmp.any(); tmp.pop_lsb()){
         int id = tmp.bsf();
         DERR << id << " ";
         pabd->pushStone(stone[id].state.getColor(), stone[id].x, stone[id].y);
         }
         DERR << endl;
         
         DERR << loop << endl;
         
         return 2;
         }
         
         template<int RINK_ONLY = 1, class board_t, class pos_t>
         int loopHitZone1_N(board_t *const pbd, const pos_t& arg, const int w, const int v_goal, const int v_theta, int r_goal){
         
         //ヒットゾーンでのシミュレーション
         struct AsInfo{
         int i;
         int r2_last;
         };
         
         AsInfo asleepList[16];
         const int aw_col = pbd->getTurnColor();
         
         int r2_hantei = (r_goal + 2 * IR_STONE_RAD)*(r_goal + 2 * IR_STONE_RAD);
         int r2_asleep = 9999;
         int as_near;
         
         //asleepList作成
         int dvtheta;
         if (w > 0){
         dvtheta = +ITHETA_CURL_MAX / 2;
         }
         else{
         dvtheta = -ITHETA_CURL_MAX / 2;
         }
         int vtheta_ex = v_theta + dvtheta;
         int NAsleep = 0;
         for (int c = 0; c<2; ++c){
         for (int n = 0; n<pbd->NActive[c]; ++n){
         //衝突可能性判定
         int dx = pbd->stone[c][n].getX() - arg.x;
         int dy = pbd->stone[c][n].getY() - arg.y;
         int r2 = XYtoR2(dx, dy);
         if (r2 < r2_hantei){
         int theta = XYtoT(dy, dx);
         //if( fabs( theta-vtheta_ex ) <= ITHETA_CURL_MAX/2 ){
         int as_num = c * 8 + n;
         asleepList[NAsleep].i = as_num;
         asleepList[NAsleep].r2_last = r2;
         ++NAsleep;
         if (r2<r2_asleep){
         r2_asleep = r2;
         as_near = as_num;
         }
         //}
         }
         }
         }
         
         //cerr<<NAsleep<<endl;
         
         FMoment aw;
         aw.set(arg.x, arg.y, v_theta, 0, v_goal, w);
         uint32 id = IVtoID(v_goal);
         
         for (;;){
         if (r2_asleep > r2_hantei){//到達しない
         stepToStop(id, &aw);
         if ((!RINK_ONLY) || isInPlayArea(aw.x, aw.y)){
         pbd->pushStoneWithState(aw_col, aw.x, aw.y);
         return 1;
         }
         else{
         return 0;
         }
         }
         else if (r2_asleep < pow(IR_COLLISION + 2 * IR_STONE_RAD, 2)){//衝突
         FMoment as;
         int c = as_near / 8;
         int n = as_near % 8;
         as.x = pbd->stone[c][n].x;
         as.y = pbd->stone[c][n].y;
         collisionAw_As(&aw, &as);
         ISBoard fsb;
         fsb.setAiterFirstColl(*pbd, as_near, aw, as);
         return loopHitZone<RINK_ONLY>(pbd, &fsb);
         
         //return simulateI_b2d<RINK_ONLY>( pbd, aw, g_timeStep );
         }
         else{//1ステップ進める
         int r_asleep = sqrt(r2_asleep) - 2 * IR_STONE_RAD;
         
         //int nr=r_goal - r_asleep * 0.8;
         int nr = sqrt(r_goal*r_goal + r_asleep*r_asleep - 2.0*r_goal*r_asleep*cos(ITHETA_CURL_MAX));
         uint32 nid = IRtoID(nr);
         //CERR<<aw.x<<","<<aw.y<<","<<aw.v*sin(aw.theta)<<","<<aw.v*cos(aw.theta)<<","<<aw.t;
         step(id, nid, &aw);
         //CERR<<" -> "<<aw.x<<","<<aw.y<<","<<aw.v*sin(aw.theta)<<","<<aw.v*cos(aw.theta)<<","<<aw.t<<endl;
         id = nid;
         r_goal = nr;
         }
         //最も近距離の石を探す
         r2_hantei = (r_goal + 2 * IR_STONE_RAD)*(r_goal + 2 * IR_STONE_RAD);
         r2_asleep = 9999;
         vtheta_ex = aw.theta + dvtheta;
         for (int i = NAsleep - 1; i >= 0; --i){
         int c = asleepList[i].i / 8;
         int n = asleepList[i].i % 8;
         int dx = pbd->stone[c][n].getX() - aw.x;
         int dy = pbd->stone[c][n].getY() - aw.y;
         int r2 = XYtoR2(dx, dy);
         if (r2<r2_hantei && r2 < asleepList[i].r2_last + 0.00001){
         
         int theta = XYtoT(dy, dx);
         //if( fabs( theta-vtheta_ex ) <= ITHETA_CURL_MAX/2 ){
         //衝突可能性あり
         //ASSERT( !((sqrt(r2)-2*IR_STONE_RAD < IR_COLLISION) && (r2>=asleepList[i].r2_last)),
         //	   cerr<<r2<<" , "<<asleepList[i].r2_last<<endl;);
         
         asleepList[i].r2_last = r2;
         if (r2<r2_asleep){
         r2_asleep = r2;
         as_near = asleepList[i].i;
         }
         continue;
         //}
         }
         //非衝突が証明された
         asleepList[i] = asleepList[--NAsleep];//リストを詰める
         }
         }
         assert(0);
         }
         
         template<int RINK_ONLY = 1, class board_t, class mst_t>
         int loopHitZone1_N(board_t *const pbd, const mst_t& arg){
         int v_goal = min(IV_MAX, XYtoR(arg.vx, arg.vy));
         int r_goal = IVtoIR(v_goal);
         int v_theta = XYtoT(arg.getVY(), arg.getVX());
         return loopHitZone1_N<RINK_ONLY>(pbd, fPosXY<>(arg.getX(), arg.getY()), arg.getW(), v_goal, v_theta, r_goal);
         }
         
         template<class board_t, class move_t>
         int simulateF(board_t *const pbd, const move_t& mv)
         {
         //高速シミュレータ
         int ret;
         
         //シミュレーション開始座標と速度を求める
         
         //停止距離計算
         int v_goal = min(IV_MAX, XYtoR(mv.getVX(), mv.getVY()));
         int r_goal = IVtoIR(v_goal);
         int v_theta = XYtoT(mv.getVY(), mv.getVX());
         if (r_goal < pbd->irontY - IR_STONE_RAD - IY_THROW){
         //どの石にも届かず停止する
         return simulateSoloF(pbd, v_goal, v_theta, mv.getSpin());
         }
         else{
         int w = mv.isLeitSpin() ? (-FANGV_ORG) : FANGV_ORG;
         if (pbd->NActive[BLACK] + pbd->NActive[WHITE] == 1){
         //tick();
         ret = loopHitZone1_1<1>(pbd, fPosXY<>(IX_THROW, IY_THROW), w, v_goal, v_theta, r_goal);
         //tock();
         }
         else{
         //tick();
         ret = loopHitZone1_N<1>(pbd, fPosXY<>(IX_THROW, IY_THROW), w, v_goal, v_theta, r_goal);
         //tock();
         }
         }
         return ret;
         }
         
         template<int SPIN, class pos_t>
         void IV_IRtoIXY(int v, int r, pos_t *const pos, int *const theta){
         //初速vのショットが投擲点から距離rの点での相対位置を求める
         //もしrまで到達しないならば、適当な値が帰るので注意
         int iv = min(v, IV_MAX);
         int id = IVtoID(iv);
         int gx, gy, x, y;
         if (SPIN == Spin::RIGHT){
         *theta = -VID_AllIntTable[id].theta;
         gx = VID_AllIntTable[id].x;
         gy = VID_AllIntTable[id].y;
         }
         else{
         *theta = VID_AllIntTable[id].theta;
         gx = -VID_AllIntTable[id].x;
         gy = VID_AllIntTable[id].y;
         }
         
         //到達点から距離ir-rの点を探す
         int ir = IVtoIR(iv);
         if (ir - r > IR_COLLISION){
         
         int r2_coll = r*(r + 2 * IR_COLLISION);
         int cnt = 0;
         
         int tr = ir - r;
         int nid = IRtoID(tr);
         
         while (1){
         //CERR<<*pos<<endl;
         //int nv=IRtoIV()
         //CERR<<ir<<" -> ";
         //この点のところまで位置を移動
         if (SPIN == Spin::RIGHT){
         x = gx - VID_AllIntTable[nid].x;
         y = gy - VID_AllIntTable[nid].y;
         ir = XYtoR2(x, y);
         }
         else{
         x = gx + VID_AllIntTable[nid].x;
         y = gy - VID_AllIntTable[nid].y;
         ir = XYtoR2(x, y);
         }
         //CERR<<sqrt(ir)<<" ( "<<r<<" )"<<" nid "<<nid<<endl;
         ++cnt;
         
         if (ir < r2_coll){
         //ASSERT( ir>r*r, cerr<<sqrt(ir)<<" ( "<<r<<" )"<<endl; );
         break;
         }
         ir = sqrt(ir);
         int nnid;
         do{
         tr += ir - r;
         //CERR<<tr<<endl;
         nnid = IRtoID(tr);
         } while (nnid == nid);
         nid = nnid;
         }
         //CERR<<cnt<<endl;
         }
         
         pos->setX(x); pos->setY(y);
         //CERR<<" ( "<<x<<" , "<<y<<" )"<<endl;
         }
         
         template<class pos_t, class move_t>
         void rotateToPassPointF(const pos_t& pos, int v, move_t *const mv){
         //指定速さで、指定された位置を通る速度ベクトルを生成
         
         //tick();
         
         const fPosXY<> tar(pos.x - IX_THROW, pos.y - IY_THROW);
         const int r = XYtoR(tar.getX(), tar.getY());
         
         fPosXY<> tp;
         int theta;
         
         if (mv->isLeitSpin()){
         IV_IRtoIXY<Spin::LEIT>(v, r, &tp, &theta);
         }
         else{
         IV_IRtoIXY<Spin::RIGHT>(v, r, &tp, &theta);
         }
         int dtheta = XYtoT(tar.y, tar.x) - XYtoT(tp.y, tp.x) + theta;
         rotate(v, 0, dtheta, &mv->y, &mv->x);
         
         //tock();
         //CERR<<mv->x<<" , "<<mv->y<<endl;
         
         }
         
         template<class pos_t, class move_t>
         void rotateToPassPointGoF(const pos_t& pos, int r, move_t *const mv){
         //指定された位置を通り、さらに指定距離分(投擲点から見て)進む速度ベクトルを生成
         const fPosXY<> tar(pos.x - IX_THROW, pos.y - IY_THROW);
         const int pr = XYtoR(tar.getX(), tar.getY());
         const int v = IRtoIV(pr + r);
         
         fPosXY<> tp;
         int theta;
         
         if (mv->isLeitSpin()){
         IV_IRtoIXY<Spin::LEIT>(v, pr, &tp, &theta);
         }
         else{
         IV_IRtoIXY<Spin::RIGHT>(v, pr, &tp, &theta);
         }
         int dtheta = XYtoT(tar.y, tar.x) - XYtoT(tp.y, tp.x) + theta;
         rotate(v, 0, dtheta, &mv->y, &mv->x);
         
         //CERR<<mv->x<<" , "<<mv->y<<endl;
         }
         */
        
    }
}

#endif // DCURLING_SIMULATION_INTSIMULATOR_H_