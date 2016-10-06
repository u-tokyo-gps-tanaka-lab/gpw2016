// デジタルカーリング
// 方策関数を方策勾配を用いた教師付き学習

#include "ayumu_dc.hpp"

#include "../simulation/fastSimulator.hpp"
#include "../simulation/primaryShot.hpp"
#include "shot/allShots.hpp"

#include "./eval/stat.hpp"
#include "./eval/evaluator.hpp"

namespace DigitalCurling{
    Peel gPeel;
    DrawRaise gDrawRaise;
    Double gDouble;
    RaiseTakeOut gRaiseTakeOut;
    TakeOut gTakeOut;
    L1Draw gL1Draw;
    
    Evaluator gEvaluator;
}

#include "logic/logic.hpp"
#include "structure/thinkBoard.hpp"

#include "move/generator.hpp"
#include "move/realizer.hpp"

#include "learn/classifier.hpp"
#include "learn/PG.hpp"

#include "initialize.hpp"

constexpr int iterations = 50;

using namespace DigitalCurling;

PlayPolicy policy;
PlayPolicyLearner learner;

int learn(int mode, const std::string& lpath){
    
    using namespace DigitalCurling;
    
    double learnRate = (mode == 0) ? 1 : 0.8;
    
    std::mt19937 mt((uint32)time(NULL));
    
    policy.setLearner(&learner);
    learner.setPolicy(&policy);
    
    // ログをファイルから読み込み
    RandomAccessor<AyumuShotLog> db;
    
    cerr << "started reading shot-log file." << endl;
    //readShotLog(DIRECTORY_SHOTLOG + "shotlog_mini.dat", &db);
    int games = readShotLog(lpath, &db, [&](const auto& shot)->bool{
        if(shot.end == END_LAST && shot.turn < 8){ return false; } // 最終盤は得点差の影響が大きいので除外
        if(abs(shot.rscore) > min(5, (shot.end + 1) * 2)){ return false; } // 得点差が大きい場合は適当な手が多いので除外
        if(!(contains(shot.player(), "umu") || contains(shot.player(), "kun") || contains(shot.player(), "GCCS"))){ // 指定したプレーヤーのみ
            return false;
        }
        return true;
    });
    cerr << db.size() << " shots in " << games << " games were imported." << endl;
    
    // ログのそれぞれの手を分類
    cerr << "started classifying each shot." << endl;
    int dist[N_STANDARD_MOVES][2] = {0};
    iterate(db, [&](auto& shot)->void{
        ThinkBoard bd;
        bd.setShotLog(shot);
        
        gEvaluator.setRewardInBound(bd.getTurnColor(), bd.getTurnColor(), bd.getRelScore(), bd.getEnd(), bd.getTurn(),
                                    [](Color c, int e, int s, int nrs)->double{
                                        return getEndWP(c, e, s, nrs);
                                    });
        gEvaluator.normalize([](int s)->double{
            return dScoreProb_EASY[s + N_COLOR_STONES];
        });
        
        MoveXY mv = classifyShot(bd, shot.chosenMove, dist);
        //MoveXY mv = classifyShotVxVy(bd, shot.chosenMove, dist); // 方策が弱かった
        shot.vmv = mv; // 自分の分類を記録
    });
    // 分類の分布表示
    cerr << "classification result = " << endl;
    for(int i = 0; i < N_STANDARD_MOVES; ++i){
        cerr << shotString[i] << " : " << dist[i][1] << " in " << dist[i][0] << endl;
    }
    
    // 特徴要素を解析して、学習の度合いを決める
    learner.setLearnParam(1, 0.00005, 0, 0.0000002, 1);
    
    cerr << "started analyzing feature." << endl;
    
    learner.initFeatureValue();
    
    BitSet32 flag(0);
    flag.set(1);
    
    iterate(db, [&](const auto& shot)->void{
        PG::learnPlayParamsShot(shot, flag, policy, &learner); // 特徴要素解析
    });
    learner.closeFeatureValue();
    learner.foutFeatureSurvey(DIRECTORY_PARAMS_OUT + "policy_feature_survey.dat");
    learner.finFeatureSurvey(DIRECTORY_PARAMS_OUT + "policy_feature_survey.dat");
    
    for(int t = 0; t < N_TURNS; ++t){
        cerr << "record : t = " << ((t < 10) ? " " : "") << t
        << " " << learner.toRecordString(t) << endl;
    }
    
    // ログのアクセス順を並べ替え
    db.initOrder();
    db.shuffle(0, db.size(), mt);
    
    // 学習とテストのショット数決定
    int shots = db.size();
    int learnShots = min((int)(shots * learnRate), shots);
    int testShots = shots - learnShots;
    
    // 学習開始
    cerr << "started learning." << endl;
    
    for (int j = 0; j < iterations; ++j){
        
        cerr << "iteration " << j << endl;
        
        learner.setLearnParam(1, 0.00005 * (0.01 + 0.99 / sqrt(double(j + 1))), 0, 0.0000001 * (0.01 + 0.99 / sqrt(double(j + 1))), 256);
        
        db.shuffle(0, learnShots, mt);
        
        BitSet32 flag(0);
        flag.set(0);
        flag.set(2);
        
        learner.initObjValue();
        
        iterateRandomly(db, 0, learnShots, [&](const auto& shot)->void{
            PG::learnPlayParamsShot(shot, flag, policy, &learner); // 学習
        });
        
        policy.fout(DIRECTORY_PARAMS_OUT + "policy_param.dat");
        foutComment(policy, DIRECTORY_PARAMS_OUT + "policy_comment.txt");
        
        for(int t = 0; t < N_TURNS; ++t){
            cerr << "Learning : t = " << ((t < 10) ? " " : "") << t
            << " " <<learner.toObjValueString(t) << endl;
        }
        if(testShots > 0){
            flag.reset(0);
            learner.initObjValue();
            
            iterateRandomly(db, learnShots, shots, [&](const auto& shot)->void{
                PG::learnPlayParamsShot(shot, flag, policy, &learner); // テスト
            });
            for(int t = 0; t < N_TURNS; ++t){
                cerr << "Test     : t = " << ((t < 10) ? " " : "") << t
                << " " << learner.toObjValueString(t) << endl;
            }
        }
    }
    cerr << "finished learning." << endl;
    
    return 0;
    
}

int main(int argc, char* argv[]){
    
    //基本的な初期化
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stdin, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
    
    int ret;
    
    int mode = 0;
    std::string lpath = DIRECTORY_SHOTLOG + "shotlog.dat";
    
    for(int c = 1; c < argc; ++c){
        if(strstr(argv[c], "-l")){
            lpath = std::string(argv[c + 1]);
        }else if(strstr(argv[c], "-t")){
            mode = 1;
        }
    }
    initAyumu(DIRECTORY_PARAMS_IN);
    
    ret = learn(mode, lpath);
    
    return ret;
}