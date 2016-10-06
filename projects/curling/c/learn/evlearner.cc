// デジタルカーリング
// 評価関数関数を方策勾配法で学習

#include "../ayumu/ayumu_dc.hpp"

#include "../simulation/fastSimulator.hpp"
#include "../ayumu/simulation/primaryShot.hpp"
#include "../ayumu/shot/allShots.hpp"

#include "../ayumu/eval/stat.hpp"
#include "../ayumu/eval/evaluator.hpp"

#include "../ayumu/structure/thinkBoard.hpp"

#include "scorePG.hpp"

#include "../ayumu/initialize.hpp"

constexpr int iterations = 50;

using namespace DigitalCurling;

ScoreEstimator estimator;
ScoreEstimatorLearner learner;

int learn(const std::string& lpath, int mode){
    
    using namespace DigitalCurling;
    
    double learnRate = (mode == 0) ? 1 : 0.9;
    
    std::mt19937 mt((uint32)time(NULL));
    
    estimator.setLearner(&learner);
    learner.setPolicy(&estimator);
    
    // ログをファイルから読み込み
    //RandomAccessor<AyumuShotLog> db;
    
    cerr << "started reading shot-log file." << endl;
    
    ShotLog *const log = new Log[N];
    FILE *const pf = fopen(lpath.c_str(), "rb");
    fread(log, sizeof(Log) * N, 1, pf);
    fclose(pf);
    
    cerr << db.size() << " shots in " << games << " games were imported." << endl;
    /*
    uint64 scoreDist[N_TURNS][SCORE_LENGTH][SCORE_LENGTH];
    uint64 score2ndDist[N_TURNS][SCORE_LENGTH][SCORE_LENGTH];
    
    for(int t = 0; t < N_TURNS; ++t){
        for(int i0 = 0; i0 < SCORE_LENGTH; ++i0){
            for(int i1 = 0; i1 < SCORE_LENGTH; ++i1){
                scoreDist[t][i0][i1] = 1;
                score2ndDist[t][i0][i1] = 1;
            }
        }
    }
    
    iterate(db, [&scoreDist, &score2ndDist](const auto& log)->void{
        //cerr << log.score;
        ThinkBoard bd;
        bd.setShotLog(log);
        //bool clean = isCleanBetterScore(log.end, )
        ++scoreDist[bd.getTurn()][StoIDX(countScore(bd.orderBits))][StoIDX(log.escore)];
        ++score2ndDist[bd.getTurn()][StoIDX(count2ndScore(bd.orderBits, flipColor(bd.getTurnColor())))][StoIDX(log.escore)];
    });
    
    cerr << "tsc -> sc distribution" << endl;
    for(int t = 0; t < N_TURNS; ++t){
        cerr << "turn = " << t << endl;
        for(int s = 0; s < SCORE_LENGTH; ++s){
            for(int ss = 0; ss < SCORE_LENGTH; ++ss){
                cerr << std::setw(5) << scoreDist[t][s][ss] << " ";
            }cerr << endl;
        }cerr << endl;
    }
    cerr << "t2sc -> sc distribution" << endl;
    for(int t = 0; t < N_TURNS; ++t){
        cerr << "turn = " << t << endl;
        for(int s = 0; s < SCORE_LENGTH; ++s){
            for(int ss = 0; ss < SCORE_LENGTH; ++ss){
                cerr << std::setw(5) << score2ndDist[t][s][ss] << " ";
            }cerr << endl;
        }cerr << endl;
    }
    
    {
        std::ofstream ofs(DIRECTORY_PARAMS_OUT + "tscore_to_escore.dat");
        for(int t = 0; t < N_TURNS; ++t){
            for(int i0 = 0; i0 < SCORE_LENGTH; ++i0){
                for(int i1 = 0; i1 < SCORE_LENGTH; ++i1){
                    ofs << log(scoreDist[t][i0][i1]) << " ";
                }
            }
        }
    }
    {
        std::ofstream ofs(DIRECTORY_PARAMS_OUT + "to2score_to_escore.dat");
        for(int t = 0; t < N_TURNS; ++t){
            for(int i0 = 0; i0 < SCORE_LENGTH; ++i0){
                for(int i1 = 0; i1 < SCORE_LENGTH; ++i1){
                    ofs << log(score2ndDist[t][i0][i1]) << " ";
                }
            }
        }
    }
    exit(1);*/
    
    // 特徴要素を解析して、学習の度合いを決める
    learner.setLearnParam(1, 0.00005, 0, 0.0000002, 1);
    
    cerr << "started analyzing feature." << endl;
    
    learner.initFeatureValue();
    
    BitSet32 flag(0);
    flag.set(1);
    
    iterate(db, [&](const auto& shot)->void{
        PG::learnEvalParamsShot(shot, flag, estimator, &learner); // 特徴要素解析
    });
    learner.closeFeatureValue();
    learner.foutFeatureSurvey(DIRECTORY_PARAMS_OUT + "estimator_feature_survey.dat");
    learner.finFeatureSurvey(DIRECTORY_PARAMS_OUT + "estimator_feature_survey.dat");
    
    // ログのアクセス順を並べ替え
    //db.initOrder();
    //db.shuffle(0, db.size(), mt);
    
    // 学習とテストのショット数決定
    int shots = db.size();
    int learnShots = min((int)(shots * learnRate), shots);
    int testShots = shots - learnShots;
    
    // 学習開始
    cerr << "started learning." << endl;
    
    for (int j = 0; j < iterations; ++j){
        
        cerr << "iteration " << j << endl;
        
        learner.setLearnParam(1, 0.00005 * (0.01 + 0.99 / sqrt(double(j + 1))), 0, 0.0000001 * (0.01 + 0.99 / sqrt(double(j + 1))), 256);
        
        //db.shuffle(0, learnShots, mt);
        
        BitSet32 flag(0);
        flag.set(0);
        flag.set(2);
        
        learner.initObjValue();
        
        iterateRandomly(db, 0, learnShots, [&](const auto& shot)->void{
            PG::learnEvalParamsShot(shot, flag, estimator, &learner); // 学習
        });
        
        estimator.fout(DIRECTORY_PARAMS_OUT + "estimator_param.dat");
        foutComment(estimator, DIRECTORY_PARAMS_OUT + "estimator_comment.txt");
        
        for(int t = 0; t < N_TURNS; ++t){
            cerr << "Learning : t = " << ((t < 10) ? " " : "") << t
            << " " <<learner.toObjValueString(t) << endl;
        }
        
        if(testShots > 0){
            flag.reset(0);
            learner.initObjValue();
            
            iterateRandomly(db, learnShots, shots, [&](const auto& shot)->void{
                PG::learnEvalParamsShot(shot, flag, estimator, &learner); // テスト
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
    for(int c = 1; c < argc; ++c){
        if(strstr(argv[c], "-t")){
            mode = 1;
        }
    }
    
    initAyumu(DIRECTORY_PARAMS_IN);
    
    ret = learn(mode);
    
    return ret;
}