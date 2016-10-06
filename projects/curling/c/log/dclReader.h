//dcl形式のログファイル読み込み
//20150913時点のログ

#include <dirent.h>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "dc.hpp"

using namespace std;
using boost::lexical_cast;
using boost::algorithm::split;

string DIRECTORY_PARAM(""),DIRECTORY_DCLOG(""),DIRECTORY_SHOTLOG("");

namespace DigitalCurling{
    namespace DCL{

#define Get if(!ifs){return pLog;}getline(ifs,str,'\n');\
if(str[str.size()-1]=='\n'){str.erase(str.size()-1);}\
if(str[str.size()-1]=='\r'){str.erase(str.size()-1);}
        
#define ToV v.clear();split(v, str, boost::algorithm::is_space(), boost::algorithm::token_compress_on);
#define ToF(str) (atof((str).c_str()))
#define ToI(str) lexical_cast<long long>(str)
   /*
        template<class log_t,class callback_t>
        int readToGameLog(log_t *const pLog,const std::string& fName,const callback_t& callback){
            
            ifstream ifs(fName);
            string str;
            vector<string> v;
            
            Get;
            pLog->setName(BLACK,trim(trim(str,'='),'\n'));
            Get;
            pLog->setTimeLimit(BLACK,trim(trim(str,'='),'\n'));
            Get;
            pLog->setName(WHITE,ToI(trim(trim(str,'='),'\n')));
            Get;
            pLog->setTimeLimit(WHITE,ToI(trim(trim(str,'='),'\n')));
            Get;
            pLog->setRandom(ToF(trim(trim(str,'='),'\n')));
            
            int ends = 10;
            int turns = 16;

            for(int e=0;e<=ends;++e){
                for(int t=0;t<turns;++t){

                    Get;
                    Get; ToV;//POSITION
                    for (int i = 0; i < 16; ++i){
                        pLog->setStonePosition(e, t, i, ToF(v[1 + i * 2]), ToF(v[1 + i * 2 + 1]));
                        pLog->setAfter
                    }
                    Get;//SETSTATE
                    Get;ToV//BESTSHOT
                    pLog->setChosenShot(e, t, i, ToF(v[1]), ToF(v[1 + 1]), ToI(v[1+2]));
                    Get; ToV;//RUNSHOT
                    pLog->setRunShot(e, t, i, ToF(v[1]), ToF(v[1 + 1]), ToI(v[1+2]));
                }

                Get; ToV;//POSITION (after last shot)
                for (int i = 0; i < 16; ++i){
                    pLog->setStonePosition(e, turns, i, ToF(v[1 + i * 2]), ToF(v[1 + i * 2 + 1]));
                }
                Get; ToV;//SCORE
                int sc = ToI(v[1 + i * 2]);
                pLog->setEndScore(e, sc);
            }
            Get; ToV;
            pLog->setTotalScore(ToI(v[1]), ToI(v[2]));
        }
        return 0;
	}
*/
        template<class callback_t>
        MoveLog* readDCLogToMoveLog(MoveLog *const pLog0, const string& fName){
            
            cerr<<fName<<endl;
            
            ifstream ifs(fName);
            
            if(!ifs){return pLog0;}
            
            MoveLog *pLog=pLog0;
            
            string str;
            vector<string> v;
            
            fPosXY<> officialPos;
            fPosXY<> ayumuPos;
            fMoveXY<> officialMove;
            fMoveXY<> ayumuMove;
            
            //何度も使う変数
            string name[2];
            fpn_t random;
            int reversed = 0;
            int relScore = 0;
            
            Get;
            Get;
            name[BLACK] = str.substr(str.find_first_of("=")+1,str.size());
            Get;
            //pLog->setTimeLimit(BLACK, trim(trim(str, '='), '\n'));
            Get;
            name[WHITE] = str.substr(str.find_first_of("=")+1,str.size());
            Get;
            //pLog->setTimeLimit(WHITE, ToI(trim(trim(str, '='), '\n')));
            Get;
            
            cerr<<str.substr(str.find_first_of("=")+1,str.size())<<endl;
            
            random = ToF(str.substr(str.find_first_of("=")+1,str.size()));
            
            int end_max = 9;
            int turn_max = TURN_FIRST;
            
            int ends = end_max - END_LAST + 1;
            int turns = turn_max - TURN_FIRST + 1;
            
            for (int e = end_max; e >= END_LAST; --e){
                for (int t = turn_max; t >= TURN_LAST; --t){
                    
                    pLog->player = name[reversed];
                    pLog->oppPlayer = name[1 ^ reversed];
                    
                    pLog->end = e;
                    pLog->turn = t;
                    pLog->relScore = relScore;
                    
                    Get;
                    Get; ToV;//POSITION
                    for (int i = 0; i < 16; ++i){
                        officialPos.set(ToF(v[1 + i * 2]), ToF(v[1 + i * 2 + 1]));
                        ayumuPos = convPosition_Official_Ayumu(officialPos);
                        if(!isInPlayArea(ayumuPos)){
                            ayumuPos=FPOSXY_THROW;
                        }
                        pLog->previousStone[i] = ayumuPos;
                        if (t != turn_max){
                            (pLog - 1)->afterStone[i] = ayumuPos;
                        }
                    }
                    Get;//SETSTATE
                    Get; ToV;//BESTSHOT
                    officialMove.set(ToF(v[1]), ToF(v[1 + 1]), ToI(v[1 + 2]));
                    ayumuMove = convMove_Official_Ayumu(officialMove);
                    pLog->chosenMove = ayumuMove;
                    Get; ToV;//RUNSHOT
                    officialMove.set(ToF(v[1]), ToF(v[1 + 1]), ToI(v[1 + 2]));
                    ayumuMove = convMove_Official_Ayumu(officialMove);
                    pLog->runMove = ayumuMove;
                    ++pLog;
                }
                
                Get;
                Get; ToV;//POSITION (after last shot)
                for (int i = 0; i < 16; ++i){
                    officialPos.set(ToF(v[1 + i * 2]), ToF(v[1 + i * 2 + 1]));
                    ayumuPos = convPosition_Official_Ayumu(officialPos);
                    if(!isInPlayArea(ayumuPos)){
                        ayumuPos=FPOSXY_THROW;
                    }
                    pLog->previousStone[i] = ayumuPos;
                    (pLog - 1)->afterStone[i] = ayumuPos;
                }
                Get; ToV;//SCORE
                int sc = ToI(v[1]);
                if (reversed){ sc = -sc; }
                
                for (auto *p = pLog - ends; p != pLog; ++p){
                    p->score = sc;
                }
                
                relScore += sc;
                reversed = (sc < 0);
                if (reversed){ relScore = -relScore; }
            }
            Get; ToV;//TOTAL SCORE
            return pLog;
        }
#undef ToI
#undef ToF
#undef ToV
#undef Get
        
        int makeMoveLog(){
            
            MoveLog shot[1024];
            
            vector<string> dclFileName;
            
            //get file name
            string path = DIRECTORY_DCLOG;
            
            //DERR<<path<<endl;
            
            DIR *pdir;
            dirent *pentry;
            
            pdir = opendir(path.c_str());
            if (pdir == nullptr){
                return -1;
            }
            do {
                pentry = readdir(pdir);
                if (pentry != nullptr){
                    dclFileName.emplace_back(path + string(pentry->d_name));
                }
            } while (pentry != nullptr);
            
            for(auto& s : dclFileName){
                cerr<<s<<" "<<endl;
            }
            
            ofstream ofs;
            ofs.open(DIRECTORY_SHOTLOG + "shotlog.dat",ios::out);
            ofs.close();
            
            for (int i = 0; i < dclFileName.size(); ++i){
                if(dclFileName.at(i).find(".dcl")!=string::npos){
                    MoveLog *last = readDCLogToMoveLog<MoveLog>(shot, dclFileName.at(i));
                    
                    for (MoveLog *pl=shot; pl!=last; ++pl){
                        ofstream ofs;
                        ofs.open(DIRECTORY_SHOTLOG + "shotlog.dat",ios::app);
                        ofs << pl->toString() << endl;
                        ofs.close();
                    }
                }
            }
            return 0;
        }
    }
}

int main(){
    
    using namespace DigitalCurling::DCL;
    
    ifstream ifs("config.txt");
    ifs>>DIRECTORY_PARAM;
    ifs>>DIRECTORY_DCLOG;
    ifs>>DIRECTORY_SHOTLOG;
    
    return makeMoveLog();
}
