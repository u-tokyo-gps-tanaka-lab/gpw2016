// 盤面の画像ログを作る

#include "../dc.hpp"

std::string DIRECTORY_SHOTLOG(""), DIRECTORY_IMAGELOG("");

namespace DigitalCurling{
    int makeImageLog(const std::string& slPath, const std::string& ilPath){
        std::vector<ShotLog> slv;
        std::array<ImageLog, 1024> ila;
        int nila = 0;
        std::ofstream ofs;
        
        if(readShotLog(slPath, &slv, [](const auto& sl)->bool{ return true; }) < 0){
            return -1;
        }
        ofs.open(ilPath, std::ios::out);
        ofs.close();
        
        for(const auto& sl : slv){
            ImageLog il;
            il.end = sl.end;
            il.turn = sl.turn;
            il.rscore = sl.rscore;
            il.escore = sl.escore;
            MiniBoard mbd;
            mbd.clearStones();
            for(int i = 0; i < N_STONES; ++i){
                if(isInPlayArea(sl.stone(i))){
                    mbd.locateNewStone(i, sl.stone(i));
                }
            }
            genBoardImage(mbd, il.img);
            ila[nila++] = il;
            if(nila % 1024 == 0){
                std::ofstream ofs;
                ofs.open(ilPath, std::ios::app);
                for(int i = 0; i < nila; ++i){
                    ofs << ila[i].toString() << endl;
                }
                nila = 0;
                ofs.close();
            }
        }
        ofs.open(ilPath, std::ios::app);
        for(int i = 0; i < nila; ++i){
            ofs << ila[i].toString() << endl;
        }
        ofs.close();
    }
}

int main(){
    
    using namespace DigitalCurling;
    
    std::ifstream ifs("config.txt");
    std::string tmp;
    ifs >> tmp >> tmp >> tmp;
    ifs >> DIRECTORY_SHOTLOG;
    ifs >> DIRECTORY_IMAGELOG;
    
    //std::string slPath = DIRECTORY_SHOTLOG + "shotlog_mini.dat";
    //std::string ilPath = DIRECTORY_IMAGELOG + "imagelog_mini.dat";
    std::string slPath = DIRECTORY_SHOTLOG + "shotlog.dat";
    std::string ilPath = DIRECTORY_IMAGELOG + "imagelog.dat";
    
    return makeImageLog(slPath, ilPath);
}
