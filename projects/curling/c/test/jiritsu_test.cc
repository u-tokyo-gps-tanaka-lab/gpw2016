/*
 jiritsu_test.cc
 Katsuki Ohto
 */

// じりつくんを真似した部分のコードのテスト
#include "../simulation/fastSimulator.hpp"
#include "../jiritsu/search.hpp"

using namespace std;

namespace DigitalCurling{
    int testJiritsuGrids(){
        
        Jiritsu::GridBoard<Jiritsu::kDrawLayer> drawLayer;
        Jiritsu::GridBoard<Jiritsu::kTakeOutLayer> takeOutLayer;
        
        cerr << "draw layer" << endl;
        cerr << drawLayer.toInfoString();
        cerr << "take out layer" << endl;
        cerr << takeOutLayer.toInfoString();
        
        return 0;
    }
        
}

int main(int argc, char* argv[]){
    
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stdin, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
    
    using namespace DigitalCurling;
    std::random_device seed;
    
    //Tester::dice.srand(seed() * (unsigned int)time(NULL));
    //Tester::ddice.srand(seed() * (unsigned int)time(NULL));

    testJiritsuGrids();
    
    return 0;
}