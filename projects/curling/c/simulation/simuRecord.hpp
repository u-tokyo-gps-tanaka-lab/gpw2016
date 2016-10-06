//デジタルカーリング
//物理演算の記録

namespace DigitalCurling{
    
    struct CollisionRecord{
        fPosXY<> pos;//where collided
        int col[2];
    };
    
    struct SimuRecord{
        
        std::vector<CollisionRecord> vcol;//collision vector
        bool usingB2D;
    };
}

