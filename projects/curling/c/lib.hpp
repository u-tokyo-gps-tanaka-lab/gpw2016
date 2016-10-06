// デジタルカーリング
// 公式のデータ構造
// by Yuuma Kitasei

#ifndef DCURLING_LIB_HPP_
#define DCURLING_LIB_HPP_

struct NETWORKINFO{
    BOOL IsConnected;
    
    char IP[IP_SIZE];
    unsigned short PORT;
    char ID[ID_SIZE];
    char PASS[PASS_SIZE];
    char NAME[NAME_SIZE];
    
    SOCKET dstSocket;
    HANDLE hThread;
    struct sockaddr_in dstAddr;
};

// 改行の削除
void DeleteNL(char *msg)
{
    char *p;
    p = msg;
    
    while( *p != 0x00 ){
        if( *p == '\n' || *p == '\r' ){
            *p = 0x00;
            break;
        }
        p++;
    }
    return;
}

// 引数の取得
bool GetArgument( char *lpResult, size_t numberOfElements, char *msg, int n )
{
    bool bRet = false;
    char *p, *q;
    
    if( msg != NULL ){
        p = msg;
        while( *p == ' ' || *p == '\r' || *p == '\n'){
            p++;
        }
        
        // ポインタを取得したい引数の先頭に合わせる
        for( int i=0; i<n; i++ ){
            while( *p != ' ' ){
                if( *p == 0x00 ){
                    return false;
                }
                p++;
            }
            while( *p == ' ' ){
                p++;
            }
        }
        
        // 取得したい引数をlpResultにコピーする
        q = strstr( p, " " );
        if( q == NULL ){
            if (strlen(p) < numberOfElements){
                strcpy( lpResult, p );
                bRet = true;
            }
            else{
                bRet = false;
            }
        }
        else{
            strncpy( lpResult, p, q-p );
            if( ( unsigned int )( q-p ) < numberOfElements ){
                lpResult[q-p] = 0x00;
            }
            bRet = true;
        }
    }
    
    return bRet;
}

struct GAMESTATE{
    int     ShotNum;        // 現在のショット数
    // ShotNum が n の場合、次に行うショットが n+1 投目になる
    
    int     CurEnd;         // 現在のエンド数
    int     LastEnd;        // 最終エンド数
    int     Score[10];      // 第1エンドから第10エンドまでのスコア
    bool    WhiteToMove;    // 手番の情報
    // WhiteToMove が 0 の場合次のショットを行うのは先手、WhiteToMove が 1 の場合次のショットを行うのは後手となる
    
    float   body[16][2];    // body[n] は n 投目のストーンの位置座標を表す
    // body[n][0] は n 投目のストーンの x 座標、body[n][1] は n 投目のストーンの y 座標を表す
    
};

#endif // DCURLING_LIB_HPP_
