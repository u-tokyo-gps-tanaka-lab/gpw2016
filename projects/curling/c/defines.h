//デジタルカーリング
//定義

#ifndef DCURLING_DEFINE_H_
#define DCURLING_DEFINE_H_

#define IP_SIZE 32
#define PORT_SIZE 16
#define ID_SIZE 32
#define PASS_SIZE 16
#define NAME_SIZE 32
#define BUFFER_SIZE 1024
#define MAX_THREAD_NUM 128


#ifdef _WIN32

using std::thread;

#else

using std::thread;

#define stricmp(s0,s1) strcasecmp(s0,s1)

#define Sleep(a) sleep(a) 

#define SOCKET_ERROR 1

typedef int SOCKET;
typedef unsigned long DWORD;
typedef pthread_t HANDLE;

typedef bool BOOL;

#define TRUE true
#define FALSE false

#endif

using std::string;
using std::array;
using std::vector;
using std::pair;
using std::tuple;
using std::swap;
using std::get;
using std::sort;



//引数による要求
namespace Requirement{
    constexpr int TEST_GAME=30;
    constexpr int PREPROCCESS=31;
}

#endif //DCURLING_DEFINE_H_