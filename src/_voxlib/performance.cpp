#include <iostream>

#include <cmath>
#include <chrono>
#include <iomanip>
#include "typses.h"

using namespace std;


class Timer
{
private:
    typedef chrono::high_resolution_clock clock_;
    typedef chrono::duration<double, ratio<1>> second_;
    chrono::time_point<clock_> beg_;
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {  return duration_cast<second_>(clock_::now() - beg_).count();  }
};
template<typename T> inline int expre(T expr) { return *reinterpret_cast<char*>(&expr)+15; }

int main(int argc, char** argv)  {
    Timer tmr;
    string indent="calibrating";
    double elpsPLus=1., elapsd;
    int total, totalOld=0;
    srand(tmr.elapsed());
    int rnd=rand();
#define  XPR(expr) expre(expr)

#define OP_TEST(typ, expr) \
    rnd=(rand()&7)+1;                         \
    totalOld = total = 0;                          \
    tmr.reset();                          \
    typ r1{}, r2{};                          \
    for (int i=0; i<500001; ++i) { \
        totalOld += bool(total); \
        rnd = (rnd+1)&7; \
        r1 = typ(i&127)*0.1+typ(0.05+1e-6*XPR(expr))+rnd+r2*(0.1+(i&1));     \
        r2 = typ(i&127)*0.1+typ(0.1+1e-6*XPR(expr))+r1*(0.1+(i&1));     \
        total +=  (XPR(expr)+XPR(expr))+(XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr));      \
        total +=  (XPR(expr)+XPR(expr))+(XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr));      \
        total -=  (XPR(expr)+XPR(expr))+(XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr));      \
        total +=  (XPR(expr)+XPR(expr))+(XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr));      \
        total *=  (XPR(expr)+XPR(expr))+(XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr)) +XPR(expr)+ (XPR(expr)+XPR(expr)+XPR(expr)+XPR(expr));      \
        total = (abs(total)&2047);    }   \
    elapsd = tmr.elapsed();          \
    std::cout<<std::setw(11)<<indent<< ": "<<std::setw(12)<<elapsd \
                 <<"    "<<std::setw(12)<<(elapsd)/elPsLoop<<"x  "<<std::setw(12)<<(elapsd-elPsLoop)/elpsPLus<<"x  "<<std::setw(10)<<total+totalOld<<"  "<<std::setw(12)<<#expr<<endl;

    // time the elpsPLus code:
    //   for loop with no elPsLoop math op
    double elPsLoop=1.;

    for (int ii=0;ii<5;++ii) {
      cout<<"calibrating to + operation, shall get 1++x           1.0x       for       r1+r2"<<endl;

      double basAvg=0., extraAvg=0.;
      for (int i=0;i<5;++i) {
        { OP_TEST(double, r2     ); extraAvg+=elapsd;}
        { OP_TEST(double, r1+r2);   basAvg+=elapsd; }
      }
      elPsLoop=extraAvg/5.;  elpsPLus=basAvg/5.-elPsLoop;
      cout<<"   basAvg: "<<basAvg/5.<<"  elPsLoop: "<<elPsLoop<<"  elpsPLus: "<<elpsPLus<<endl<<endl;
    }

    indent="";
    cout<< "\n  elPsLoop: "<<elPsLoop<<"  elpsPLus: "<<elpsPLus<<endl;

    // time various floating point operations.
    //   subtracts off the elPsLoop time to give
    //   a better approximation of the cost
    //   for just the specified operation
    cout<< "\n     indent:    elapsd        elps/elpsLoop   elps/elpsPLus     hash     expressiom "<<endl;
    indent="loop";
    { OP_TEST(double, 1.0                 ); }
    { OP_TEST(double, r1                  ); }
    { OP_TEST(double, r2                  ); }

    indent="";
    { OP_TEST(double, r1 + 10.1           ); }
    { OP_TEST(double, r1 + r2             ); }
    { OP_TEST(double, r1 + 0.1            ); }
    { OP_TEST(double, r1 + r2             ); }
    { OP_TEST(dbl3,   magSqr(r1*r1[0])    ); }
    { OP_TEST(dbl3,   magSqr(r1*r1.x)     ); }
    { OP_TEST(double, r1 + 10.1           ); }
    { OP_TEST(double, r2 + r1 + r2 + rnd  ); }
    { OP_TEST(double, r1 - r2             ); }
    { OP_TEST(double, r1 * r2             ); }
    { OP_TEST(double, r1 / r2             ); }
    //{ OP_TEST(int,    r1 % r2             ); }
    //{ OP_TEST(int,    r1 % r2             ); }
    { OP_TEST(double, fabs(r1)            ); }
    { OP_TEST(double, tanh(r1)            ); }
    { OP_TEST(double, atanh(r1)           ); }
    { OP_TEST(double, sqrt(r1+rnd)        ); }
    { OP_TEST(double, sqrt(r1+rnd)        ); }
    { OP_TEST(double, sqrt((sqrt(r1+rnd)))); }
    { OP_TEST(double, cbrt(r1)            ); }
    { OP_TEST(double, cbrt((cbrt(r1)))    ); }
    { OP_TEST(double, sin(r1)             ); }
    { OP_TEST(double, cos(r1)             ); }
    { OP_TEST(double, tan(r1)             ); }
    { OP_TEST(double, atan(r1)            ); }
    { OP_TEST(double, atan(r1)            ); }
    { OP_TEST(double, signbit(r1)         ); }
    { OP_TEST(double, exp(r1)             ); }
    { OP_TEST(double, log(r1)             ); }
    { OP_TEST(double, log10(r1)           ); }
    { OP_TEST(double, max(-1.*r1,r2)      ); }
    { OP_TEST(double, max(r2,r1)          ); }
    { OP_TEST(double, min(r1,r2)          ); }
    { OP_TEST(double, min(r2,r1)          ); }
    { OP_TEST(double, min(r2,max(r1,1.))  ); }
    { OP_TEST(double, min(r1,max(r2,r1))  ); }
    { OP_TEST(double, min(r2,max(r1,-1.)) ); }
    { OP_TEST(double, pow(r2,r1)          ); }
    { OP_TEST(double, pow(r1,0.2)         ); }
    { OP_TEST(double, pow(r1,5)           ); }
    { OP_TEST(double, pow(r1,2)           ); }
}
