#ifndef UP_RELACSATION_METHOD_H
#define UP_RELACSATION_METHOD_H
#include <iostream>
using namespace std;
class up_relacsation_method
{
    double **A;
    double *b;
    double *x;
    int n;
    double eps;
    long int  Nmax;
    double w;
    void GetStartX();
public:
    up_relacsation_method(int n, double **a,double *b,double eps_,long int  Nmax_,double w);
    void solve();
    void print_x();
};

#endif // UP_RELACSATION_METHOD_H
