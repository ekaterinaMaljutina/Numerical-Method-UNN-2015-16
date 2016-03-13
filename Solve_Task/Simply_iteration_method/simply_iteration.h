#ifndef SIMPLY_ITERATION_H
#define SIMPLY_ITERATION_H

#include <iostream>
using namespace std;

class Simply_iteration
{
private:
    double** A;
    double*  x;
    double*  f;
    int n;
    double eps;  
    long int  Nmax;
    double lymda;

    void print_x();
    void GetStartX();

public:
    Simply_iteration(int n, double** A_, double* f_,double lymda, double eps, double Nmax);
    void solve();
    double get_solve(int i);
    void print();
};

#endif // SIMPLY_ITERATION_H
