
#ifndef SEIDEL_METHOD_H
#define SEIDEL_METHOD_H

#include <iostream>
using namespace std;

class seidel_method
{
    double **A;
    double *b;
    double *x;
    int n;
    double eps;
    long int  Nmax;
    void print_x();
    void GetStartX();
public:
    seidel_method(int n, double **a,double *b_, double eps_, long int Nmax_ );
    void solve();
    void print();
    double get_solve(int i);

};

#endif // SEIDEL_METHOD_H
