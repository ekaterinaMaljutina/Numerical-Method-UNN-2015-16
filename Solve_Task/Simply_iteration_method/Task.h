#ifndef TASK_H
#define TASK_H

#include "cmath"

class Task{

public:
    Task(int num):numder_task(num){}

    double f(double x , double y){
        if (numder_task == 0)
            return y*y*exp(x*y) + x*x*exp(x*y);
        else
            return cosh(x*x*y);
    }
    double mu_1(double x){
        if (numder_task == 0)
            return 1;
        else
            return -x*(x+1);
    }
    double mu_2(double x){
        if (numder_task==0)
            return exp(x);
        else
            return -x*(x+1);
    }
    double mu1(double y){
        if (numder_task==0)
            return exp(-y);
        else
            return sin(M_PI * y);
    }
    double mu2(double y){
        if (numder_task ==0 )
            return 1;
        else
            return fabs( sin(M_PI*2*y));
    }


    double cut(double x, double y){
        return exp(x*y);
    }
private:
    int numder_task;
    // numder_task : 0 - test task
    //               1 - basic task
};


#endif // TASK_H
