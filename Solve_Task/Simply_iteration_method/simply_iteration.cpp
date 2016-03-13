#include "simply_iteration.h"
#include "cmath"
Simply_iteration::Simply_iteration(int n, double **A_, double *f_, double lymda, double eps, double Nmax):
                    n(n), eps(eps), Nmax(Nmax), lymda(lymda)
{
    A = new double*[n];
        for(int i=0;i<n;i++){
            A[i] = new double[n];
        }
        f = new double[n];
        x = new double [n];
        for (int i=0; i<n; i++){
            f[i]=f_[i];
            for (int j=0; j<n; j++){
                A[i][j] = A_[i][j];
            }
        }
        GetStartX();
}

void Simply_iteration::print_x(){
    for (int j=0;j<n;j++){
            cout<<x[j]<<" ";
        }
        cout<<endl;
}

void Simply_iteration::solve(){
    double* result_vector = new double [n];


    double eps_max = 0,
            eps_cur = 0;
    int count =0;
    double x_new;
    bool flag = true;
    while ( flag ){
         eps_max=0;
        /*cout<<endl;
        print_x();
        cout<<endl;*/
        //считаем A*x[i]
        for(int i=0; i<n; i++)
            {
                result_vector[i]=0;
                for(int j=0; j<n; j++)
                {
                    //cout<<"A[i][j]*x[j] = "<<A[i][j]<<"*"<<x[j]<<"="<<A[i][j]*x[j]<<endl;
                    result_vector[i] += A[i][j]*x[j];
                    //cout<<"result["<<i<<"]="<<result_vector[i]<<endl;
                }
                //cout<<endl;
                //cout<<"result["<<i<<"]="<<result_vector[i]<<endl;
            }

        /*cout<<"A*x[i]:"<<endl;
        for (int j=0;j<n;j++){
                cout<<result_vector[j]<<" ";
            }
        cout<<endl;*/

        for (int i=0;i<n;i++){
            x_new = x[i] + lymda * (f[i] - result_vector[i] );
            //cout<<"x[i] - x_new =  "<<x[i]<<" - "<<x_new<<"  = "<<fabs(x[i] - x_new)<<endl;
            eps_cur = fabs(x[i] - x_new);
            if (eps_cur > eps_max)
                eps_max = eps_cur;
            x[i] = x_new;
        }
        //print_x();
        count++;
        if (eps_max < eps || count>=Nmax)
            flag = false;
    }

    cout<< "method is stop"<<endl;
    cout<<"Nmax="<<Nmax<<"  eps="<<eps<<" behind="<<count<<endl;
    /*cout<<"solve: "<<endl;
    if (n<30)
        print_x();*/
    cout<< "with accuracy = "<<eps_max<<endl;
}

void Simply_iteration::GetStartX(){
    //начальное приближение
    for (int i=0;i<n;i++)
            x[i]=0;

}


double Simply_iteration::get_solve(int i){
    return x[i];
}

void Simply_iteration::print(){
    for(int i=0;i<n;i++){
        for (int j=0;j<n;j++)
            cout<<A[i][j]<<" ";
        cout<<"    "<<f[i]<<endl;
    }
}
