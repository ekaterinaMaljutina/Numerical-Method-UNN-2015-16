#include "Seidel_method.h"
#include "cmath"

seidel_method::seidel_method(int n,double **a,double *b_, double eps_, long int Nmax_ ): n(n), eps(eps_), Nmax(Nmax_){
    A = new double*[n];
    for(int i=0;i<n;i++){
        A[i] = new double[n];
    }
    b = new double[n];
    x = new double [n];
    for (int i=0; i<n; i++){
        b[i]=b_[i];
        for (int j=0; j<n; j++){
            A[i][j] = a[i][j];
        }
    }
    GetStartX();
}
void seidel_method::GetStartX(){
    for (int i=0;i<n;i++){
        x[i]=0;
    }
}

void seidel_method::solve(){
    int count=0;
    bool stop = false;
    double eps_max=0,eps_cur = 0;
    double x_old,x_new;
    while(!stop){
        eps_max=0;
        for (int i=0;i<n;i++){
            x_old = x[i];
            x_new = b[i];
            for (int j=0;j<n;j++){
                if (i!=j)
                    x_new = x_new - A[i][j]*x[j];
            }
            x_new = x_new/A[i][i];
            eps_cur = fabs(x_old - x_new);
            if (eps_cur > eps_max)
                eps_max = eps_cur;
            x[i] = x_new;
        }
        count++;
        //print_x();
        if ((eps_max<eps) || (count>=Nmax))
                stop = true;
        }

    cout<< "method is stop"<<endl;
    cout<<"Nmax="<<Nmax<<"  eps="<<eps<<" behind="<<count<<endl;
    cout<<"solve: "<<endl;
    if (n<10)
        print_x();
    cout<< "with accuracy = "<<eps_max<<endl;
}

void seidel_method::print(){
    for(int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            cout<<A[i][j]<<" ";
        }
        cout<<"    "<<b[i]<<endl;
    }
}

void seidel_method::print_x(){
    for (int j=0;j<n;j++){
        cout<<x[j]<<" ";
    }
    cout<<endl;
}
double seidel_method::get_solve(int i){
    return x[i];
}
