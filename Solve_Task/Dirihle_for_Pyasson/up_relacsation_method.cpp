#include "up_relacsation_method.h"
#include "cmath"
up_relacsation_method::up_relacsation_method(int n, double **a, double *b_, double eps_, long Nmax_,double w): n(n), eps(eps_),Nmax(Nmax_),w(w) {
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
void up_relacsation_method::GetStartX(){
    for (int i=0;i<n;i++){
        x[i]=0;
    }
}

void up_relacsation_method::solve(){
    int S=0,i,j;
    double eps_max=0.0,eps_cur=0.0;
    double x_old, x_new;
    bool stop = false;
    while(!stop){
        eps_max =0.0;
        for(i=0;i<n;i++){
            x_old = x[i];
            x_new = w*b[i] + (1-w)*A[i][i]*x[i];
            for(j=0;j<n;j++)
                if (i!=j)
                    x_new = x_new -w*A[i][j]*x[j];
            x_new = x_new/A[i][i];
            eps_cur = fabs( x_old - x_new);
            if(eps_cur > eps_max)
                eps_max = eps_cur;
            x[i]=x_new;
        }
        print_x();
        S++;
        if  ((eps_max<eps) || (Nmax<=S))
            stop = true;
    }
    cout<< "method is stop"<<endl;
    cout<<"Nmax="<<Nmax<<"  eps="<<eps<<" behind="<<S<<endl;
    cout<<"solve: "<<endl;
    print_x();
    cout<< "with accuracy = "<<eps_max<<endl;

}

void up_relacsation_method::print_x(){
    for (int j=0;j<n;j++){
        cout<<x[j]<<" ";
    }
    cout<<endl;
}
