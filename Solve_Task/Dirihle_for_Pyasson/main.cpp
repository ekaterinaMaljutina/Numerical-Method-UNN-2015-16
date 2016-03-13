#include <iostream>
#include <cstring>
#include "Seidel_method.h"
#include "up_relacsation_method.h"
#include "cmath"
using namespace std;
void print(double **a,double *b,int n){
    for(int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            cout<<a[i][j]<<" ";
        }

        cout<<"    "<<b[i]<<endl;
    }
}
void print(double **a,int n){
    for(int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            cout<<a[i][j]<<" ";
        }

        cout<<endl;
    }
}
double f(double x , double y){
//return x*x + y*y -1;
//return 4;
    return y*y*exp(x*y) + x*x*exp(x*y);
}
double mu_1(double x){
//return -x*x;
//return 0;
    return 1;
}
double mu_2(double x){
//return -x*x;
//return 0;
    return exp(x);
}
double mu1(double y){
    return exp(-y);
//return 0.5*(y*y-1);
//return -y*y;
}
double mu2(double y){
    return 1;
//return 0.5*(y*y-1);
//return -y*y;
}

double solve(double x, double y){
    return exp(x*y);
}

int main()
{
    cout<< " 1- seidel method "<<endl<<" 2 - upper relaxation method"<<endl;
    int method_num;
    cin>>method_num;
    cout<<"Ok"<<endl;
    long int n,m;
    cout<<" n - spet to x"<<endl<<" m - step to y"<<endl;
    cout<< "n = ";
    cin>>n;
    cout<<"m = ";
    cin>>m;
    cout<<"Ok"<<endl;

    double a = -1.0 , b=0.0,
            c = 0.0,  d = 1.0;
    //long int n=5;
    //long int m=6;
    double h = (b-a)/n;
    double k = (d-c)/m;

    double h2 = 1/(h*h);
    double k2 = 1/(k*k);
    double a2 = -2 * (h2 + k2);
    long int size = (n-1)*(m-1);

    double *x = new double[n+1];
   // memset(x,0,sizeof(x));
    double *y = new double [m+1];
   // memset(y,0,sizeof(y));
    for (int i=0;i<n+1;i++)
        x[i]= a + i*h;
    for (int i=0;i<m+1;i++)
        y[i] = c + i*k;
    double *F = new double [size];
   // memset(F,0,sizeof(F));
    double **A = new double*[ size ];
    for (int i=0; i <size; i++ ){
        A[i] =new double [size];
      //  memset(A[i],0,sizeof(A[i]));
    }
    int j=1;
    //заполняем А
    for (int i=0; i< size;i++){
        A[i][i] = a2;
        if ( i - (n-1) >=0)
            A[i][i-(n-1)] = k2;
        if (i+ (n-1) <size)
            A[i][i+(n-1)] =k2;
        if (i-1>=0 && i % (n-1) != 0)
            A[i][i-1]=h2;
        if (i+1<size && i % (n-1) != n-2)
            A[i][i+1]= h2;

        //заполняем вектор F
       /* if (i % (n-1) == 0)
            j++;
        cout<<x[i % (n-1)]<<"  "<<y[j]<<endl;
        F[i] = f(x[i % (n-1)],y[j]);
        cout<<F[i]<<endl;*/
    }
     //print(A,F,size);
    int tmp=0;
    for (int i=0;i < m-1;i++){
        for (int j = 0; j< n-1 ;j++){
            F[tmp] = f(x[j+1],y[i+1]);//f(x[i+1],y[j+1]);
            //cout<<"i="<<i+1<<" "<<"j="<<j+1<<" "<<F[tmp]<<endl;
            if ( i==0 || i==m-2 ){
                if (i == 0){
                    if (j==0)
                        F[tmp]  = F[tmp] - k2 * mu_1(x[j+1]) - h2*mu1(y[i+1]);
                    else
                           if (j==n-2)
                               F[tmp] = F[tmp] - k2 * mu_1(x[j+1]) - h2*mu2(y[i+1]);
                           else
                               F[tmp] = F[tmp]  - k2*mu_1(x[j+1]);

                }

                if (i==m-2){
                    if (j==0)
                        F[tmp]  = F[tmp] - k2 * mu_2(x[j+1]) - h2*mu1(y[i+1]);
                    else
                           if (j==n-2)
                               F[tmp] = F[tmp] - k2 * mu_2(x[j+1]) - h2*mu2(y[i+1]);
                           else
                               F[tmp] = F[tmp]  - k2*mu_2(x[j+1]);
                }
            }
            else{
                if (j==0)
                    F[tmp]  = F[tmp]- h2*mu1(y[i+1]);
                else
                       if (j==n-2)
                           F[tmp] = F[tmp] - h2*mu2(y[i+1]);
            }
            //cout<<"i="<<i+1<<" "<<"j="<<j+1<<" "<<F[tmp]<<endl;
            tmp++;
        }
    }
    //print(A,size);
    if (n<10)
        print(A,F,size);
    double eps;
    cout<<"eps = ";
    cin>>eps;
    long int Nmax;
    cout<<"Nmax = ";
    cin>>Nmax;
    cout<< "wait solution"<<endl;
    if (method_num == 1){
        seidel_method method (size,A,F,eps,Nmax);
        method.solve();
        double max = 0, tmp;
        int k=0;
        for (int i=0;i<size;i++){
            if ( i%(n-1)==0)
                k++;
            tmp = abs( method.get_solve(i) - solve(x[i % (n-1) +1], y[k]) );
            if (tmp>max)
                max = tmp;
            //cout<<"max|u-v|="<<max<<endl;
        }
        cout<<"max|u-v|="<<max<<endl;
    }
    if (method_num == 2){
        cout<<" upper relaxation method use param W, please ask it (0 to 2)"<<endl;
        cout<<"w = ";
        double w;
        cin>>w;
        up_relacsation_method method(size,A,F,eps,Nmax,w);
        method.solve();
    }
    cout<<endl;
    return 0;
}

