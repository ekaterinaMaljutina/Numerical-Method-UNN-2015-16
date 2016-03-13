#include <iostream>
#include <cstring>
#include "cmath"
#include "simply_iteration.h"
#include "Task.h"
#include "Task_with_cup.h"
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

double solve(double x, double y){
    return exp(x*y);
}

int main()
{

    int area;
    cout << " 0 - rectangle"<<endl \
         <<" 1 - rectangle with cut "<<endl;
    cin>>area;


    double a = -1.0 , b=0.0,
            c = 0.0,  d = 1.0,
            x_r = (a+b)/2, y_r = (c+d)/2;
    int size;
    double *x ,*y,
            *F;
    double **A;

    long int n,m;
    double h,k,a2;

    int numder_task;

    if (!area){
        //задача ствиться в прямоугольной области

        cout<< " 0-  test tesk "<<endl<<" 1 - basic task "<<endl;
        numder_task;
        cin>>numder_task;

        Task task(numder_task);


        cout<<" n - spet to x"<<endl<<" m - step to y"<<endl;
        cout<< "n = ";
        cin>>n;
        cout<<"m = ";
        cin>>m;


        //шаги по сетке
        h = (b-a)/n;
        k = (d-c)/m;

        //элементы в матрице
        double h2 = 1/(h*h);
        double k2 = 1/(k*k);
        a2 = -2 * (h2 + k2);
        size = (n-1)*(m-1);

        x = new double[n+1];
        memset(x,0,sizeof(x));
        y = new double [m+1];
        memset(y,0,sizeof(y));

        //запонляем значение узлов сетки
        for (int i=0;i<n+1;i++)
            x[i]= a + i*h;
        for (int i=0;i<m+1;i++)
            y[i] = c + i*k;

        F = new double [size];
        memset(F,0,sizeof(F));

        A = new double*[ size ];
        for (int i=0; i <size; i++ ){
            A[i] =new double [size];
            memset(A[i],0,sizeof(A[i]));
        }

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
          }
        int tmp=0;
        for (int i=0;i < m-1;i++){
            for (int j = 0; j< n-1 ;j++){
                F[tmp] = task.f(x[j+1],y[i+1]);
                if ( i==0 || i==m-2 ){
                    if (i == 0){
                        if (j==0)
                            F[tmp]  = F[tmp] - k2 * task.mu_1(x[j+1]) - h2*task.mu1(y[i+1]);
                        else
                               if (j==n-2)
                                   F[tmp] = F[tmp] - k2 * task.mu_1(x[j+1]) - h2*task.mu2(y[i+1]);
                               else
                                   F[tmp] = F[tmp]  - k2*task.mu_1(x[j+1]);

                    }

                    if (i==m-2){
                        if (j==0)
                            F[tmp]  = F[tmp] - k2 * task.mu_2(x[j+1]) - h2*task.mu1(y[i+1]);
                        else
                               if (j==n-2)
                                   F[tmp] = F[tmp] - k2 * task.mu_2(x[j+1]) - h2*task.mu2(y[i+1]);
                               else
                                   F[tmp] = F[tmp]  - k2*task.mu_2(x[j+1]);
                    }
                }
                else{
                    if (j==0)
                        F[tmp]  = F[tmp]- h2*task.mu1(y[i+1]);
                    else
                           if (j==n-2)
                               F[tmp] = F[tmp] - h2*task.mu2(y[i+1]);
                }
                tmp++;
            }
        }
        cout<<endl;
        //print(A,size);
        if (n<10)
            print(A,F,size);

        /* тестирование самого  метода
        double** A = new double* [3];
        for (int i=0;i<3;i++)
            A[i] = new double [3];
        A[0][0] = A[1][1] = A[1][2] = A[2][0] = A[2][2] = 1;
        A[0][1] = A[0][2] = A[2][1] = 0;
        A[1][0] = 2;

        double b[3] = {1,2,1};

        print(A,b,3);

        Simply_iteration sympl_iter_method(3,A,b,0.5,0.00001,30);
        cout<<endl;
        sympl_iter_method.solve();
        cout<<endl;
        */
    }
    else
    {

        //задача в прямоугольной области с вырезом

        cout<< " 0-  test tesk "<<endl<<" 1 - basic task "<<endl;

        cin>>numder_task;

        Task task(numder_task);

        bool flag_x = true,
                flag_y = true;
        //подбирем нужные нам сетки, чтоб попадали в границы вырезанной обдасти
        cout<<" n - spet to x"<<endl<<" m - step to y"<<endl;
        while(flag_x || flag_y){
                cout<< "n = ";
                cin>>n;
                cout<<"m = ";
                cin>>m;
                //шаги по сетке
                h = (b-a)/n;
                k = (d-c)/m;

                x = new double[n+1];
                memset(x,0,sizeof(x));
                y = new double [m+1];
                memset(y,0,sizeof(y));

                //запонляем значение узлов сетки
                for (int i=0;i<n+1;i++)
                    x[i]= a + i*h;
                for (int i=0;i<m+1;i++)
                    y[i] = c + i*k;
                for (int i=0;i<n+1;i++)
                    if (x[i] == x_r)
                        flag_x = false;
                for (int i=0;i<m+1;i++)
                    if (y[i] == y_r)
                        flag_y =false;
                if (flag_x || flag_y){
                    delete[] x;
                    delete[] y;
                }
        }


        size = (n-1)*(m/2 - 1) + (n/2 - 1)*(m/2 - 1) + (n/2 - 1);


        A = new double* [ size ];
        for (int i=0; i < size; i++ ){
            A[i] =new double [size];
            for (int j=0;j<size;j++){
                A[i][j]=0;
            }
        }
        //print(A,size);
        //элементы в матрице
        double h2 = 1/(h*h);
        double k2 = 1/(k*k);
        a2 = -2 * (h2 + k2);

        long int size_basic = (n-1)*(m/2 - 1) + (n/2 - 1);

        //cout<<endl<<size<<"   "<<size_basic<<endl;

        //Заполняем матрицу А
        for (int i=0; i< size_basic ;i++){
            A[i][i] = a2;
            if ( i - (n-1) >=0)
                A[i][i-(n-1)] = k2;
            if (i+ (n-1) < size_basic)
                A[i][i+(n-1)] =k2;
            if (i-1>=0 && i % (n-1) != 0)
                A[i][i-1]=h2;
            if (i+1< size_basic && i % (n-1) != n-2)
                A[i][i+1]= h2;
          }
        //print(A,size);
        int j=0;
        int t = n/2 - 1;
        int t1 = n/2 - 1;

        int size_dop = (n/2 - 1)*(m/2 - 1);
        for (int i=size_basic ; i< size ;i++){
            A[i][i] = a2;
            if (t>0){
                A[i - t - j][i] = k2;
                A[i][i -t - j] = k2;
                t--;
            }
            if ( i + t1  < size ){
                A[i + t1][i] = k2;
                A[i][i + t1] = k2;
            }
            j++;
          }
        t = n/2 - 1;
        t1 =m/2 - 1;
        for (int i=0; i< size_dop ;i++){
            if ( i % t == 0)
                A[i + size_basic][i + size_basic +1 ]=h2;
            else
                if ( i % t == (t-1)){
                    A[i + size_basic ][i + size_basic  -1 ]=h2;
                }

                else{
                    A[i + size_basic ][i + size_basic + 1 ]=h2;
                    A[i + size_basic ][i + size_basic  -1 ]=h2;
                }
        }
        //print(A,size);


        F = new double [size];
        memset(F,0,sizeof(F));

        for (int i=0; i<n+1;i++){
               cout<<x[i]<<" ";
        }
        cout<<endl;

        for (int i=0; i<n+1;i++){
               cout<<y[i]<<" ";
        }
        cout<<endl;

        int tmp=0;
        // 1 квадрат
        for (int i=0; i<t;i++){
            for (int j=0;j<t1;j++){
                F[j + tmp] = (-1) * task.f(x[j+1],y[i+1]);
                cout<<x[j+1]<<"  "<<y[i+1]<<endl;
                cout<<"  "<<F[j]<<endl;
                if (tmp == 0)
                    F[j + tmp] -= k2 *task.mu_1(x[j+1]);
                if (j == 0)
                    F[j + tmp] -= h2 * task.mu1(y[i+1]);

            }
            tmp += (2*t+1);
        }
        // 2 квадрат
        tmp = 0;
        for (int i=0; i<t1;i++){
            for (int j=t+2;j<n;j++){
                F[j-1 + tmp] = (-1) * task.f(x[j],y[i+1]);
                cout<<x[j]<<"  "<<y[i+1]<<endl;
                cout<<"  "<<F[j-1]<<endl;
                if ( tmp == 0 ){
                    F[j - 1 + tmp] -= k2 * task.mu_1(x[j]);
                }

                if ( j == n-1 && i!=t1-1){
                     F[j - 1 + tmp] -= h2 *task.mu2(y[i+1]);
                }

                if ( j == n-1 && i==t1-1){
                    F[j - 1 + tmp] -= h2 *task.mu2(y[i+1]);
                    F[j - 1 + tmp] -= k2 *task.cut(x[j],y[i+1]);
                }
                if ( i==t1-1 && j<n-1){
                    F[j - 1 + tmp] -= k2 *task.cut(x[j],y[i+1]);
                }
            }
            tmp += (2*t+1);

             // вертикальная грань выреза

             F[(2*i +1)*t + i] = (-1) * task.f(x[t+1],y[i+1]);
             cout<<(2*i +1)*t + i<<endl;
             if (i==0){
                 F[(2*i +1)*t + i] -= k2* task.mu_1(x[t+1]);
             }
             if (i == t1-1){
                 F[(2*i +1)*t + i] -= k2* task.cut(x[t+1],y[t+1]);
             }

             // горизонтальная грань
             F[(n-1) * t1 + i] = (-1) * task.f(x[i+1],y[t+1]);
             if (i ==0 ){
                 F[(n-1) * t1 + i] -= h2*task.mu1(y[t+1]);
             }
             if (i == t1-1)
                F[(n-1) * t1 + i] -= h2* task.cut(x[i+1],y[t+1]);
        }
        //3 квадрат
        tmp=(n-1)*(m/2 - 1) + (n/2-1);
        //cout<<tmp<<endl;
        for(int i=t1 + 2 ;i<m;i++){
            for (int j=0; j<t;j++){
                F[j+ tmp] = (-1) * task.f(x[j+1],y[i]);
                //cout<<x[j+1]<<"  "<<y[i]<<endl;
                //cout<<"  "<<F[j + tmp]<<endl;
                if (j==0){
                    //cout<<"   "<< j + tmp<<endl;
                    F[j + tmp] -= h2 *task.mu1(y[i]);
                }
                if ( j == t-1){
                    F[j + tmp] -= h2 *task.cut(x[j+1],y[i]);
                }
                if (i==m-1){
                    F[j + tmp] -= k2* task.mu_2(x[j+1]);
                }
            }
            tmp += t;
        }
        print(A,F,size);
    }


    delete x;
    delete y;

    //мы решаем систему Ах=в
    // но в А собственные числа отрицательные
    // поэтому решаем другую систему Су=в
    // где С = -А, у = -х

    //для этого исправим матрицу А (поменяем знаки)
    for (int i=0;i<size;i++)
        for (int j=0;j<size;j++)
            if (A[i][j]!=0)
                A[i][j] *= (-1);
    cout<<endl;
    if (n<10)
        print(A,F,size);

    //для применимости итерацианного метода решения СЛАУ
    // метода простой итерации с оптимальным параметром нужно знать собственные числа
    // в нашем случаем мы их не знаем
    // но мы можем взять оценку диапозона собст чисел
    // при этом все собственные числа матрицы С вещественные и положительные
    // поэтому по кругам Гиршгорина можно провести оценку
    // вычислим параметр метода простой итераци, используя оценку собственных чисел
    double lymda = 1/a2*(-1);
    cout<<"lymda = "<<lymda<<endl;

    // теперь можно применять метод простой итерации для подсчета решения

    double eps;
    cout<<"eps = ";
    cin>>eps;
    long int Nmax;
    cout<<"Nmax = ";
    cin>>Nmax;
    cout<< "wait solution"<<endl;

    Simply_iteration sympl_iter_method (size,A,F,lymda,eps,Nmax);
    sympl_iter_method.solve();

    // затем после подсчета вектра у, вычислим исходный х = -у
    if (n<10){
        cout<<endl<<"solve:"<<endl;
        for (int i=0;i<size;i++)
            cout<<sympl_iter_method.get_solve(i)*(-1)<<" ";
        cout<<endl;
    }

    if (numder_task ==0 ){
        double max = 0, tmp;
        int k=0;
        for (int i=0;i<(m/2-1)*(n-1);i++){
            if ( i%(n-1)==0)
                k++;
            tmp = abs( sympl_iter_method.get_solve(i)*(-1) - solve(x[i % (n-1) +1], y[k]) );
            if (tmp>max)
                max = tmp;
        }
        cout<<endl<<"max|u-v|="<<max<<endl;
    }


    for( int i = 0; i < size; i++ )
            delete [] A[i];
    delete [] A;
    delete F;

    return 0;
}

