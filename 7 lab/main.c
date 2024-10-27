#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

// 7 1 1 ГУ 2
// y`` + y = 2x - pi [0,pi]
// точное решение y = 2x-pi+pi*cos(x) + sin(x)
// y(0) = A
// y(b) + y`(b) = B
// A = 0 
// B = 1

#define MAX_ITER 2000000

typedef struct FuncRes {
    double Point;
    double Func, Func1, Func2;
    double Der, Der1, Der2;
} FuncRes;

FuncRes Results[MAX_ITER];

double Func1(double x, double y, double y1) {
    return -1 * y + 2 * x - M_PI;
}

double Func2(double x, double y, double y1) {
    return -1 * y;
}

double Sol(double x) {
    return 2 * x - M_PI + M_PI * cos(x) + sin(x);
}

int Superposition(double a, double b, double y0, double yd0, double y10, double y1d0, double y20, double y2d0, double h, double eps) {
    int N = 0;
    bool solved1 = false, solved2 = false, solved3 = false; // Флаги завершения для каждого уравнения
    for (int i = 0; i < MAX_ITER; i++) {
        Results[i].Point = Results[i].Func = Results[i].Der = Results[i].Func1 = Results[i].Der1 = Results[i].Func2 = Results[i].Der2 = 0;
    }
    while (a <= b) {
        double prev = (N - 1 >= 0 ? Results[N - 1].Func : y0);
        double prev_d = (N - 1 >= 0 ? Results[N - 1].Der : yd0);
        double curr = prev, curr_d = prev_d;  
        if (!solved1) {
            EulerStep(&curr, &curr_d, h, a, Func1);
            if (eps != INFINITY && fabs(curr - prev) * 100000 < eps) {
                solved1 = true;
            }
        }
        Results[N].Func = curr;
        Results[N].Der = curr_d;
        double prev1 = (N - 1 >= 0 ? Results[N - 1].Func1 : y10);
        double prev_d1 = (N - 1 >= 0 ? Results[N - 1].Der1 : y1d0);
        double curr1 = prev1, curr_d1 = prev_d1;
        if (!solved2) {
            EulerStep(&curr1, &curr_d1, h, a, Func2);
            
            if (eps != INFINITY && fabs(curr1 - prev1) * 100000 < eps) {
                solved2 = true;
            }
        }
        Results[N].Func1 = curr1;
        Results[N].Der1 = curr_d1;
        double prev2 = (N - 1 >= 0 ? Results[N - 1].Func2 : y20);
        double prev_d2 = (N - 1 >= 0 ? Results[N - 1].Der2 : y2d0);
        double curr2 = prev2, curr_d2 = prev_d2;
        if (!solved3) {
            EulerStep(&curr2, &curr_d2, h, a, Func2);
            
            if (eps != INFINITY && fabs(curr2 - prev2) * 100000 < eps) {
                solved3 = true;
            }
        }
        
        Results[N].Func2 = curr2;
        Results[N].Der2 = curr_d2;
        Results[N].Point = a;
        N++;
        if (solved1 && solved2 && solved3) {
            break;
        }
        
        a += h;
    }
    return N;
}


void research() {
    double tmp1, tmp2;
    /*double h = 0.01;
    FILE* func_h1 = fopen("h1_x_y_err.txt", "w");
    FILE* func_h2 = fopen("h2_x_y_err.txt", "w");
    int N1 = Superposition(0, M_PI, 0, 1, 0, 0, 0, 0, h, INFINITY);
    tmp1 = (1 - Results[N1 - 1].Func - Results[N1 - 1].Func2) / (Results[N1 - 1].Func1 - Results[N1 - 1].Func2);
    tmp2 = (Results[N1 - 1].Func + Results[N1 - 1].Func1 - 1) / (Results[N1 - 1].Func1 - Results[N1 - 1].Func2);
    for (int i = 0; i < N1; i++) {
        double arg1, arg2, arg3;
        arg1 = Results[i].Point;
        arg2 = Results[i].Func + tmp1 * Results[i].Func1 + tmp2 * Results[i].Func2;
        arg3 = fabs(Results[i].Func + tmp1 * Results[i].Func1 + tmp2 * Results[i].Func2 - Sol(Results[i].Point));
        fprintf(func_h1, "%.15lf %.15lf %.15lf\n", arg1, arg2, arg3);
    }
    h = 0.0001;
    int N2 = Superposition(0, M_PI, 1, 1, 1, 1, 1, 1, h, INFINITY);
    tmp1 = (1 - Results[N2 - 1].Func - Results[N2 - 1].Func2) / (Results[N2 - 1].Func1 - Results[N2 - 1].Func2);
    tmp2 = (Results[N2 - 1].Func + Results[N2 - 1].Func1 - 1) / (Results[N2 - 1].Func1 - Results[N2 - 1].Func2);
    for (int i = 0; i < N2; i++) {
        double arg1, arg2, arg3;
        arg1 = Results[i].Point;
        arg2 = Results[i].Func + tmp1 * Results[i].Func1 + tmp2 * Results[i].Func2;
        arg3 = fabs(Results[i].Func + tmp1 * Results[i].Func1 + tmp2 * Results[i].Func2 - Sol(Results[i].Point));
        fprintf(func_h2, "%.15lf %.15lf %.15lf\n", arg1, arg2, arg3);
    }*/
    double max, aux;
    int N;
    FILE* RES_err_eps = fopen("error_eps.txt", "w");
    for (double eps = 1e-5; eps <= 0.1; eps *= 10) {
        max = 0;
        N = Superposition(0, M_PI, 0, 1, 1, 1, 1, 1, 0.1, eps);
        printf("N при eps = %.15lf: %d\n", eps, N);
        tmp1 = (1 - Results[N - 1].Func - Results[N - 1].Func2) / (Results[N - 1].Func1 - Results[N - 1].Func2);
        tmp2 = (Results[N - 1].Func + Results[N - 1].Func1 - 1) / (Results[N - 1].Func1 - Results[N - 1].Func2);
        for (int i = 0; i < N; i++) {
            aux = fabs(Results[i].Func + tmp1 * Results[i].Func1 + tmp2 * Results[i].Func2 - Sol(Results[i].Point));
            if (aux > max) {
                max = aux;
            }
        }
        fprintf(RES_err_eps, "%.15lf\n", max);
    }
}

int main() {
    research();
    return 0;
}
