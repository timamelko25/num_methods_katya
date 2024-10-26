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
    double prev, prev_d, curr, curr_d, tmp0, der, tmp1, prev_1;
    double prev1, prev_d1, curr1, curr_d1, tmp01, der1, tmp11, prev_11;
    double prev2, prev_d2, curr2, curr_d2, tmp02, der2, tmp12, prev_12;
    int N = 0, i;
    for (i = 0; i < MAX_ITER; i++) {
        Results[i].Point = Results[i].Func = Results[i].Der = Results[i].Func1 = Results[i].Der1 = Results[i].Func2 = Results[i].Der2 = 0;
    }
    while (a <= b) {
        prev_d = (N - 1 >= 0 ? Results[N - 1].Der : yd0);
        prev = (N - 1 >= 0 ? Results[N - 1].Func : y0);
        tmp0 = prev + h * prev_d;
        der = prev_d + h * Func1(a, prev, prev_d);
        prev_1 = tmp0;

        prev_d1 = (N - 1 >= 0 ? Results[N - 1].Der1 : y1d0);
        prev1 = (N - 1 >= 0 ? Results[N - 1].Func1 : y10);
        tmp01 = prev1 + h * prev_d1;
        der1 = prev_d1 + h * Func2(a, prev1, prev_d1);
        prev_11 = tmp01;

        prev_d2 = (N - 1 >= 0 ? Results[N - 1].Der2 : y2d0);
        prev2 = (N - 1 >= 0 ? Results[N - 1].Func2 : y20);
        tmp02 = prev2 + h * prev_d2;
        der2 = prev_d2 + h * Func2(a, prev2, prev_d2);
        prev_12 = tmp02;

        if (eps != INFINITY) {
            while (N < MAX_ITER) {
                tmp0 = prev + h / 2 * prev_d;
                tmp1 = prev_d + h / 2 * Func1(a, prev, prev_d);
                curr = tmp0 + h / 2 * tmp1;
                curr_d = tmp1 + h / 2 * Func1(a + h / 2, tmp0, tmp1);

                tmp01 = prev1 + h / 2 * prev_d1;
                tmp11 = prev_d1 + h / 2 * Func2(a, prev1, prev_d1);
                curr1 = tmp01 + h / 2 * tmp11;
                curr_d1 = tmp11 + h / 2 * Func2(a + h / 2, tmp01, tmp11);

                tmp02 = prev2 + h / 2 * prev_d2;
                tmp12 = prev_d2 + h / 2 * Func2(a, prev2, prev_d2);
                curr2 = tmp02 + h / 2 * tmp12;
                curr_d2 = tmp12 + h / 2 * Func2(a + h / 2, tmp02, tmp12);

                if (fabs(curr - prev_1) * 100000 < eps && fabs(curr1 - prev_11) * 100000 < eps && fabs(curr2 - prev_12) * 100000 < eps) {
                    Results[N].Func = curr;
                    Results[N].Der = curr_d;
                    Results[N].Func1 = curr1;
                    Results[N].Der1 = curr_d1;
                    Results[N].Func2 = curr2;
                    Results[N].Der2 = curr_d2;
                    Results[N].Point = a;
                    //fprintf(F_a, "%.15lf\n", Results[N].Point);
                    //fprintf(F_h, "%.15lf\n", h);
                    N++;
                    break;
                }
                h /= 2;
                prev_1 = tmp0;
                prev_11 = tmp01;
                prev_12 = tmp02;
            }
            y0 = curr - (curr - prev_1), yd0 = curr_d, a += h;
            y10 = curr1 - (curr1 - prev_11), y1d0 = curr_d1;
            y20 = curr2 - (curr2 - prev_12), y2d0 = curr_d2;
        }
        else {
            Results[N].Func = tmp0;
            Results[N].Der = der;
            Results[N].Func1 = tmp01;
            Results[N].Der1 = der1;
            Results[N].Func2 = tmp02;
            Results[N].Der2 = der2;
            Results[N].Point = a;
            N++;
            a += h;
        }
    }
    //fclose(F_a);
    //fclose(F_h);
    return N;
}

void research() {
    double tmp1, tmp2;
    double h = 0.01;
    FILE* func_h1 = fopen("h1_x_y_err.txt", "w");
    FILE* func_h2 = fopen("h2_x_y_err.txt", "w");
    int N1 = Superposition(0, M_PI, 0.0, 0, 0, 0.0, 0.0, 0, h, INFINITY);
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
    int N2 = Superposition(0, M_PI, 0.0, 0, 0, 0.0, 0.0, 0, h, INFINITY);
    tmp1 = (1 - Results[N2 - 1].Func - Results[N2 - 1].Func2) / (Results[N2 - 1].Func1 - Results[N2 - 1].Func2);
    tmp2 = (Results[N2 - 1].Func + Results[N2 - 1].Func1 - 1) / (Results[N2 - 1].Func1 - Results[N2 - 1].Func2);
    for (int i = 0; i < N2; i++) {
        double arg1, arg2, arg3;
        arg1 = Results[i].Point;
        arg2 = Results[i].Func + tmp1 * Results[i].Func1 + tmp2 * Results[i].Func2;
        arg3 = fabs(Results[i].Func + tmp1 * Results[i].Func1 + tmp2 * Results[i].Func2 - Sol(Results[i].Point));
        fprintf(func_h2, "%.15lf %.15lf %.15lf\n", arg1, arg2, arg3);
    }
    double max, aux;
    int N;
    FILE* RES_err_eps = fopen("error_eps.txt", "w");
    for (double eps = 1e-5; eps <= 0.1; eps *= 10) {
        max = 0;
        N = Superposition(0, M_PI, 0, 0, 0, 0, 0, 0, 0.1, eps);
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
