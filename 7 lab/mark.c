#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define pi M_PI
#define MAX_ITER 200000

typedef struct FuncRes {
    double Point;
    double Func, Func1, Func2;
    double Der, Der1, Der2;
} FuncRes;

FuncRes Results[MAX_ITER];

double Sol(double x) {
    return 2 * x - M_PI + M_PI * cos(x) + sin(x);
}

double Func1(double x, double y, double y1, int k) {
    double fluct = (k / 100.0) * ((double)(rand()) / RAND_MAX * 2 - 1.0);
    return (1 + fluct) * (-1 * y + 2 * x - M_PI);
}

double Func2(double x, double y, double y1, int k) {
    double fluct = (k / 100.0) * ((double)(rand()) / RAND_MAX * 2 - 1.0);
    return (1 + fluct) * (-1 * y);
}

int Euler(double a, double b, double y0, double yd0, double y10, double y1d0, double y20, double y2d0, double h, double eps, int k) {
    double prev, prev_d, curr, curr_d, tmp0, der, tmp1, prev_1;
    double prev1, prev_d1, curr1, curr_d1, tmp01, der1, tmp11, prev_11;
    double prev2, prev_d2, curr2, curr_d2, tmp02, der2, tmp12, prev_12;
    int N = 0, i;
    //FILE* F_h = fopen("h.txt", "w");
    //FILE* F_a = fopen("a.txt", "w");
    for (i = 0; i < MAX_ITER; i++) {
        Results[i].Point = Results[i].Func = Results[i].Der = Results[i].Func1 = Results[i].Der1 = Results[i].Func2 = Results[i].Der2 = 0;
    }
    while (a <= b) {
        prev_d = (N - 1 >= 0 ? Results[N - 1].Der : yd0);
        prev = (N - 1 >= 0 ? Results[N - 1].Func : y0);
        tmp0 = prev + h * prev_d;
        der = prev_d + h * Func1(a, prev, prev_d, k);
        prev_1 = tmp0;

        prev_d1 = (N - 1 >= 0 ? Results[N - 1].Der1 : y1d0);
        prev1 = (N - 1 >= 0 ? Results[N - 1].Func1 : y10);
        tmp01 = prev1 + h * prev_d1;
        der1 = prev_d1 + h * Func2(a, prev1, prev_d1, k);
        prev_11 = tmp01;

        prev_d2 = (N - 1 >= 0 ? Results[N - 1].Der2 : y2d0);
        prev2 = (N - 1 >= 0 ? Results[N - 1].Func2 : y20);
        tmp02 = prev2 + h * prev_d2;
        der2 = prev_d2 + h * Func2(a, prev2, prev_d2, k);
        prev_12 = tmp02;

        if (eps != INFINITY) {
            while (1) {
                tmp0 = prev + h / 2 * prev_d;
                tmp1 = prev_d + h / 2 * Func1(a, prev, prev_d, k);
                curr = tmp0 + h / 2 * tmp1;
                curr_d = tmp1 + h / 2 * Func1(a + h / 2, tmp0, tmp1, k);

                tmp01 = prev1 + h / 2 * prev_d1;
                tmp11 = prev_d1 + h / 2 * Func2(a, prev1, prev_d1, k);
                curr1 = tmp01 + h / 2 * tmp11;
                curr_d1 = tmp11 + h / 2 * Func2(a + h / 2, tmp01, tmp11, k);

                tmp02 = prev2 + h / 2 * prev_d2;
                tmp12 = prev_d2 + h / 2 * Func2(a, prev2, prev_d2, k);
                curr2 = tmp02 + h / 2 * tmp12;
                curr_d2 = tmp12 + h / 2 * Func2(a + h / 2, tmp02, tmp12, k);

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

void Task1() {
    double C1;
    double C2;
    double h1 = 0.01;
    double h2 = 0.005;
    FILE* RES_x_h1 = fopen("draw_h1_x.txt", "w");
    FILE* RES_y_h1 = fopen("draw_h1_y.txt", "w");
    FILE* RES_x_h2 = fopen("draw_h2_x.txt", "w");
    FILE* RES_y_h2 = fopen("draw_h2_y.txt", "w");
    FILE* RES_err_h1_y = fopen("draw_err_h1_y.txt", "w");
    FILE* RES_err_h2_y = fopen("draw_err_h2_y.txt", "w");
    int N1 = Euler(0, pi, 0, 0, 0, 0, 0, 0, h1, INFINITY, 0);
    C1 = (1 - Results[N1 - 1].Func - Results[N1 - 1].Func2) / (Results[N1 - 1].Func1 - Results[N1 - 1].Func2);
    C2 = (Results[N1 - 1].Func + Results[N1 - 1].Func1 - 1) / (Results[N1 - 1].Func1 - Results[N1 - 1].Func2);
    for (int i = 0; i < N1; i++) {
        fprintf(RES_x_h1, "%.15lf\n", Results[i].Point);
        fprintf(RES_y_h1, "%.15lf\n", Results[i].Func + C1 * Results[i].Func1 + C2 * Results[i].Func2);
        fprintf(RES_err_h1_y, "%.15lf\n", fabs(Results[i].Func + C1 * Results[i].Func1 + C2 * Results[i].Func2 - Sol(Results[i].Point)));
    }
    int N2 = Euler(0, pi, 0, 0, 0, 0, 0, 0, h2, INFINITY, 0);
    C1 = (1 - Results[N2 - 1].Func - Results[N2 - 1].Func2) / (Results[N2 - 1].Func1 - Results[N2 - 1].Func2);
    C2 = (Results[N2 - 1].Func + Results[N2 - 1].Func1 - 1) / (Results[N2 - 1].Func1 - Results[N2 - 1].Func2);
    for (int i = 0; i < N2; i++) {
        fprintf(RES_x_h2, "%.15lf\n", Results[i].Point);
        fprintf(RES_y_h2, "%.15lf\n", Results[i].Func + C1 * Results[i].Func1 + C2 * Results[i].Func2);
        fprintf(RES_err_h2_y, "%.15lf\n", fabs(Results[i].Func + C1 * Results[i].Func1 + C2 * Results[i].Func2 - Sol(Results[i].Point)));
    }
}

void Task2() {
    double C1;
    double C2;
    double max, aux;
    int N;
    FILE* RES_err_eps = fopen("error_eps.txt", "w");
    for (double eps = 1e-5; eps <= 0.1; eps *= 10) {
        max = 0;
        N = Euler(0, pi, 0, 0, 0, 0, 0, 0, 0.1, eps, 0);
        C1 = (1 - Results[N - 1].Func - Results[N - 1].Func2) / (Results[N - 1].Func1 - Results[N - 1].Func2);
        C2 = (Results[N - 1].Func + Results[N - 1].Func1 - 1) / (Results[N - 1].Func1 - Results[N - 1].Func2);
        for (int i = 0; i < N; i++) {
            aux = fabs(Results[i].Func + C1 * Results[i].Func1 + C2 * Results[i].Func2 - Sol(Results[i].Point));
            if (aux > max) {
                max = aux;
            }
        }
        fprintf(RES_err_eps, "%.15lf\n", max);
    }
    N = Euler(0, pi, 0, 0, 1, 0, 0, 1, 0.5, 1e-3, 0);
}

void main(void) {
    Task1();
    Task2();
}