#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define MAX_ITER 2000

double df(double x, double y) {
    return pow(x, 2) + 2 * x + y / (x + 2);
}

double f(double x) {
    return 0.5 * x + 1 + pow(x, 3) / 2 + pow(x, 2);
}

double RK_Without_H1(double (*f)(double, double), double a, double y0, double y, double h) {
    double h1, h2, h3, h4;
    h1 = y;
    h2 = f(a + h / 2.0, y0 + h / 2.0 * h1);
    h3 = f(a + h / 2.0, y0 + h / 2.0 * h2);
    h4 = f(a + h, y0 + h * h3);
    return y0 + h / 6.0 * (h1 + 2 * h2 + 2 * h3 + h4);
}

double** RK(double (*f)(double, double), double a, double b, double Y0, double eps) {
    double x0 = a, y0 = Y0;
    double h = (b - a) / 1;
    double* x = (double*)calloc(1, sizeof(double)), * y = (double*)calloc(1, sizeof(double)), * segment = (double*)calloc(1, sizeof(double));
    int iter = 0;
    x[iter] = x0, y[iter] = y0;
    do {
        double tmp = f(x0, y0);
        double y1 = RK_Without_H1(df, x0, y0, tmp, h);
        double y2 = RK_Without_H1(df, x0, y0, tmp, h / 2);
        tmp = f(x0 + h / 2, y2);
        double y3 = RK_Without_H1(df, x0 + h / 2, y2, tmp, h / 2);

        while (fabs(y3 - y1) / 15.0 > eps) {
            h /= 2;
            y1 = y2;
            double tmp1 = f(x0, y0);
            y2 = RK_Without_H1(df, x0, y0, tmp1, h / 2);
            tmp1 = f(x0 + h / 2, y2);
            y3 = RK_Without_H1(df, x0 + h / 2, y2, tmp1, h / 2);
        }
        iter++;
        x0 += h;
        y0 = y1;
        x = (double*)realloc(x, (iter + 1) * sizeof(double)), y = (double*)realloc(y, (iter + 1) * sizeof(double)), segment = (double*)realloc(segment, (iter + 1) * sizeof(double));
        x[iter] = x0;
        y[iter] = y0;
        segment[iter] = h;
    } while (x0 <= b && iter < MAX_ITER);
    double** res = (double**)malloc(4 * sizeof(double*));
    res[0] = x;
    res[1] = y;
    res[2] = segment;
    double* iter_ptr = (double*)malloc(sizeof(double));
    *iter_ptr = iter;
    res[3] = iter_ptr;
    return res;
}

double* euler_method(double (*f)(double, double), double y0, double a, double b, int n) {
    double h = (b - a) / (n - 1);
    double* y = (double*)malloc(n * sizeof(double));
    y[0] = y0;
    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = y[i] + h * f(a + i * h, y[i]);
    }
    return y;
}

double** predictor_corrector(double (*f)(double, double), double y0, double a, double b, int n) {
    double h = (b - a) / (n - 1);
    double* y1 = (double*)malloc(n * sizeof(double));
    double* y2 = (double*)malloc(n * sizeof(double));
    y1[0] = y0;
    double x = a;
    double* tmp = euler_method(f, y0, a, a + 2 * h, 2);
    double* x_arr = (double*)malloc(n * sizeof(double));
    x_arr[0] = a; x_arr[1] = a + h;
    y1[0] = tmp[0], y1[1] = tmp[1];
    y2[0] = tmp[0], y2[1] = tmp[1];
    for (int i = 2; i < n; i++) {
        x_arr[i] = a + i * h;
        y1[i] = y1[i - 1] + h / 2 * (3 * f(x_arr[i - 1], y1[i - 1]) - f(x_arr[i - 2], y1[i - 2]));
        y2[i] = y1[i - 1] + h / 2 * (f(x_arr[i], y1[i]) + f(x_arr[i - 1], y1[i - 1]));
    }

    double** res = (double**)malloc(2 * sizeof(double*));
    res[0] = x_arr;
    res[1] = y2;
    return res;
}


double find_max_error(double* x, int size) {
    double max = x[0];
    for (int i = 1; i < size; i++) {
        if (x[i] > max) {
            max = x[i];
        }
    }
    return max;
}

double find_error(double (*f)(double), double** tmp1, int size) {
    double* error = (double*)calloc(size, sizeof(double));
    for (int i = 0; i < size; i++) {
        error[i] = fabs(tmp1[1][i] - f(tmp1[0][i]));
    }
    return find_max_error(error, size);
}
void write_to_file2(char const file_name[], double* x1, double* x2, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.15f %.15f\n", x1[i], x2[i]);
    }
    fclose(file);
}

void write_to_file3(char const file_name[], double* x1, double* x2, double* x3, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.15f %.15f %.15f\n", x1[i], x2[i], x3[i]);
    }
    fclose(file);
}

void write_to_file4(char const file_name[], double* x1, double* x2, double* x3, double* x4, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.15f %.15f %.15f %.15f\n", x1[i], x2[i], x3[i], x4[i]);
    }
    fclose(file);
}

int main() {
    setlocale(LC_ALL, "Rus");

    int count_iter = 100;

    double* point1 = (double*)calloc(count_iter, sizeof(double));
    double* max_error1 = (double*)calloc(count_iter, sizeof(double));
    double* point2 = (double*)calloc(count_iter, sizeof(double));
    double* max_error2 = (double*)calloc(count_iter, sizeof(double));
    FILE* file = fopen("data.txt", "w");
    for (int i = 2; i < 100; i++) {
        double a = -1, b = 0;
        double h = (b - a) / i;
        double* x1 = (double*)calloc(i, sizeof(double));
        double* y1 = (double*)calloc(i, sizeof(double));
        double* error = (double*)calloc(i, sizeof(double));
        double y0 = 1;
        x1[0] = a, y1[0] = y0, error[0] = 0;
        int k = 0;
        double err1 = 0;
        for (double j = a + h; j < b; j += h) {
            k += 1;
            double numerical_solution = RK_Without_H1(df, j, y0, df(j, y0), h);
            x1[k] = j;
            y1[k] = numerical_solution;
            error[k] = fabs(y1[k] - f(j));
            y0 = numerical_solution;
            if (err1 < error[k]) {
                err1 = error[k];
                max_error1[i - 2] = err1;
                point1[i - 2] = j;
            }
        }
        y0 = 1, a = -1, b = 0;
        double** res = predictor_corrector(df, y0, a, b, i);
        double* x2 = res[0];
        double* y2 = res[1];
        double err2 = 0;
        for (int j = 0; j < i; j++) {
            if (err2 < fabs(y2[j] - f(x2[j]))) {
                err2 = fabs(y2[j] - f(x2[j]));
                max_error2[i - 2] = err2;
                point2[i - 2] = x2[j];
            }
        }
        printf("%d\n", i);
        fprintf(file, "%lf %lf %lf %lf %lf\n", h, point1[i-2], max_error1[i-2], point2[i-2], max_error2[i-2]);

    }
    fclose(file);
    return 0;
}
