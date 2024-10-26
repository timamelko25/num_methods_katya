#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

// 6 1
//  (2x+1)y’=4x+2y, (0, 4), y(a) = 1 Точное решение: y=(2x+1)ln(|2x+1|)+1

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

    int count_iter = 7;
    double a = -1, b = 0, y0 = 1;
    double h1 = 0.1;
    double h2 = 0.01;

    double* x = (double*)calloc((b - a) / h1 + 1, sizeof(double));
    double* y = (double*)calloc((b - a) / h1 + 1, sizeof(double));
    double* err = (double*)calloc((b - a) / h1 + 1, sizeof(double));
    double* sol = (double*)calloc((b - a) / h1 + 1, sizeof(double));
    int k = 0;
    x[0] = a, y[0] = y0, err[0] = 0, sol[0] = f(a);
    for (double i = a+h1; i < b; i += h1) {
        k += 1;
        double exact_solution = f(i);
        double tmp = df(i, y0);
        double numerical_solution = RK_Without_H1(df, i, y0, tmp, h1);
        double error = fabs(exact_solution - numerical_solution);
        x[k] = i;
        y[k] = numerical_solution;
        err[k] = error;
        sol[k] = exact_solution;
        y0 = numerical_solution;
    }

    write_to_file4("h1_x_y_err.txt", x, sol, y, err, k);

    y0 = 1;
    k = 0;
    double* x2 = (double*)calloc((b - a) / h2 + 1, sizeof(double));
    double* y2 = (double*)calloc((b - a) / h2 + 1, sizeof(double));
    double* err2 = (double*)calloc((b - a) / h2 + 1, sizeof(double));
    double* sol2 = (double*)calloc((b - a) / h2 + 1, sizeof(double));
    x2[0] = a, y2[0] = y0, err2[0] = 0, sol2[0] = f(a);
    for (double i = a+h2; i < b; i += h2) {
        k += 1;
        double exact_solution = f(i);
        double tmp = df(i, y0);
        double numerical_solution = RK_Without_H1(df, i, y0, tmp, h2);
        double error = fabs(exact_solution - numerical_solution);
        x2[k] = i;
        y2[k] = numerical_solution;
        err2[k] = error;
        sol2[k] = exact_solution;
        y0 = numerical_solution;
    }

    write_to_file4("h2_x_y_err.txt", x2, sol2, y2, err2, k);

    double* error = (double*)calloc(count_iter, sizeof(double));
    double* segment = (double*)calloc(count_iter, sizeof(double));
    double eps = 0.1;
    y0 = 1;

    for (int i = 0; i < count_iter; i++) {
        double** res3 = RK(df, a, b, y0, eps / 100);
        error[i] = find_error(f, res3, (int)(*res3[3]));
        segment[i] = find_max_error(res3[2], (int)(*res3[3]));
        eps *= 0.1;
    }

    write_to_file2("segment-error.txt", segment, error, count_iter);

    free(error);
    free(segment);

    return 0;
}
