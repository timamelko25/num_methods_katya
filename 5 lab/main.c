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
    return pow(x,2) + 2 * x  + y / (x + 2); 
}

double f(double x) {
    return 0.5 * x + 1 + pow(x,3) / 2 + pow(x,2);
}

double RK_Step(double (*f)(double, double),double x, double y0, double h) {
    double h1, h2, h3, h4;
    h1 = f(x, y0);
    h2 = f(x + h / 2.0, y0 + h / 2.0 * h1);
    h3 = f(x + h / 2.0, y0 + h / 2.0 * h2);
    h4 = f(x + h, y0 + h * h3);
    return y0 + (h * (h1 + 2 * h2 + 2 * h3 + h4)) / 6;
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
    double* x = (double*)calloc(1, sizeof(double)), *y = (double*)calloc(1, sizeof(double)), *segment = (double*)calloc(1, sizeof(double));
    int iter = 0;
    x[iter] = x0, y[iter] = y0;
    do{
        double y1 = RK_Step(df, x0, y0, h);
        double y2 = RK_Step(df, x0, y0, h/2);
        double y3 = RK_Step(df, x0+h/2, y2, h/2);
        while(fabs(y3 - y1) / 15.0 > eps) {
            h = h / 2;
            y1 = y2;
            y2 = RK_Step(df, x0, y0, h/2);
            y3 = RK_Step(df, x0+h/2, y2, h/2);
        }
        iter++;
        x0 += h;
        y0 = y1;
        x = (double*)realloc(x, (iter+1) * sizeof(double)), y = (double*)realloc(y, (iter+1) * sizeof(double)), segment = (double*)realloc(segment, (iter+1) * sizeof(double));
        x[iter] = x0;
        y[iter] = y0;
        segment[iter] = h;
    }while(x0 <= b && iter < MAX_ITER);
    double **res = (double**)malloc(4 * sizeof(double*));
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

int main() {
    setlocale(LC_ALL, "Rus");

    int count_iter = 7;
    double a = -1, b = 0, y0 = 1;
    double h1 = 0.1;
    double h2 = 0.01;

    FILE *F1, *F2;
    F1 = fopen("results_h1.txt", "w");
    F2 = fopen("results_h2.txt", "w");

    for (double i = a; i < b; i += h1) {
        double exact_solution = f(i);
        double numerical_solution = RK_Step(df, i, y0, h1);
        double error = fabs(exact_solution - numerical_solution); 

        fprintf(F1, "%.5lf %.5lf %.5lf %.5lf\n", i, exact_solution, numerical_solution, error);
        y0 = numerical_solution; 
    }

    y0 = 1;

    for (double i = a; i < b; i += h2) {
        double exact_solution = f(i); 
        double numerical_solution = RK_Step(df, i, y0, h2); 
        double error = fabs(exact_solution - numerical_solution); 

        fprintf(F2, "%.5lf %.5lf %.5lf %.5lf\n", i, exact_solution, numerical_solution, error);
        y0 = numerical_solution;
    }

    fclose(F1);
    fclose(F2);

    double* error = (double*)calloc(count_iter, sizeof(double));
    double* segment = (double*)calloc(count_iter, sizeof(double));
    double eps = 0.1;
    y0 = 1;

    for (int i = 0; i < count_iter; i++) {
        double** res3 = RK(df, a, b, y0, eps/100);
        error[i] = find_error(f, res3, (int)(*res3[3]));
        segment[i] = find_max_error(res3[2], (int)(*res3[3]));
        eps *= 0.1;
    }

    write_to_file2("segment-error.txt", segment, error, count_iter);

    free(error);
    free(segment);

    return 0;
}
