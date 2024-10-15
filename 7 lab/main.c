#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

// 7 1 1 ГУ 2
// 

double left_part(double x) {
    return 4 * x * x + 3; 
}

double right_part(double x) {
    return exp(-x * x); 
}

double func(double x) {
    return exp(-x * x); 
}

double* uniform_grid(double a, double b, int n) {
    double* x = (double*)malloc(n * sizeof(double));
    double h = (b - a) / (n - 1); 

    for (int i = 0; i < n; i++) {
        x[i] = a + i * h;
    }
    return x;
}

double* smooth_grid(double a, double b, int n, double factor) {
    double* x = (double*)malloc((n+1) * sizeof(double));
    int i;
    double h = ((b - a) * 0.55 - (b - a) * 0.45) / (n - (n * factor) * 2);
    double h1 = (b - a) * 0.45 / (n * factor);

    for (i = 0; i <= (n * factor); i++)
        x[i] = a + i * h1;
    for (i = (n * factor) + 1; i <= n - (n * factor); i++)
        x[i] = (b - a) * 0.45 + (i - (n * factor)) * h;
    for (i = n - (n * factor) + 1; i <= n; i++)
        x[i] = (b - a) * 0.55 + (i - n + (n * factor)) * h1;
    return x;
}

void thomas_solve(double* a, double* b, double* c, double* d, double* y, int N) {
    double* c_star = (double*)malloc(N * sizeof(double));
    double* d_star = (double*)malloc(N * sizeof(double));

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < N; i++) {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    y[N - 1] = d_star[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        y[i] = d_star[i] - c_star[i] * y[i + 1];
    }

    free(c_star);
    free(d_star);
}

void FDM(double* x, double* y, double a, double b, double alpha, double beta, int N) { 
    double h = (b - a) / (N - 1);

    double* A = (double*)malloc((N - 2) * sizeof(double));
    double* B = (double*)malloc((N - 2) * sizeof(double));
    double* C = (double*)malloc((N - 2) * sizeof(double));
    double* F = (double*)malloc((N - 2) * sizeof(double));

    for (int i = 0; i < N - 2; ++i) {
        double xi = x[i + 1];
        A[i] = 1 - (h * 4 * xi) / 2;
        B[i] = -2 + h * h * left_part(xi);
        C[i] = 1 + (h * 4 * xi) / 2;
        F[i] = h * h * right_part(xi);
    }

    F[0] -= (1 - (h * 4 * x[1]) / 2) * alpha;
    F[N - 3] -= (1 + (h * 4 * x[N - 2]) / 2) * beta;

    thomas_solve(A, B, C, F, y + 1, N - 2);

    y[0] = alpha;
    y[N - 1] = beta;

    free(A);
    free(B);
    free(C);
    free(F);
}

double** FDM_EPS(double a, double b, double alpha, double beta, double eps) {
    int n = 10;
    double *x, *y, *x2, *y2;
    double mid_x = (a + b) / 2.0;
    double max_error = eps + 1;

    while (max_error > eps) {
        x = uniform_grid(a, b, n);
        y = (double*)calloc(n, sizeof(double));
        FDM(x, y, a, b, alpha, beta, n);

        x2 = uniform_grid(a, b, 2 * n);
        y2 = (double*)calloc(2 * n, sizeof(double));
        FDM(x2, y2, a, b, alpha, beta, 2 * n);

        int mid_index_n = n / 2;
        int mid_index_2n = n;

        max_error = fabs(y[mid_index_n] - y2[mid_index_2n]);

        if (max_error <= eps) {
            break;
        }

        free(x);
        free(y);
        free(x2);
        free(y2);
        n *= 2;
    }
    printf("%d\n", n);
    double** res = (double**)malloc(3 * sizeof(double*));
    res[0] = x2;
    res[1] = y2;
    res[2] = (double*)(intptr_t)(2 * n);

    return res;
}
double* find_error(double* x, double* y, int size) {
    double* error = (double*)calloc(size, sizeof(double));
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        error[i] = fabs(func(x[i]) - y[i]);
    }
    return error;
}

double max(double* x, int size) {
    double max = x[0];
    for (int i = 1; i < size; i++) {
        if (x[i] > max) {
            max = x[i];
        }  
    }
    return max;
}

double find_error2(double (*f)(double), double** tmp1, int size) {
    double* error = (double*)calloc(size, sizeof(double));
    for (int i = 0; i < size; i++) {
        error[i] = fabs(tmp1[1][i] - f(tmp1[0][i]));
    }
    return max(error, size);
}

void write_to_file3(char const file_name[], double* x1, double* x2, double* x3, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.15f %.15f %.15f\n", x1[i], x2[i], x3[i]);
    }
    fclose(file);
}

int main() {
    double a = 0, b = 1;
    double A = func(a), B = func(b);

    int size_4 = 4, size_8 = 8;
    double* x_4 = uniform_grid(a, b, size_4);
    double* y_4 = (double*)malloc(size_4 * sizeof(double));

    double* x_8 = uniform_grid(a, b, size_8);
    double* y_8 = (double*)malloc(size_8 * sizeof(double));

    FDM(x_4, y_4, a, b, A, B, size_4);
    FDM(x_8, y_8, a, b, A, B, size_8);

    double* error1 = find_error(x_4, y_4, size_4);
    double* error2 = find_error(x_8, y_8, size_8);

    write_to_file3("x1-y1-error1.txt", x_4, y_4, error1, size_4);
    write_to_file3("x2-y2-error2.txt", x_8, y_8, error2, size_8);

    double** res = FDM_EPS(a, b, A, B, 0.0001);

    double* x = res[0];
    double* y = res[1];
    int n = (intptr_t)res[2];

    double* error = find_error(x, y, n);
    write_to_file3("x3-y3-error3.txt", x, y, error, n);

    double eps = 0.1;
    FILE* file = fopen("error3.txt", "w");
    for(int i = 0; i < 5; i++){
        double** res = FDM_EPS(a, b, A, B, eps);
        double* x = res[0];
        double* y = res[1];
        int n = (intptr_t)res[2];
        fprintf(file, "%.16f\n",find_error2(func, res, res[2]));
        eps /= 10;
    }


    free(x_4);
    free(y_4);
    free(x_8);
    free(y_8);
    free(error1);
    free(error2);


    return 0;
}
