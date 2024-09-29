#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// 1 4 
// x^5-2.2x^3+0.5x^2-7x [-2;2]


#define result -10.8

double func(double x) {
    return pow(x, 5) - 2.2 * pow(x, 3) + 0.5 * pow(x, 2) - 7 * x;
}

double integrate_rec_l(double(*f)(double),double prev_S, int n, double a, double b) {
	double S;
	double h = (double)(b - a) / n;

	S = 0;
	for (int i = 2;i <= n;i += 2) {
		S += f(a + (i-1)*h);
	}

	return h * S + (0.5) * prev_S;
}

double* runge_write(double(*f)(double), double a, double b, double eps) {
	double S, S2 = 0;
	int n = 1;
	int counter = 0;
	S2 = (b - a) * f(a);

	do {

		counter++;

		S = S2;
		n *= 2;

		S2 = integrate_rec_l(f, S, n, a, b);

	} while (fabs(S2 - S) > eps);
    double *res = (double*)calloc(3, sizeof(double));
    res[0] = S2;
    res[1] = (b - a) / n;
    res[2] = counter;
    return res;
}

double* find_error(double* x, int size, double res) {
    double* error = (double*)calloc(size, sizeof(double));
    for (int i = 0; i < size; i++) {
        error[i] = fabs(res - x[i]);
    }
    return error;
}

void write_to_file3(char const file_name[], int* x1, double* x2, double* x3, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%d %.12f %.12f\n", x1[i], x2[i], x3[i]);
    }
    fclose(file);
}

int main() {
    setlocale(LC_ALL, "Rus");
    double eps = 0.1;
    double a = 0, b = 2;
    double count_iterations = 7;
    int* iterations = (int*)calloc(count_iterations, sizeof(int));
    double* results = (double*)calloc(count_iterations, sizeof(double));
    double* segments = (double*)calloc(count_iterations, sizeof(double));

    for (int i = 0; i < count_iterations; i++) {
        double *res = runge_write(func, a, b, eps);
        results[i] = res[0];
        segments[i] = res[1];
        iterations[i] = res[2];
        eps *= 0.1;
    }

    double* error = find_error(results, count_iterations, result);

    write_to_file3("iterations-error-segments.txt", iterations, error, segments, count_iterations);
    return 0;
}