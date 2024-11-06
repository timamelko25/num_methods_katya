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

#define PI 3.14159265358979323851


typedef struct {
	double y1;
	double y2;
}results;

double f(double x) {
	return 1 - sin(x);
}

double p(double x) {
	return cos(x);
}
double q(double x) {
	return sin(x);
}

double solution(double x) {
    return sin(x);
}

double F(double x, results y) {
	return f(x) - p(x) * y.y2 - q(x) * y.y1;
}

double F_0(double x, results y) {
	return -p(x) * y.y2 - q(x) * y.y1;
}


results F_vect(double x, results y) {
	results res = { y.y2, F(x,y) };
	return res;
}
results F_vect_0(double x, results y) {
	results res = { y.y2,F_0(x,y) };
	return res;
}// y(0) = A
// y(b) + y`(b) = B
results* find_initials(double A, double alpha_0, double alpha_1) {
	results u = { alpha_0*A/(pow(alpha_0,2)+ pow(alpha_1,2)),alpha_1*A / (pow(alpha_0,2) + pow(alpha_1,2)) };
	results v = { alpha_1,-alpha_0 };

	results* res= malloc(sizeof(results)*2);
	res[0] = u;
	res[1] = v;
	return res;
}


results* euler_method(results y0, double a, double b, results(*F)(double, results), int n) {
	double h = (b - a) / (n - 1);
	results* y = malloc(sizeof(results) * n);
	y[0] = y0;
	for (int i = 0;i < n - 1;i++) {
		results res = F(a + i * h, y[i]);
		y[i + 1].y1 = y[i].y1 + h * res.y1;
		y[i + 1].y2 = y[i].y2 + h * res.y2;
	}
	return y;
}

double* solve(results* initials,double betta_0,double betta_1, double B,double a,double b,results(*F)(double, results), results(*F_0)(double, results),int n) {
	results u0 = initials[0];
	results v0 = initials[1];

	results* u = euler_method(u0, a, b, F, n);
	results* v = euler_method(v0, a, b, F_0, n);
	double* y = malloc(sizeof(double) * n);
	double C = (B - betta_0 * u[n - 1].y1 - betta_1 * u[n - 1].y2) / (betta_0 * v[n - 1].y1 + betta_1 * v[n - 1].y2);

	for (int i = 0;i < n;i++) {
		y[i] = u[i].y1 + C * v[i].y1;
	}
	return y;
}


double** solve_runge(results* initials, double betta_0, double betta_1, double B, double a, double b, results(*F)(double, results), results(*F_0)(double, results), double eps, int* n_, int* len_h){
	double* x = malloc(sizeof(double));
	double* h_all = malloc(sizeof(long double));
	double* h_made = malloc(sizeof(long double));

	long double h = (b - a) / 2;
	x[0] = a;
	h_all[0] = h;

	results res;
	double c1,c2;

	results *u1 = malloc(sizeof(results)*2), * v1 = malloc(sizeof(results) * 2), u2, v2;
	u1[0] = initials[0];
	v1[0] = initials[1];
	u2 = initials[0];
	v2 = initials[1];

	int i = 1, k = 1, n = 10000;

	results* u = malloc(sizeof(results));
	results* v = malloc(sizeof(results));
	u[0] = initials[0];
	v[0] = initials[1];
	eps = pow(eps, 3.0/2);
	double difu,difv,dif;
	while (x[i - 1] + h <= b) {
		h_all = realloc(h_all, sizeof(long double) * (k + 1));
		h_all[k] = h;

		u1[1] = euler_method(u1[0], x[i - 1], x[i - 1] + h, F, 3)[2];
		v1[1] = euler_method(v1[0], x[i - 1], x[i - 1] + h, F_0, 3)[2];

		u2 = euler_method(u1[0], x[i - 1], x[i - 1] + h, F, 5)[4];
		v2 = euler_method(v1[0], x[i - 1], x[i - 1] + h, F_0, 5)[4];

		difu = fabs(u1[1].y1 - u2.y1);
		difv = fabs(v1[1].y1 - v2.y1);
	
		if (difu > eps/100 || difv > eps/100) {
			h /= 2;
		}
		else if (difu < eps / 500 && difv < eps / 500) {
			x = realloc(x, sizeof(double) * (i + 1));
			x[i] = x[i - 1] + h;
			u = realloc(u, sizeof(results) * (i + 1));
			u[i] = u2;
			v = realloc(v, sizeof(results) * (i + 1));
			v[i] = v2;
			h_made = realloc(h_made, sizeof(long double) * (i + 1));
			h_made[i - 1] = h;
			i += 1;
			u1[0] = u2;
			v1[0] = v2;
			if (x[i] + h * 2 <= b) {
				h *= 2;
			}
		}
		else {
			x = realloc(x, sizeof(double) * (i + 1));
			x[i] = x[i - 1] + h;
			u = realloc(u, sizeof(results) * (i + 1));
			u[i] = u2;
			v = realloc(v, sizeof(results) * (i + 1));
			v[i] = v2;
			h_made = realloc(h_made, sizeof(long double) * (i + 1));
			h_made[i-1] = h;
			i += 1;
			u1[0] = u2;
			v1[0] = v2;
		}
		k++;
	}
	h_all[k - 1] = h;
	h_made[i - 1] = h;
	*len_h = k;
	*n_ = i;

	double C = (B - betta_0 * u[i - 1].y1 - betta_1 * u[i - 1].y2) / (betta_0 * v[i - 1].y1 + betta_1 * v[i - 1].y2);
	double* y = malloc(sizeof(double) * i);
	for (int l = 0;l < i;l++) {
		y[l] = u[l].y1 + C * v[l].y1;
	}


	double** res_all = malloc(sizeof(double*) * 4);
	res_all[0] = x;
	res_all[1] = y;
	res_all[2] = h_all;
	res_all[3] = h_made;

	return res_all;
}

double* find_err(double* res, int n, double a, double h){
    double *err = (double*)calloc(n, sizeof(double));
    for (int i = 0;  i < n; i++){
        err[i] = fabs(res[i] - solution(a+i*h));
    }
    return err;
}

double find_err2(double **res, int n){
	double err = 0; 
	double max = 0;
	for (int i = 0; i < n; i++){
		if (max < fabs(solution(res[0][i]) - res[1][i])){
			max = fabs(solution(res[0][i]) - res[1][i]);
		}
	}
	err = max;
	return err;
}

void write_to_file(char const *filename, double *x1, double *x2, double *x3, int n){
    FILE* file = fopen(filename, "w");
    for (int i = 0; i < n; i++){
        fprintf(file, "%lf %lf %lf\n", x1[i], x2[i], x3[i]);
    }
    fclose(file);
}

void write_to_file2(char const *filename, double *x1, int n){
    FILE* file = fopen(filename, "w");
    for (int i = 0; i < n; i++){
        fprintf(file, "%lf\n", x1[i]);
    }
    fclose(file);
}


int main() {
	double A = 0, B = 1, alpha_0 = 1, alpha_1 = 0, betta_0 = 1, betta_1 = 1, a =0, b=PI;
	
	results* initials = find_initials(A, alpha_0, alpha_1);

	int n = 3;
	double h = (b - a) / (n - 1);
	double* y = solve(initials, betta_0, betta_1, B, a, b, F_vect, F_vect_0, n);

    double *error1 = find_err(y, n, a, h);
    double *x = malloc(sizeof(double) * n);
    for (int i = 0;i < n;i++) {
        x[i] = a + i * h;
    }

    write_to_file("h1_x_y_err.txt", x, y, error1, n);


	n = 100;
	h = (b - a) / (n - 1);
	y = solve(initials, betta_0, betta_1, B, a, b, F_vect, F_vect_0, n);

    double *error2 = find_err(y, n, a, h);
    double *x2 = malloc(sizeof(double) * n);
    for (int i = 0;i < n;i++) {
        x2[i] = a + i * h;
    }  
    write_to_file("h2_x_y_err.txt", x2, y, error2, n);

	double eps = 1;
	int n_=0, len_h=0;
	double** res;
    char path[100];
	double *err3 = (double*)calloc(6, sizeof(double));
	for (int i = 0; i < 6;i++) {
		eps /= 10;
		
		FILE* file, * file_h;
		sprintf(path, "runge_%i.txt", i);
		file = fopen(path, "w");
		sprintf(path, "runge_h_%i.txt", i);
		file_h = fopen(path, "w");

		res = solve_runge(initials, betta_0, betta_1, B, a, b, F_vect, F_vect_0, eps, &n_, &len_h);

		err3[i] = find_err2(res, n_);

		for (int i = 0;i < n_;i++) {
			fprintf(file, "%.16lf %.16lf %.20lf %.20lf\n", res[0][i], res[1][i], res[3][i], fabs(solution(res[0][i]) - res[1][i]));
		}
		for (int k = 0;k < len_h;k++) {
			fprintf(file_h, "%.20lf\n", res[2][k]);
		}
		free(res[0]);
		free(res[1]);
		free(res[2]);
		free(res[3]);
		free(res);
		fclose(file);
		fclose(file_h);
	}

	write_to_file2("error_eps.txt", err3, 6);


	return 0;
}