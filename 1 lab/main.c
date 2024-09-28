#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


// 23 variant 1 3
// y=lg(x-1) + 0.5x
// y=x^3+0.4x^2+0.6|x|-1.6

// x y5 y7 y10 

double f1 (double x) {
    return (2*x-log10(x));
}

double f1_derivative_6 (double x) {
    return 120 / (log(10)* pow(x,6));
}

double sgn(double x) {
    if (x > 0) {
        return 1;
    } else if (x < 0) {
        return -1;
    } else {
        return 0;
    }
}

double f2 (double x) {
    return sgn(x)*pow(x,4)-18*pow(x,2) + 6;
}

double factorial(int n) {
    double result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

double* chebish_mesh(double a, double b, int size) {
    double* tmp = (double*)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        tmp[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * M_PI / (2 * size));
    }
    return tmp;
}

double* simple_mesh(double a, double b, int size) {
    double* tmp = (double*)malloc(size*sizeof(double));
    double tmp1 = (b-a)/(size-1);
    for (int i = 0; i < size; i++) {
        tmp[i] = a + i * tmp1;
    }
    return tmp;
}

double lagranj_poly(double x, double* x_arr, double* y_arr, int size) {
    double result  = 0;
    for(int i = 0; i < size; i++) {
        double tmp1 = y_arr[i];
        for (int j = 0; j < size; j++) {
            if (i != j) {
                tmp1 *=(x - x_arr[j]) / (x_arr[i] - x_arr[j]);
            }
        }
        result += tmp1;
    }
    return result;
}

void calculate_func(double (*f)(double), double *x, double *y, int size){
    for(int i = 0; i < size; i++) {
        y[i] = f(x[i]);
    }
}

void write_to_file1(char const file_name[], double* x, double* y, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%lf %.25f\n", x[i], y[i]);
    }
    fclose(file);
}

void write_to_file2(char const file_name[], double* x, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.20f\n", x[i]);
    }
    fclose(file);
}

void write_to_file3(char const file_name[], int* x, double *y, int size) {
    FILE* file = fopen(file_name, "w");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.20f %.20f\n", x[i], y[i]);
    }
    fclose(file);
}

double* find_error(double* poly, double* y, int size) {
    double *tmp = (double*)calloc(size, sizeof(double));
    for(int i = 0; i < size; i++) {
        tmp[i] = fabs(poly[i] - y[i]);
    }
    return tmp;
}

double find_theor_error_func1(double x, double *x_arr, double M, int size) {
    double product = 1.0;
    for (int i = 0; i < size; i++) {
        product *= fabs(x - x_arr[i]);
    }

    double error = M * product / factorial(size + 1);
    return error;
}


double* find_differance(double *x, double *y, int size) {
    double *tmp = (double*)calloc(size, sizeof(double));
    for(int i = 0; i < size; i++) {
        tmp[i] = fabs(x[i] - y[i]);
    }
    return tmp;
}

double* find_differance_reverse(double *x, double *y, int size) {
    double *tmp = (double*)calloc(size, sizeof(double));
    for(int i = 0; i < size; i++) {
        tmp[i] = fabs(y[i] - x[i]);
    }
    return tmp;
}

double find_max(double *x, int size){
    double max = x[0];
    for(int i = 0; i < size; i++) {
        if(x[i] > max) {
            max = x[i];
        }
    }
    return max;
}

void research(double (*f)(double), double* (*mesh)(double,double,int), char const file_name[], double a, double b, int size) {
    double *x_mesh = mesh(a, b, size);
    double *y_mesh = (double*)calloc(size, sizeof(double));

    double *x = (double*)calloc(10000, sizeof(double));
    x = mesh(a, b, 10000);
    double *y = (double*)calloc(10000, sizeof(double));

    calculate_func(f, x_mesh, y_mesh, size);
    calculate_func(f, x, y, 10000);

    double *lagranj = (double*)calloc(10000, sizeof(double));

    for(int i = 0; i < 10000; i++) {
        lagranj[i] = lagranj_poly(x[i], x_mesh, y_mesh, size);
    }

    double *error = find_error(lagranj, y, 10000);

    char str[100];
    sprintf(str, "%s_%d.txt", file_name, size);
    write_to_file1(str, x_mesh, y_mesh, size);

    sprintf(str, "%s_%d_10000.txt", file_name, size);
    write_to_file1(str, x, lagranj, 10000);
    
    sprintf(str, "%s_error_%d.txt", file_name, size);
    write_to_file1(str, x, error, 10000);
}

void research_theor_error(double a, double b, int size) {
    double *x_mesh = simple_mesh(a, b, size);
    double *y_mesh = (double*)calloc(size, sizeof(double));

    double *x = (double*)calloc(10000, sizeof(double));
    x = simple_mesh(a, b, 10000);
    double *y = (double*)calloc(10000, sizeof(double));

    calculate_func(f1, x_mesh, y_mesh, size);
    calculate_func(f1, x, y, 10000);

    double *lagranj = (double*)calloc(10000, sizeof(double));

    for(int i = 0; i < 10000; i++) {
        lagranj[i] = lagranj_poly(x[i], x_mesh, y_mesh, size);
    }

    double *theor_error = (double*)calloc(10000, sizeof(double));
    double *max = (double*)calloc(10000, sizeof(double));
    for(int i = 0; i < 10000; i++) {
        max[i] = f1_derivative_6(x[i]);
    }
    double M = find_max(max, 10000);
    for (int i = 0; i < 10000; i++) {
        theor_error[i] = fabs(find_theor_error_func1(x[i], x_mesh, M, size));
    }
    write_to_file1("theor_error.txt", x, theor_error, 10000);
}

void research_nodes(double (*f)(double), double* (*mesh)(double,double,int), char const file_name[], double a, double b, int place1, int place2) {
    int *n_nodes = (int*)calloc(96, sizeof(int));
    for(int i = 0; i < 96; i++) {
        n_nodes[i] = 5+i;
    }

    double *lagranj = (double*)calloc(10000, sizeof(double));
    double *max_error = (double*)calloc(96, sizeof(double));

    double *point1 = (double*)calloc(96, sizeof(double));
    double *point2 = (double*)calloc(96, sizeof(double));
    

    for(int i = 0; i < 96; i++) {
        double *x_mesh = mesh(a, b, n_nodes[i]);
        double *y_mesh = (double*)calloc(n_nodes[i], sizeof(double));
        calculate_func(f, x_mesh, y_mesh, n_nodes[i]);

        double *x = mesh(a, b, 10000);
        double *y = (double*)calloc(10000, sizeof(double));
        calculate_func(f, x, y, 10000);
        for (int j = 0; j < 10000; j++) {
            lagranj[j] = lagranj_poly(x[j], x_mesh, y_mesh, n_nodes[i]);
    
        }
        double* tmp = find_differance(lagranj, y, 10000);
        max_error[i] = find_max(tmp, 10000);
        point1[i] = fabs(lagranj[place1]-y[place1]);
        point2[i] = fabs(lagranj[place2]-y[place2]);
    }
    char str[100];
    sprintf(str, "%s_max_error.txt", file_name);
    write_to_file2(str, max_error, 96);

    sprintf(str, "%s_point1.txt", file_name);
    write_to_file2(str, point1, 96);

    sprintf(str, "%s_point2.txt", file_name);
    write_to_file2(str, point2, 96);
}

int main() {
    setlocale(LC_ALL, "Rus");
    // func1
    double a1 = 1, b1 = 10;
    research(f1, simple_mesh, "func1_mesh", a1, b1, 5);
    research(f1, simple_mesh, "func1_mesh", a1, b1, 7);
    research(f1, simple_mesh, "func1_mesh", a1, b1, 10);
    research(f1, chebish_mesh, "func1_chebish", a1, b1, 5);
    research(f1, chebish_mesh, "func1_chebish", a1, b1, 7);
    research(f1, chebish_mesh, "func1_chebish", a1, b1, 10);

    // func 2
    double a2 = 0, b2 = 2;
    research(f2, simple_mesh, "func2_mesh", a2, b2, 5);
    research(f2, simple_mesh, "func2_mesh", a2, b2, 7);
    research(f2, simple_mesh, "func2_mesh", a2, b2, 10);
    research(f2, chebish_mesh, "func2_chebish", a2, b2, 5);
    research(f2, chebish_mesh, "func2_chebish", a2, b2, 7);
    research(f2, chebish_mesh, "func2_chebish", a2, b2, 10);

    research_theor_error(a1, b1, 5);

    research_nodes(f1, simple_mesh, "func1_mesh", a1, b1, 100, 5000);
    research_nodes(f1, chebish_mesh, "func1_chebish", a1, b1, 9366, 5001);
    research_nodes(f2, simple_mesh, "func2_mesh", a2, b2, 100, 5000);
    research_nodes(f2, chebish_mesh, "func2_chebish", a2, b2, 9366, 5001);
    return 0;
}