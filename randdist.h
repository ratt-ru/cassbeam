#ifndef __GAUSS__
#define __GAUSS__

#include <stdlib.h>
#include "vector.h"

#define rand_int(n) 	(1+(int)(((float)(n))/(RAND_MAX+1.0)))

#define hypot(x, y)	sqrt((x)*(x) + (y)*(y))

double rand_one();
double rand_pm_one();
double rand_gauss();
double rand_rice(double m, double sigma);
double rice_moment(double m, double sigma, int moment);
int rice_params2(const Vector x, const Vector y, double *mag, double *sigma);
int rice_params2a(const Vector x, const Vector y, double *mag, double *sigma);
int rice_params(const Vector data, double *mag, double *sigma);
int rice_params3(const Vector data, double *mag, double *sigma);
void rand_disc(double *a, double *b);
void rand_shell(double *x, double *y, double *z);
void rand_sphere(double *x, double *y, double *z);
double Vectorrandomelement(const Vector v);

#endif
