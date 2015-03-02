#ifndef __MATHFUNC_H__
#define __MATHFUNC_H__

#include <math.h>

#ifndef M_SQRT3
#define M_SQRT3 1.73205075688772935
#endif

#ifndef C_EULER
#define C_EULER 0.57721566490153286
#endif


double pythag(double a, double b);

double bessel_J0(double x);
double bessel_J1(double x);
double bessel_Y0(double x);
double bessel_Y1(double x);
double bessel_I0(double x);
double logbessel_I0(double x);
double bessel_I1(double x);
double bessel_K0(double x);
double bessel_K1(double x);

double cubic(double p, double q, double r);
int quarticroots(double a, double b, double c, double d, double *rts);

#endif
