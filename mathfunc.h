/* cassbeam - a Cassegrain antenna simulator
    Copyright (C) August 18, 2003  Walter Brisken

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see
    <http://www.gnu.org/licenses/>.
 */
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
