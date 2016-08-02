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
