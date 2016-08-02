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
#ifndef __POLYGON_H__
#define __POLYGON_H__

#include "vector.h"

/* returns 1 if x0, y0 is inside C(x, y),
 *	   0 if x0, y0 is on C(x, y),
 *	  -1 if x0, y0 is outside C(x, y).
 */
int polygonside(const Vector x, const Vector y, double x0, double y0);

double polygonperimeter(const Vector x, const Vector y);

/* positive for counter clockwise polygon */
double polygonarea(const Vector x, const Vector y);

#endif
