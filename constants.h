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
#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <math.h>

/* X_Y = multiplicitive constant required to convert a number in X units to
 * a number in Y units.
 *
 * 39.37 * INCH_METER = 1.0 
 */

#define METER_INCH	39.37008
#define INCH_METER	(1.0/METER_INCH)

#define NS_METER	0.299792458	/* Exact */
#define METER_NS	(1.0/NS_METER)

#define DEG_RAD		M_PI/180.0
#define RAD_DEG		180.0/M_PI

#endif
