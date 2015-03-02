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
