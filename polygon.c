#include <glib.h>
#include <math.h>
#include "polygon.h"

/* returns 	-1 if segment passes through 0,0, 
 *		 0 if segment doesn't intersect +x axis
 *		 1 if segment touches +x axis over interval [0,1) 
 */
static int testsegment(double x0, double y0, double x1, double y1)
{
	double dx, dy, t;

	/* Check for pathological cases: */

	if(x0 == 0.0 && y0 == 0.0) return -1;

	if(y0 == y1)
	{
		if(y0 == 0.0)
		{
			if(x0*x1 < 0.0) return -1;
			else return 0;
		}
		else return 0;
	}

	/* Now test the common cases: */

	if(y0*y1 < 0.0)
	{
		dx = x1-x0;
		dy = y1-y0;
		t = -y0/dy;
		if(t < 0.0 || t >= 1.0) return 0;
		t = x0 + t*dx;
		if(t > 0.0) return 1;
		if(t == 0.0) return -1;
	}

	return 0;
}

/* returns 1 if x0, y0 is inside C(x, y),
 *	   0 if x0, y0 is on C(x, y),
 *	  -1 if x0, y0 is outside C(x, y).
 */
int polygonside(const Vector x, const Vector y, double x0, double y0)
{
	int i, n, a, b=0;

	g_assert(x);
	g_assert(y);

	n = VectorSize(x);
	g_assert(VectorSize(y) == n);

	for(i = 0; i < n; i++)
	{
		if(i == 0) 
			a = testsegment(x[n-1]-x0, y[n-1]-y0, 
					x[0]  -x0, y[0]  -y0);
		else 
			a = testsegment(x[i-1]-x0, y[i-1]-y0,
					x[i]  -x0, y[i]  -y0);
		if(a == -1) return 0; /* point on polygon */
		b += a;
	}

	if(b%2 == 1) return 1;
	return -1;
}

double polygonperimeter(const Vector x, const Vector y)
{
	int i, n;
	double x0, y0, dx, dy, p=0.0;

	g_assert(x);
	g_assert(y);

	n = VectorSize(x);
	g_assert(VectorSize(y) == n);

	x0 = x[n-1];
	y0 = y[n-1];

	for(i = 0; i < n; i++)
	{
		dx = x[i] - x0;
		x0 = x[i];
		dy = y[i] - y0;
		y0 = y[i];
		p += sqrt(dx*dx + dy*dy);
	}

	return p;
}

/* positive for counter clockwise polygon */
double polygonarea(const Vector x, const Vector y)
{
	int i, n;
	double x0, y0, x1, y1, a=0.0;

	g_assert(x);
	g_assert(y);

	n = VectorSize(x);
	g_assert(VectorSize(y) == n);

	x0 = x[n-1];
	y0 = y[n-1];

	for(i = 0; i < n; i++)
	{
		x1 = x[i];
		y1 = y[i];

		a += (y0+y1)*(x0-x1);

		x0 = x1;
		y0 = y1;
	}

	return 0.5*a;
}
