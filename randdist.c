#include <math.h>
#include <glib.h>
#include <stdio.h>
#include "randdist.h"
#include "vector.h"
#include "image-vector.h"
#include "mathfunc.h"


/* return random number in [0, 1] */

double rand_one()
{
	return rand()/(double)RAND_MAX;
}


/* return random number in [-1, 1] */

double rand_pm_one()
{
	return 2.0*rand()/(double)RAND_MAX - 1.0;
}


/* return a random number with gaussian statistics */

double rand_gauss()
{
	static int iset=0;
	static double gset;
	double fac, rsq, v1, v2;

	if(iset)
	{
		iset = 0;
		return gset;
	}
	else
	{
		iset = 1;
		do
		{
			v1 = 2.0*rand()/(double)RAND_MAX - 1.0;
			v2 = 2.0*rand()/(double)RAND_MAX - 1.0;
			rsq = v1*v1 + v2*v2;
		} while(rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		return v2*fac;
	}
}


/* return a random number with Rice statistics 
 * for a good reference, see 
 * http://www.ruca.ua.ac.be/visielab/theses/sijbers/html/node22.html */

double rand_rice(double m, double sigma)
{
	double fac, x, y, v1, v2, rsq;
	
	do
	{
		v1 = 2.0*rand()/(double)RAND_MAX - 1.0;
		v2 = 2.0*rand()/(double)RAND_MAX - 1.0;
		rsq = v1*v1 + v2*v2;
	} while(rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	
	x = fac*v1*sigma + m;
	y = fac*v2*sigma;

	return sqrt(x*x + y*y);
}

double rice_moment(double m, double sigma, int moment)
{
	double r2;
	switch(moment)
	{
		case 0:
			return 1.0;
		case 1:
			r2 = m*m/(sigma*sigma);
			return sigma * sqrt(0.5*M_PI) * exp(-0.25*r2)*
				( (1.0 + 0.5*r2)*bessel_I0(0.25*r2) +
				  0.5*r2*bessel_I1(0.25*r2) );
		case 2:
			return m*m + 2.0*sigma*sigma;
		case 3:
			r2 = m*m/(sigma*sigma);
			return sigma*sigma*sigma*sqrt(0.5*M_PI)*exp(-0.25*r2)*
				( (3.0+3.0*r2+0.5*r2*r2)*bessel_I0(0.25*r2) +
				  (2.0*r2 + 0.5*r2*r2)*bessel_I1(0.25*r2) );
		case 4:
			return m*m*m*m + 8.0*sigma*sigma*(m*m + sigma*sigma);
		default:
			g_assert(0);
	}
	return -1.0;
}

int rice_params2a(const Vector x, const Vector y, double *mag, double *sigma)
{
	const int Q = 100;
	Vector ave, dif;
	int p, q, i, n;
	Matrix M;
	double v, m, mm, s, ss, dm, ds, dx, dy;

	g_assert(x);
	g_assert(y);
	n = VectorSize(x);
	g_assert(n == VectorSize(y));

	if(n < 2) return 0;

	/* try to get reasonable starting values, accurate at snr=infinity */
	m = s = ss = 0.0;
	for(i = 0; i < n; i++)
	{
		v = x[i]*x[i] + y[i]*y[i];
		s  += sqrt(v);
		ss += v;
	}

	s  /= n;
	ss /= n;
	
	m = s;
	s = sqrt(ss - s*s);

	ave = newVector(n-1);
	dif = newVector(n-1);
	for(i = 0; i < n-1; i++)
	{
		dx = x[i]+x[i+1];
		dy = y[i]+y[i+1];
		ave[i] = 0.5*sqrt(dx*dx + dy*dy);
		dx = x[i]-x[i+1];
		dy = y[i]-y[i+1];
		dif[i] = sqrt(dx*dx + dy*dy);
	}

	M = newMatrix(Q, Q);
	zeroMatrix(M);
	v = 3.0/Q;

	dm = m*v;
	ds = s*v;

	for(p = 0; p < Q; p++) for(q = 0; q < Q; q++)
	{
		mm = (p+1)*dm;
		ss = (q+1)*ds;

		for(i = 0; i < n-1; i++)	
		{
			M[p][q] += log(ave[i]/(ss*ss));
			M[p][q] += log(dif[i]/(ss*ss));
			M[p][q] += logbessel_I0(2.0*mm*ave[i]/(ss*ss));
			M[p][q] +=  (-1.0/(2.0*ss*ss))*
				(0.5*dif[i]*dif[i] + 2.0*ave[i]*ave[i] 
					+ 2.0*mm*mm);
		}
	}

	Matrixmaxvalue(M, &p, &q);
	//printf("Inital peak up = %f %f\n", (p+1)*dm, (q+1)*ds);

	Matrixpeakup(M, p, q, 1, 1, &mm, &ss, 0, 0, 0);
	if(fabs(mm-p) >= 1.0 || fabs(ss-q) >= 1.0)
	{
		m = (p+1)*dm;
		s = (q+1)*ds;
	}
	else
	{
		m = (mm+1.0)*dm;
		s = (ss+1.0)*ds;
	}

	//printf("Second peak up = %f %f\n", m, s);

	applyfunctoMatrix(M, cbrt);
	saveMatrixaspgm(M, "Pmap.pgm");
	deleteMatrix(M);

	deleteVector(dif);
	deleteVector(ave);

	if(mag) *mag = m;
	if(sigma) *sigma = s;

	return 1;
}

int rice_params2(const Vector x, const Vector y, double *mag, double *sigma)
{
	const int niter = 500;
	Vector ave, dif;
	int i, n, iter;
	double L, L0, v, m, mm, s, ss, dm, ds, dx, dy;

	g_assert(x);
	g_assert(y);
	n = VectorSize(x);
	g_assert(n == VectorSize(y));

	if(n < 2) return 0;

	/* try to get reasonable starting values, accurate at snr=infinity */
	m = s = ss = 0.0;
	for(i = 0; i < n; i++)
	{
		v = x[i]*x[i] + y[i]*y[i];
		s  += sqrt(v);
		ss += v;
	}

	s  /= n;
	ss /= n;
	
	m = s;
	s = sqrt(ss - s*s);

	ave = newVector(n-1);
	dif = newVector(n-1);
	for(i = 0; i < n-1; i++)
	{
		dx = x[i]+x[i+1];
		dy = y[i]+y[i+1];
		ave[i] = 0.5*sqrt(dx*dx + dy*dy);
		dx = x[i]-x[i+1];
		dy = y[i]-y[i+1];
		dif[i] = sqrt(dx*dx + dy*dy);
	}

	L0 = 0;

	v = 0.8;

	for(i = 0; i < n-1; i++)	
	{
		L0 += log(ave[i]/(s*s));
		L0 += log(dif[i]/(s*s));
		L0 += logbessel_I0(2.0*m*ave[i]/(s*s));
		L0 += (-1.0/(2.0*s*s))*
			(0.5*dif[i]*dif[i] + 2.0*ave[i]*ave[i] 
				+ 2.0*m*m);
	}

	for(iter = 0; iter < niter; iter++)
	{
		L = 0;

		while( (dm = v*m*rand_gauss()) + m <= 0.0);
		while( (ds = v*s*rand_gauss()) + s <= 0.0);
		v *= 0.998;

		ss = s+ds;
		mm = m+dm;
		
		for(i = 0; i < n-1; i++)	
		{
			L += log(ave[i]/(ss*ss));
			L += log(dif[i]/(ss*ss));
			L += logbessel_I0(2.0*mm*ave[i]/(ss*ss));
			L += (-1.0/(2.0*ss*ss))*
				(0.5*dif[i]*dif[i] + 2.0*ave[i]*ave[i] 
					+ 2.0*mm*mm);
		}

		if(L > L0)
		{
			L0 = L;
			s = ss;
			m = mm;
		}
	}

	deleteVector(dif);
	deleteVector(ave);

	if(mag) *mag = m;
	if(sigma) *sigma = s;

	return 1;
}


int rice_params(const Vector data, double *mag, double *sigma)
{
	double m, s, mm, ss, dm, ds, v, M1, M2, d1, d2;
	int i, j, iter, p=20, q=2, a, b;
	Matrix M, S;

	M = newMatrix(2*p+1, 2*p+1);
	S = refsubMatrix(M, q, q, 2*p-q, 2*p-q);

	/* try to get reasonable starting values, accurate at snr=infinity */
	VectormeanRMS(data, &m, &s);

	M1 = log(Vectormoment(data, 1)/VectorSize(data));
	M2 = log(Vectormoment(data, 2)/VectorSize(data));
	
	for(iter = 0; iter < 50; iter++)
	{
		if(m < 0.0) m = -m;
		if(s < 0.0) s = -s;
		v = 0.3/(p*sqrt(iter+4.0));
		dm = m*v;
		ds = s*v;
		printf("Iter = %d m = %f, s = %f", iter, m, s);
		for(i = -p; i <= p; i++) for(j = -p; j <= p; j++)
		{
			d1 = log(rice_moment(m+i*dm, s+j*ds, 1)) - M1;
			d2 = log(rice_moment(m+i*dm, s+j*ds, 2)) - M2;
			M[i+p][j+p] = d1*d1 + d2*d2;
		}

		Matrixminvalue(S, &a, &b);		

		/* go from S coords to M coords */
		a+=q;
		b+=q;

		printf("  minvalue = %d %d", a, b);
	
		v = Matrixpeakup(M, a, b, q, q, &mm, &ss, 0, 0, 0);

		printf("  peakup = %f %f", mm, ss);

		if(fabs(mm - a) > 1.9) mm = a;
		if(fabs(ss - b) > 1.9) ss = b;

		printf("  peak = %f\n", v);

		m += (mm-p)*dm;
		s += (ss-p)*ds;
	}
	deleteMatrix(M);
	deleteMatrix(S);

	if(mag) *mag = m;
	if(sigma) *sigma = s;

	return 1;
}



/* return random (a, b) such that a^2 + b^2 <= 1 */

void rand_disc(double *a, double *b)
{
	for(;;)
	{
		*a = 2.0*rand()/(double)RAND_MAX - 1.0;
		*b = 2.0*rand()/(double)RAND_MAX - 1.0;
		if((*a)*(*a) + (*b)*(*b) <= 1.0) return;
	}
}


/* return random (x, y, z) such that x*x + y*y + z*z = 1 */

void rand_shell(double *x, double *y, double *z)
{
	double r;

	for(;;)
	{
		*x = 2.0*rand()/(double)RAND_MAX - 1.0;
		*y = 2.0*rand()/(double)RAND_MAX - 1.0;
		*z = 2.0*rand()/(double)RAND_MAX - 1.0;
		r = (*x)*(*x) + (*y)*(*y) + (*z)*(*z);
		if(r <= 1.0 && r > 0.0) break;
	}
	r = sqrt(r);
	*x /= r;
	*y /= r;
	*z /= r;
}

void rand_sphere(double *x, double *y, double *z)
{
	for(;;)
	{
		*x = 2.0*rand()/(double)RAND_MAX - 1.0;
		*y = 2.0*rand()/(double)RAND_MAX - 1.0;
		*z = 2.0*rand()/(double)RAND_MAX - 1.0;
		if((*x)*(*x) + (*y)*(*y) + (*z)*(*z) <= 1.0) return;
	}
}

double Vectorrandomelement(const Vector v)
{
	return v[rand_int(VectorSize(v))];
}
