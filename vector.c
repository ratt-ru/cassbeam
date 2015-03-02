#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include "vector.h"
#include "randdist.h"

static int totalvectordata = 0;

int gettotalVectordata()
{
	return totalvectordata;
}

/* routine that given matrix a and vector b, repaces a by its inverse, */
/* and replaces b with [a]^-1 * b */
/* Returns 1 on success, 0 on failure */
/* Adapted from Numerical Recipes */

int gaussjordan(Matrix a, Vector b)
{
        double big, dum, pivinv, tmp;
        static int *indxc, *indxr, *ipiv;
	static int oldn = 0;
 	int i, icol = 0, irow = 0, j, k, l, ll;
	int m, n, p;

	g_assert(a != 0);
	g_assert(b != 0);
 
	n = MatrixSize1(a);
	m = MatrixSize2(a);
	p = Vectorlength(b);

	g_assert(m == n);
	g_assert(m == p);

	if(n != oldn)
	{
		oldn = n;
		if(oldn > 0)
		{
			g_free(indxc);
			g_free(indxr);
			g_free(ipiv);
		}
		if(n > 0)
		{
			indxc = g_new(int, n);
       			indxr = g_new(int, n);
        		ipiv  = g_new(int, n);
		}
	}
	if(n == 0) return 1;
        for(j = 0; j < n; j++) ipiv[j] = 0;
        for(i = 0; i < n; i++) 
	{
                big = 0.0;
                for(j = 0; j < n; j++) if(ipiv[j] != 1) for(k=0; k < n; k++)
		{
			if(ipiv[k] == 0)
			{
				if(fabs(a[j][k]) >= big) 
				{
					big = fabs(a[j][k]);
					irow=j;
					icol=k;
				}
			} 
			else if(ipiv[k] > 1) return -1;
		}
                ipiv[icol]++;
                if(irow != icol)
		{
                        for(l=0; l < n; l++) 
			{
				tmp = a[irow][l];
				a[irow][l] = a[icol][l];
				a[icol][l] = tmp;
			}
			tmp = b[irow];
			b[irow] = b[icol];
			b[icol] = tmp;
                }
                indxr[i] = irow;
                indxc[i] = icol;
                if(a[icol][icol] == 0.0) return 0;	
                pivinv = 1.0/a[icol][icol];
                a[icol][icol] = 1.0;
                for(l = 0; l < n; l++) a[icol][l] *= pivinv;
                b[icol] *= pivinv;
                for(ll = 0; ll < n; ll++) if (ll != icol)
		{
			dum = a[ll][icol];
			a[ll][icol] = 0.0;
			for(l = 0; l < n; l++) a[ll][l] -= a[icol][l]*dum;
			b[ll] -= b[icol]*dum;
		}
        }

	return 1;
}

/* Routine to fit an (n-1)-th order polynomial to (X, Y) data 
 * eg. if order=3, then the 3 element output vector contains
 * coefficients for x^0, x^1, and X^2
 */
Vector polyfit(const Vector X, const Vector Y, int order)
{
	Matrix M;
	Vector V;
	Vector f;
	int i, j, k, n;
	int res;

	if(!VectorSizeSame(X, Y))
	{
		fprintf(stderr, "polyfit : vectors not of same size\n");
		return 0;
	}
	n = VectorSize(X);
	if(n < order)
	{
		fprintf(stderr, "polyfit : system undetermined\n");
		return 0;
	}

	V = newVector(order);
	zeroVector(V);

	M = newMatrix(order, order);
	zeroMatrix(M);

	f = newVector(order);

	for(k = 0; k < n; k++)
	{
		f[0] = 1;
		for(i = 1; i < order; i++) f[i] = f[i-1]*X[k];

		for(i = 0; i < order; i++) V[i] += Y[k]*f[i];
		for(i = 0; i < order; i++) for(j = 0; j < order; j++)
			M[i][j] += f[i]*f[j];
	}

	res = gaussjordan(M, V);

	deleteMatrix(M);
	deleteVector(f);

	if(res < 0)
	{
		deleteVector(V);
		fprintf(stderr, "polyfit : system not invertable\n");
		return 0;
	}

	return V;
}

Vector weightedpolyfit(const Vector X, const Vector Y, const Vector W, 
	int order)
{
	Matrix M;
	Vector V;
	Vector f;
	int i, j, k, n;
	int res;

	if(!VectorSizeSame(X, Y) || !VectorSizeSame(X, W))
	{
		fprintf(stderr, "weightedpolyfit : vectors not of same size\n");
		return 0;
	}
	n = VectorSize(X);
	if(n < order)
	{
		fprintf(stderr, "weightedpolyfit : system undetermined\n");
		return 0;
	}

	V = newVector(order);
	zeroVector(V);

	M = newMatrix(order, order);
	zeroMatrix(M);

	f = newVector(order);

	for(k = 0; k < n; k++)
	{
		f[0] = 1;
		for(i = 1; i < order; i++) f[i] = f[i-1]*X[k];

		for(i = 0; i < order; i++) V[i] += W[k]*Y[k]*f[i];
		for(i = 0; i < order; i++) for(j = 0; j < order; j++)
			M[i][j] += W[k]*f[i]*f[j];
	}

	res = gaussjordan(M, V);

	deleteMatrix(M);
	deleteVector(f);

	if(res < 0)
	{
		deleteVector(V);
		fprintf(stderr, "weightedpolyfit : system not invertable\n");
		return 0;
	}

	return V;
}

/* Routine to fit an (n-1)-th order polynomial to (X[i], Y[i-lag]) data 
 * eg. if order=3, then the 3 element output vector contains
 * coefficients for x^0, x^1, and X^2
 */
Vector polyfitwithlag(const Vector X, const Vector Y, int order, int lag)
{
	Matrix M;
	Vector V;
	Vector f;
	int i, j, k, n;
	int res;
	int A, B;

	if(!VectorSizeSame(X, Y))
	{
		fprintf(stderr, "polyfit : vectors not of same size\n");
		return 0;
	}
	n = VectorSize(X)-abs(lag);
	if(n < order || order < 1)
	{
		fprintf(stderr, "polyfit : system undetermined\n");
		return 0;
	}

	if(lag > 0) 
	{
		A = lag; 
		B = VectorSize(X);
	}
	else 
	{
		A = 0;
		B = n;
	}
	
	V = newVector(order);
	zeroVector(V);

	M = newMatrix(order, order);
	zeroMatrix(M);

	f = newVector(order);

	for(k = A; k < B; k++)
	{
		f[0] = 1;
		for(i = 1; i < order; i++) f[i] = f[i-1]*X[k];

		for(i = 0; i < order; i++) V[i] += Y[k-lag]*f[i];
		for(i = 0; i < order; i++) for(j = 0; j < order; j++)
			M[i][j] += f[i]*f[j];
	}

	res = gaussjordan(M, V);

	deleteMatrix(M);
	deleteVector(f);

	if(res < 0)
	{
		deleteVector(V);
		fprintf(stderr, "polyfit : system not invertable\n");
		return 0;
	}

	return V;
}



/* Vector routines */

Vector newVector(int n)
{
	Vector v;

	v = g_new(VectorType, n+1);
	v++;
	VectorSize(v) = n;

	totalvectordata += n;

	return v;
}

Vector newVectorfromarray(int n, const double *arr)
{
	Vector v;
	int i;

	v = g_new(VectorType, n+1);
	v++;
	VectorSize(v) = n;

	totalvectordata += n;

	for(i = 0; i < n; i++) v[i] = arr[i];

	return v;
}

Vector newVectorfromintVector(intVector v)
{
	Vector u;
	int i, n;
	
	n = intVectorSize(v);

	u = newVector(n);
	
	for(i = 0; i < n; i++) u[i] = v[i];

	return u;
}

void zeroVector(Vector v)
{
	int i, m;

	m = VectorSize(v);
	for(i = 0; i < m; i++) v[i] = 0.0;
}

void fillVector(Vector v, double value)
{
	int i, m;

	m = VectorSize(v);
	for(i = 0; i < m; i++) v[i] = value;
}

Vector dupVector(const Vector v)
{
	int i, m;
	Vector v2;

	m = VectorSize(v);
	v2 = newVector(m);
	for(i = 0; i < m; i++) v2[i] = v[i];

	return v2;
}

void deleteVector(Vector v)
{
	totalvectordata -= VectorSize(v);
	g_free(v-1);
}

int Vectorisfinite(const Vector v)
{
	int i, m;

	m = VectorSize(v);
	for(i = 0; i < m; i++) if(!finite(v[i])) return 0;

	return 1;
}

double dotVectors(const Vector v1, const Vector v2)
{
	int i, m, n;
	double sum = 0;

	m = VectorSize(v1);
	n = VectorSize(v2);
	g_assert(m == n);

	for(i = 0; i < m; i++) sum += v1[i]*v2[i];
	
	return sum;
}

void subfromVector(Vector a, const Vector b)
{
	int i, m, n;

	m = VectorSize(a);
	n = VectorSize(b);
	g_assert(m == n);

	for(i = 0; i < m; i++) a[i] -= b[i];
}

void addtoVector(Vector a, const Vector b)
{
	int i, m;

	g_assert(VectorSizeSame(a, b));
	
	m = VectorSize(a);

	for(i = 0; i < m; i++) a[i] += b[i];
}

void addscaledVectortoVector(Vector a, const Vector b, double scale)
{
	int i, m;

	if(scale == 0.0) return;
	
	g_assert(VectorSizeSame(a, b));
	
	m = VectorSize(a);

	for(i = 0; i < m; i++) a[i] += scale*b[i];
}

Vector addVectors(const Vector a, const Vector b)
{
	int i, m;
	Vector Sum;

	g_assert(VectorSizeSame(a, b));
	
	m = VectorSize(a);

	Sum = newVector(m);
	for(i = 0; i < m; i++) Sum[i] = a[i] + b[i];

	return Sum;
}

Vector subVectors(const Vector a, const Vector b)
{
	int i, m;
	Vector Dif;

	g_assert(VectorSizeSame(a, b));
	
	m = VectorSize(a);

	Dif = newVector(m);
	for(i = 0; i < m; i++) Dif[i] = a[i] - b[i];

	return Dif;
}

void copytoVector(Vector a, const Vector b)
{
	int i, m;

	g_assert(a);
	g_assert(b);
	g_assert(VectorSizeSame(a, b));

	m = VectorSize(a);

	for(i = 0; i < m; i++) a[i] = b[i];
}

void scaleVector(Vector v, double f)
{
	int i, m;

	m = VectorSize(v);

	for(i = 0; i < m; i++) v[i] *= f;
}

void biasVector(Vector v, double f)
{
	int i, m;

	m = VectorSize(v);

	for(i = 0; i < m; i++) v[i] += f;
}

void squareVector(Vector v)
{
	int i, m;

	m = VectorSize(v);

	for(i = 0; i < m; i++) v[i] = v[i]*v[i];
}

void normalizeVector(Vector v)
{
	scaleVector(v, 1.0/sqrt(Vectorsumsquare(v)));
}

void printVector(const Vector v)
{
	int i, m;

	m = VectorSize(v);

	for(i = 0; i < m; i++)
	{
		if(i != 0) printf(" ");
		printf("%e", v[i]);
	}
	printf("\n");
}

void applyfunctoVector(Vector v, double (*func)(double x))
{
	int i, m;

	m = VectorSize(v);

	for(i = 0; i < m; i++) v[i] = func(v[i]);
}

VectorType Vectormax(const Vector v)
{
	VectorType max;
	int i, m;

	m = VectorSize(v);

	max = v[0];
	for(i = 1; i < m; i++) if(v[i] > max) max = v[i];

	return max;
}

VectorType Vectormin(const Vector v)
{
	VectorType min;
	int i, m;

	m = VectorSize(v);

	min = v[0];
	for(i = 1; i < m; i++) if(v[i] < min) min = v[i];

	return min;
}

int Vectorfindmax(const Vector v)
{
	VectorType max;
	int i, m, best;

	m = VectorSize(v);

	max = v[0];
	best = 0;
	for(i = 1; i < m; i++) if(v[i] > max) 
	{
		max = v[i];
		best = i;
	}

	return best;
}

int Vectorfindmin(const Vector v)
{
	VectorType min;
	int i, m, best;

	m = VectorSize(v);

	min = v[0];
	best = 0;
	for(i = 1; i < m; i++) if(v[i] < min) 
	{
		min = v[i];
		best = i;
	}

	return best;
}

void Vectorminmax(const Vector v, VectorType *min, VectorType *max)
{
	int i, m;

	m = VectorSize(v);

	*min = *max = v[0];

	for(i = 1; i < m; i++)
	{
		if(v[i] < *min) *min = v[i];
		if(v[i] > *max) *max = v[i];
	}
}

VectorType interpolateVector(const Vector V, double index)
{
	int N, i;
	double f;

	g_assert(V);

	N = VectorSize(V);
	if(index <= 0.0) return V[0];
	if(index >= N-1) return V[N-1];

	i = (int)index;
	f = index - i;

	return (1.0-f)*V[i] + f*V[i+1];
}

#if 0
/* determine the interpolated index for a value */
double uninterpolateVector(const Vector V, double value)
{
	double max, min;
	int i, N;

	g_assert(V);
	
	N = VectorSize(V);
	if(
}
#endif

void saveVectorasascii(const Vector v, const char *filename)
{
	FILE *out;
	int i, m;

	g_assert(v);
	out = fopen(filename, "w");
	g_assert(out);
	m = VectorSize(v);

	for(i = 0; i < m; i++) fprintf(out, "%d %f\n", i, v[i]);

	fclose(out);
}

void saveVectorasformattedascii(const Vector v, const char *filename,
	const char *format)
{
	FILE *out;
	int i, m;

	g_assert(v);
	out = fopen(filename, "w");
	g_assert(out);
	m = VectorSize(v);

	for(i = 0; i < m; i++) 
	{
		fprintf(out, "%d ", i);
		fprintf(out, format, v[i]);
		fprintf(out, "\n");
	}

	fclose(out);
}

void Vectorsavebinary(const Vector v, const char *filename)
{
	FILE *out;
	int m;

	g_assert(v);
	if(strcmp(filename, "-") == 0) out = stdout;
	else out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "Vectorsavebinary : Error opening %s for write\n", filename);
		return;
	}
	m = VectorSize(v);
	fwrite(&m, sizeof(int), 1, out);
	fwrite(v, sizeof(VectorType), m, out);
	if(strcmp(filename, "-") != 0) fclose(out);
}

Vector Vectorloadbinary(const char *filename)
{
	FILE *in;
	int m;
	Vector v;

	if(strcmp(filename, "-") == 0) in = stdin;
	else in = fopen(filename, "r");
	if(!in)
	{
		fprintf(stderr, "Vectorloadbinary : Error opening %s for read\n", filename);
		return 0;
	}
	if(fread(&m, sizeof(int), 1, in) < 1)
	{
		fprintf(stderr, "Vectorloadbinary : Error reading from %s\n", filename);
		return 0;
	}
	g_assert(m > 0);
	v = newVector(m);
	g_assert(v);
	if(fread(v, sizeof(VectorType), m, in) < m)
	{
		fprintf(stderr, "WARNING: Vectorloadbinary : partial vector read from %s\n", filename);
	}
	if(strcmp(filename, "-") != 0) fclose(in);
	return v;
}

void VectormeanRMS(const Vector v, double *mean, double *rms)
{
	double s=0.0, ss=0.0;
	int i, m;
	
	g_assert(v);
	m = VectorSize(v);

	for(i = 0; i < m; i++)
	{
		s += v[i];
		ss += v[i]*v[i];
	}
	s  /= (float)m;	
	ss /= (float)m;
	*mean = s;
	*rms = sqrt(ss-s*s);
}

double Vectorsum(const Vector v)
{
	double s=0.0;
	int i, m;
	
	g_assert(v);
	m = VectorSize(v);

	for(i = 0; i < m; i++) s += v[i];

	return s;	
}

double Vectorsumsquare(const Vector v)
{
	double s=0.0;
	int i, m;
	
	g_assert(v);
	m = VectorSize(v);

	for(i = 0; i < m; i++) s += v[i]*v[i];

	return s;	
}

double Vectormean(const Vector v)
{
	double s=0.0;
	int i, m;
	
	g_assert(v);
	m = VectorSize(v);

	for(i = 0; i < m; i++) s += v[i];

	return s / (float)m;	
}

double Vectormeansquare(const Vector v)
{
	double s=0.0;
	int i, m;
	
	g_assert(v);
	m = VectorSize(v);

	for(i = 0; i < m; i++) s += v[i]*v[1];

	return s / (float)m;	
}

double VectorRMS(const Vector v)
{
	double mean, rms;
	
	VectormeanRMS(v, &mean, &rms);

	return rms;
}

double Vectormoment(const Vector v, int moment)
{
	int i, n;
	double s = 0.0, val;

	g_assert(v);
	g_assert(moment >= 0);
	
	n = VectorSize(v);
	if(moment == 0) return n;

	for(i = 0; i < n; i++)
	{
		val = pow(fabs(v[i]), moment);
		if(moment % 2 && v[i] < 0.0) s -= val;
		else s += val;
	}

	return s;
}

/* New vector is factor smaller than src */
Vector rebinVector(const Vector src, int factor)
{
	Vector v;
	int i, n, m;
	
	n = VectorSize(src);
	m = n/factor;
	n = m*factor;
	
	v = newVector(m);
	zeroVector(v);

	for(i = 0; i < n; i++) v[i/factor] += src[i];

	scaleVector(v, 1.0/(float)factor);

	return v;
}

void Vectorboxcar(const Vector src, Vector dest, int s)
{
	int i = 0, get, put, m=0;
	double sum = 0.0;
	int n;

	if(s <= 1) 
	{
		if(s < 1) fprintf(stderr, "Vectorboxcar : Warning : s < 1\n");
		copytoVector(dest, src);
		return;
	}

	n = VectorSize(src);
	get = 0;
	put = -s;

	for(;;)
	{
		if(get >= 0 && get < n)
		{
			sum += src[get];
			m++;
		}
		if(put >= 0 && put < n)
		{
			sum -= src[put];
			m--;
		}
		get++;
		put++;
		i = (get + put)/2;
		if(i < 0) continue;
		if(i >= n) break;
		dest[i] = sum/(double)m;
	}	
}

double VectorAllanStandardDeviation(const Vector v, int lag, int pad, int step)
{
	double ss = 0.0;
	double delta;
	int i = pad;
	int n;
	int m = 0;

	g_assert(v);

	n = VectorSize(v);

	if(step < 1) step = 1;

	for(i = pad; i < n - pad - lag; i += step)
	{
		delta = v[i] - v[i+lag];
		ss += delta*delta;
		m++;
	}

	return sqrt(0.5*ss/(double)m);
}

intVector ComputeAllanLags(int n, double power)
{
	intVector lags;
	double S;
	int lag, lastlag, m;
	
	lastlag = m = 0;
	for(S = 1.0; S < n/3.0; S*=power)
	{
		lag = (int)(S + 0.5);
		if(lag == lastlag) continue;
		m++;
		lastlag = lag;
	}
	lags = newintVector(m);
	lastlag = 0;
	m = 0;
	for(S = 1.0; S < n/3.0; S*=power)
	{
		lag = (int)(S + 0.5);
		if(lag == lastlag) continue;
		lags[m++] = lag;
		lastlag = lag;
	}
	
	return lags;
}

Vector VectorAllanStandardDeviations(const Vector v, const intVector lags)
{
	Vector allanvar;
	Vector boxcar;
	int n, l, i, s;

	n = VectorSize(v);
	l = intVectorSize(lags);

	boxcar = newVector(n);
	allanvar = newVector(l);

	for(i = 0; i < l; i++)
	{
		s = lags[i];
		Vectorboxcar(v, boxcar, s);
		allanvar[i] = VectorAllanStandardDeviation(
			boxcar, s, (s+1)/2, (s+3)/4);
	}

	deleteVector(boxcar);

	return allanvar;
}

int stringtoVector(const char *str, Vector v)
{
	int n=0, l, i=0;
	double d;
	gchar *s;

	s = g_strdup(str);

	for(i = 0; s[i]; i++) 
	{
		if(s[i] >= '0' && s[i] <= '9') continue;
		if(s[i] == '.' || s[i] == '+' || s[i] == '-') continue;
		if((s[i] == 'e' || s[i] == 'E') &&
		   (s[i+1] == '+' || s[i+1] == '-')) continue;
		s[i] = ' ';
	}
	i = 0;

	while(sscanf(s+n, "%lf%n", &d, &l) > 0)
	{
		if(v) if(i < VectorSize(v)) v[i] = d;
		i++;
		n+=l;
	}

	g_free(s);

	return i;
}

Vector newVectorfromstring(const char *str)
{
	Vector v;
	int n;
	
	n = stringtoVector(str, 0);
	if(n == 0) return 0;
	
	v = newVector(n);
	stringtoVector(str, v);

	return v;
}

char *Vectortostring(const Vector v)
{
	GString *str;
	gchar *out;
	char temp[100];
	int i, n;

	g_assert(v);
	
	n = VectorSize(v);
	sprintf(temp, "%f", v[0]);
	str = g_string_new(temp);

	if(n > 1) for(i = 1; i < n; i++)
	{
		sprintf(temp, ", %f", v[i]);
		str = g_string_append(str, temp);
	}

	out = str->str;

	return out;
}

Vector *Vectorcolumnsfromfile(const char *filename, int mincols, int *ncol)
{
	FILE *in;
	char format[] = "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
	double d[16];
	int nc = 0, nr = 0, n, i, row;
	Vector *data;
	char l[1000];
	
	in = fopen(filename, "r");
	if(!in) 
	{
		if(ncol != 0) *ncol = 0;
		return 0;
	}
	for(;;)
	{
		fgets(l, 999, in);
		if(feof(in)) break;
		if(l[0] == '#' || l[0] > 57 || l[0] == '*') continue;
		n = sscanf(l , format, 	d   , d+1 , d+2 , d+3 ,
					d+4 , d+5 , d+6 , d+7 ,
					d+8 , d+9 , d+10, d+11,
					d+12, d+13, d+14, d+15);
		if(n >= mincols)
		{
			if(nc == 0) nc = n;
			if(n >= nc) nr++;
		}
	}
	fclose(in);

	data = g_new(Vector, nc+1);
	for(i = 0; i < nc; i++) data[i] = newVector(nr);
	data[nc] = 0;	/* mark the end of the list.  For code that ignores ncol */

	in = fopen(filename, "r");
	for(row = 0; row < nr;)
	{
		fgets(l, 999, in);
		if(feof(in)) break;
		if(l[0] == '#' || l[0] > 57) continue;
		n = sscanf(l , format, 	d   , d+1 , d+2 , d+3 ,
					d+4 , d+5 , d+6 , d+7 ,
					d+8 , d+9 , d+10, d+11,
					d+12, d+13, d+14, d+15);
		if(n >= nc)
		{
			if(n >= nc)
			{
				for(i = 0; i < nc; i++) data[i][row] = d[i];
				row++;
			}
		}
	}
	fclose(in);

	if(ncol != 0) *ncol = nc;

	return data;
}

int Vectorcolumnstofile(const Vector *data, const char *filename, int ncol)
{
	FILE *out;
	int i, j, n;

	if(ncol <= 0) for(ncol = 0; data[ncol]; ncol++);
	
	for(i = 0; i < ncol; i++)
	{
		if(!data[i])
		{
			fprintf(stderr, "Vectorcolumnstofile: null data "
					"column %d\n", i);
			return -1;
		}
		if(i > 0) if(!VectorSizeSame(data[0], data[i]))
		{
			fprintf(stderr, "Vectorcolumnstofile: vector size "
					"mismatch, cols %d %d\n", 0, i);
			return -1;
		}
	}

	n = VectorSize(data[0]);

	out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "Vectorcolumnstofile: error opening file %s\n",
			filename);
		return -1;
	}
	for(j = 0; j < n; j++)
	{
		for(i = 0; i < ncol; i++) fprintf(out, "%10.5f ", data[i][j]);
		fprintf(out, "\n");
	}
	fclose(out);

	return 0;
}

int Vectorcolumnstofilewithformat(const Vector *data, const char *filename, int ncol,
	const char *format)
{
	FILE *out;
	int i, j, n;

	if(ncol <= 0) for(ncol = 0; data[ncol]; ncol++);
	
	for(i = 0; i < ncol; i++)
	{
		if(!data[i])
		{
			fprintf(stderr, "Vectorcolumnstofile: null data "
					"column %d\n", i);
			return -1;
		}
		if(i > 0) if(!VectorSizeSame(data[0], data[i]))
		{
			fprintf(stderr, "Vectorcolumnstofile: vector size "
					"mismatch, cols %d %d\n", 0, i);
			return -1;
		}
	}

	n = VectorSize(data[0]);

	out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "Vectorcolumnstofile: error opening file %s\n",
			filename);
		return -1;
	}
	for(j = 0; j < n; j++)
	{
		for(i = 0; i < ncol; i++) 
		{
			fprintf(out, format, data[i][j]);
			fprintf(out, " ");
		}
		fprintf(out, "\n");
	}
	fclose(out);

	return 0;
}

int Vectorcolumnstofilewithmask(const Vector *data, const Vector mask, 
	const char *filename, int ncol)
{
	FILE *out;
	int i, j, n;

	if(ncol <= 0) for(ncol = 0; data[ncol]; ncol++);
	
	for(i = 0; i < ncol; i++)
	{
		if(!data[i])
		{
			fprintf(stderr, "Vectorcolumnstofile: null data "
					"column %d\n", i);
			return -1;
		}
		if(i > 0) if(!VectorSizeSame(data[0], data[i]))
		{
			fprintf(stderr, "Vectorcolumnstofile: vector size "
					"mismatch, cols %d %d\n", 0, i);
			return -1;
		}
	}

	n = VectorSize(data[0]);

	out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "Vectorcolumnstofile: error opening file %s\n",
			filename);
		return -1;
	}
	for(j = 0; j < n; j++)
	{
		if(mask[j] < 0.5) continue;
		for(i = 0; i < ncol; i++) fprintf(out, "%10.5f ", data[i][j]);
		fprintf(out, "\n");
	}
	fclose(out);

	return 0;
}

int Vectorcolumnstofilewithgaps(const Vector *data, const char *filename, 
	int ncol, int refcol, double gap)
{
	FILE *out;
	int i, j, n;

	if(ncol <= 0) for(ncol = 0; data[ncol]; ncol++);
	
	for(i = 0; i < ncol; i++)
	{
		if(!data[i])
		{
			fprintf(stderr, "Vectorcolumnstofile: null data "
					"column %d\n", i);
			return -1;
		}
		if(i > 0) if(!VectorSizeSame(data[0], data[i]))
		{
			fprintf(stderr, "Vectorcolumnstofile: vector size "
					"mismatch, cols %d %d\n", 0, i);
			return -1;
		}
	}

	n = VectorSize(data[0]);

	out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "Vectorcolumnstofile: error opening file %s\n",
			filename);
		return -1;
	}
	for(j = 0; j < n; j++)
	{
		if(j > 0)
		{
			if(data[refcol][j] > data[refcol][j-1]+gap)
				fprintf(out, "\n");
		}
		for(i = 0; i < ncol; i++) fprintf(out, "%10.5f ", data[i][j]);
		fprintf(out, "\n");
	}
	fclose(out);

	return 0;
}

Vector subVector(const Vector v, int n1, int n2)
{
	Vector u;
	int i, m;
	
	g_assert(v);
	m = VectorSize(v);
	
	if(n1 < 0 || n2 < n1 || n2 > m)
	{
		fprintf(stderr, "subVector : range error\n");
		return 0;
	}
	
	u = newVector(n2-n1+1);
	for(i = n1; i <= n2; i++) u[i-n1] = v[i];

	return u;
}

Vector newunitVector(const Vector v)
{
	Vector u;
	double l;

	l = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	if(l == 0.0) return 0;
	u = dupVector(v);
	u[0] /= l;
	u[1] /= l;
	u[2] /= l;

	return u;
}

int Vectorisincreasing(const Vector v)
{
	int i, N;

	g_assert(v);

	N = VectorSize(v);
	if(N < 2) return 1;

	for(i = 1; i < N; i++) if(v[i] < v[i-1]) return 0;
	return 1;
}

int Vectorisdecreasing(const Vector v)
{
	int i, N;

	g_assert(v);

	N = VectorSize(v);
	if(N < 2) return 1;

	for(i = 1; i < N; i++) if(v[i] > v[i-1]) return 0;
	return 1;
}

/* Matrix Routines */

/* newMatrix returns matrix with dimensions [n][m] */

Matrix newpaddedMatrix(int n, int m, int rowpad)
{
	Matrix M;
	MatrixTypePointer D;
	int j, w;

	w = m+rowpad;
	D = g_new(MatrixType, w*n+MATRIXVALUESPERBLOCK)+MATRIXVALUESPERBLOCK; 
	M = g_new(MatrixTypePointer, ((n>1) ? (n+MATRIXSTUBINDICES) : 
					      (2+MATRIXSTUBINDICES)));
	M += MATRIXSTUBINDICES;
	MatrixStub(M)->datastart = D;
	MatrixStub(M)->n = n;
	MatrixStub(M)->m = m;
	MatrixStub(M)->rowpad = rowpad;
	Matrixrefcount(M) = 1;
	for(j = 0; j < n; j++) M[j] = D + (j*w);  /* Index array for [][] */
	if(n == 1) M[1] = D + w;

	return M;
}

void deleteMatrix(Matrix M)
{
	MatrixStubType *ms;
	
	g_assert(M);
	Matrixrefcount(M)--;
	ms = MatrixStub(M);
	if(Matrixrefcount(M) <= 0) 
	{
		if(ms->datastart-MATRIXVALUESPERBLOCK == 0)
			fprintf(stderr, "ACK!!! matrix data gone!\n");
		else
			g_free(ms->datastart-MATRIXVALUESPERBLOCK);
	}
	g_free(ms);
}

Matrix refsubMatrix(const Matrix M, int n1, int m1, int n2, int m2)
{
	Matrix N;
	int j, m, n, w;

	g_assert(M);

	n = MatrixSize1(M);
	m = MatrixSize2(M);

	if(n1 < 0) n1 = 0;
	if(n2 < 0) n2 = n-1;
	if(m1 < 0) m1 = 0;
	if(m2 < 0) m2 = m-1;

	if(n2 < n1 || m2 < m1 || m2 >= m || n2 >= n)
	{
		fprintf(stderr, "subMatrix: limit error\n");
		return 0;
	}

	w = Matrixstride(M);
	N = g_new(MatrixTypePointer, ((n>1) ? (n+MATRIXSTUBINDICES) : 
					      (2+MATRIXSTUBINDICES)));
	N += MATRIXSTUBINDICES;
	MatrixStub(N)->datastart = MatrixStub(M)->datastart;
	MatrixStub(N)->n = n2-n1+1;
	MatrixStub(N)->m = m2-m1+1;
	MatrixStub(N)->rowpad = w-(m2-m1+1);
	Matrixrefcount(M)++;
	for(j = 0; j <= n2-n1; j++) N[j] = M[n1] + m1 + (j*w);  
	if(n2 == n1) N[1] = N[0]+w;

	return N;
}

void Matrixsavebinary(const Matrix M, const char *filename)
{
	FILE *out;
	int m, n, p;

	g_assert(M);
	if(strcmp(filename, "-") == 0) out = stdout;
	else out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "Matrixsavebinary : Error opening %s for write\n", filename);
		return;
	}
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	p = Matrixrowpad(M);
	fwrite(&n, sizeof(int), 1, out);
	fwrite(&m, sizeof(int), 1, out);
	fwrite(&p, sizeof(int), 1, out);
	fwrite(M[0], sizeof(MatrixType), n*(m+p), out);
	if(strcmp(filename, "-") != 0) fclose(out);
}

Matrix Matrixloadbinary(const char *filename)
{
	FILE *in;
	int m, n, p;
	Matrix M;

	if(strcmp(filename, "-") == 0) in = stdin;
	else in = fopen(filename, "r");
	if(!in)
	{
		fprintf(stderr, "Matrixloadbinary : Error opening %s for read\n", filename);
		return 0;
	}
	if(!fread(&n, sizeof(int), 1, in) ||
	   !fread(&m, sizeof(int), 1, in) ||
	   !fread(&p, sizeof(int), 1, in))
	{
		fprintf(stderr, "Matrixloadbinary : Error reading from %s .\n", filename);
		return 0;
	}
	M = newpaddedMatrix(n, m, p);
	g_assert(M);
	if(fread(M[0], sizeof(MatrixType), n*(m+p), in) < n*(m+p))
		fprintf(stderr, "WARNING: Matrixloadbinary : partial matrix loaded from %s\n", filename);

	return M;
}

void zeroMatrix(Matrix M)
{
	int i, j, m, n;
	MatrixTypePointer Mrow;

	n = MatrixSize1(M);
	m = MatrixSize2(M);
	for(j = 0; j < n; j++) 
	{
		Mrow = M[j];
		for(i = 0; i < m; i++) Mrow[i] = 0.0;
	}
}

Matrix dupMatrix(const Matrix M)
{
	int i, j, m, n;
	Matrix M2;

	n = MatrixSize1(M);
	m = MatrixSize2(M);
	M2 = newpaddedMatrix(n, m, Matrixrowpad(M));
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) M2[j][i] = M[j][i];

	return M2;
}

void printMatrix(const Matrix M)
{
	int i, j, n, m;
	
	if(!M)
	{
		printf("printMatrix: NULL\n");
		return;
	}
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	
	printf("Matrix (%d by %d)\n", n, m);
	
	for(j = 0; j < n; j++)
	{
		for(i = 0; i < m; i++) printf("%12f ", M[j][i]);
		printf("\n");
	}
}

Matrix repadMatrix(const Matrix M, int pad)
{
	int i, j, m, n;
	Matrix M2;

	n = MatrixSize1(M);
	m = MatrixSize2(M);
	M2 = newpaddedMatrix(n, m, pad);
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) M2[j][i] = M[j][i];

	return M2;
}

MatrixType Matrixmax(const Matrix M)
{
	MatrixType max;
	int i, j, m, n;
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	max = M[0][0];
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) 
		if(M[j][i] > max) max = M[j][i];

	return max;
}

MatrixType Matrixmaxvalue(const Matrix M, int *bestj, int *besti)
{
	MatrixType max;
	int i, j, m, n, a, b;
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	max = M[0][0];
	a = b = 0;
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) 
		if(M[j][i] > max) 
		{
			max = M[j][i];
			a = i;
			b = j;
		}
	if(bestj) *bestj = b;
	if(besti) *besti = a;

	return max;
}

MatrixType Matrixmin(const Matrix M)
{
	MatrixType min;
	int i, j, m, n;
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	min = M[0][0];
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) 
		if(M[j][i] < min) min = M[j][i];

	return min;
}

MatrixType Matrixminvalue(const Matrix M, int *bestj, int *besti)
{
	MatrixType min;
	int i, j, m, n, a, b;
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	min = M[0][0];
	a = b = 0;
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) 
		if(M[j][i] < min) 
		{
			min = M[j][i];
			a = i;
			b = j;
		}
	if(bestj) *bestj = b;
	if(besti) *besti = a;

	return min;
}

void Matrixminmax(const Matrix M, MatrixType *min, MatrixType *max)
{
	int i, j, m, n;
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	*min = *max = M[0][0];

	for(i = 0; i < m; i++) for(j = 0; j < n; j++)
	{
		if(M[j][i] < *min) *min = M[j][i];
		if(M[j][i] > *max) *max = M[j][i];
	}
}

/* Find best fit to M(y, x) = Ayy (y-y0)^2 + Axy (x-x0) (y-y0) + Axx (x-x0)^2 +C
   to matrix M, using a patch of size (2 py + 1)*(2 px + 1) centered at y , x .
   The value C is returned */
double Matrixpeakup(const Matrix D, int y, int x, int py, int px,
	double *y0, double *x0, double *Ayy, double *Axy, double *Axx)
{
	int m, n, i, j, a, b, u, v;
	Vector V;
	Matrix M;
	double f[6], C;
	int res;

	g_assert(D);
	n = MatrixSize1(D);
	m = MatrixSize2(D);

	V = newVector(6);
	M = newMatrix(6, 6);
	zeroVector(V);
	zeroMatrix(M);

	f[5] = 1.0;

	for(b = -py; b <= py; b++) 
	{
		j = b + y;
		if(j < 0) j += n;
		else if(j >= n) j -= n;
		f[2] = b*b;
		f[4] = b;
		
		for(a = -px; a <= px; a++)
		{
			i = a + x;
			if(i < 0) i += m;
			else if(i >= m) i -= m;
			f[0] = a*a;
			f[1] = a*b;
			f[3] = a;
			for(u = 0; u < 6; u++) V[u] += f[u]*D[j][i];
			for(v = 0; v < 6; v++) for(u = 0; u < 6; u++)
				M[v][u] += f[u]*f[v];
		}
	}

	res = gaussjordan(M, V);
	g_assert(res == 1);

	deleteMatrix(M);

	if(Ayy) *Ayy = V[2];
	if(Axy) *Axy = V[1];
	if(Axx) *Axx = V[0];

	C = 4.0*V[0]*V[2] - V[1]*V[1];

	if(y0) *y0 = (float)y + (V[1]*V[3] - 2.0*V[0]*V[4])/C;
	if(x0) *x0 = (float)x + (V[1]*V[4] - 2.0*V[2]*V[3])/C;

	C = V[5];

	deleteVector(V);

	return C;
}

void scaleMatrix(Matrix M, double f)
{
	int i, j, m, n;
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	
	for(j = 0; j < n; j++) for(i = 0; i < m; i++) M[j][i] *= f;
}

void addscaledMatrixtoMatrix(Matrix a, const Matrix b, double scale)
{
	int i, j, m, n;

	if(scale == 0.0) return;

	g_assert(MatrixSizeSame(a, b));

	n = MatrixSize1(a);
	m = MatrixSize2(a);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		a[j][i] += scale*b[j][i];
}

void biasMatrix(Matrix M, double f)
{
	int i, j, m, n;
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	
	for(j = 0; j < n; j++) for(i = 0; i < m; i++) M[j][i] += f;
}

void addtoMatrix(Matrix M, const Matrix N)
{
	int i, j, m, n;
	
	g_assert(MatrixSizeSame(M, N));
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] += N[j][i];
}

void subfromMatrix(Matrix M, const Matrix N)
{
	int i, j, m, n;
	
	g_assert(MatrixSizeSame(M, N));
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] -= N[j][i];
}

Matrix addMatrices(const Matrix M, const Matrix N)
{
	int i, j, m, n;
	Matrix Sum;
	
	g_assert(MatrixSizeSame(M, N));
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	Sum = newMatrix(n, m);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		Sum[j][i] = M[j][i] + N[j][i];

	return Sum;
}

void copytoMatrix(Matrix M, const Matrix N)
{
	int i, j, m, n;
	
	g_assert(MatrixSizeSame(M, N));
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] = N[j][i];
}

void applyfunctoMatrix(Matrix M, double (*func)(double x))
{
	int i, j, m, n;
	
	g_assert(M);
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] = func(M[j][i]);
}

Matrix Matrixproduct(Vector v2, Vector v1)
{
	Matrix M;
	int i, j, m, n;

	g_assert(v2);
	g_assert(v1);

	M = newMatrix(Vectorlength(v2), Vectorlength(v1));
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] = v2[j]*v1[i];
	
	return M;
}

Matrix dupsubMatrix(const Matrix M, int n1, int m1, int n2, int m2)
{
	Matrix N;
	int i, j, m, n;

	g_assert(M);

	n = MatrixSize1(M);
	m = MatrixSize2(M);

	if(n1 < 0) n1 = 0;
	if(n2 < 0) n2 = n-1;
	if(m1 < 0) m1 = 0;
	if(m2 < 0) m2 = m-1;

	if(n2 < n1 || m2 < m1 || m2 >= m || n2 >= n)
	{
		fprintf(stderr, "dupsubMatrix: limit error\n");
		return 0;
	}

	N = newMatrix(n2-n1+1, m2-m1+1);
	
	for(j = n1; j <= n2; j++) for(i = m1; i <= m2; i++)
		N[j-n1][i-m1] = M[j][i];

	return N;
}

/* Copy region (m1, n1)-(m2, n2) from M to region (m3, n3) in N */
void copysubMatrix(Matrix N, const Matrix M, int n1, int m1, int n2, int m2,
	int n3, int m3)
{
	int i, j, m, n;

	g_assert(M);
	g_assert(N);

	n = MatrixSize1(M);
	m = MatrixSize2(M);

	if(n1 < 0) n1 = 0;
	if(n2 < 0) n2 = n-1;
	if(m1 < 0) m1 = 0;
	if(m2 < 0) m2 = m-1;

	if(n2 < n1 || m2 < m1 || m2 >= m || n2 >= n ||
		n3 >= MatrixSize1(N) || m3 >= MatrixSize2(N))
	{
		fprintf(stderr, "copysubMatrix: limit error\n");
		return;
	}

	if(n3 < 0) n3 = MatrixSize1(N)-(n2-n1);
	if(m3 < 0) m3 = MatrixSize2(N)-(m2-m1);

	if(n3 < 0)
	{
		n1 -= n3;
		n3 = 0;
	}
	if(m3 < 0)
	{
		m1 -= m3;
		m3 = 0;
	}
	
	if(n3+n2-n1 >= MatrixSize1(N))
		n2 = n1-n3+MatrixSize1(N)-1;
	if(m3+m2-m1 >= MatrixSize2(N))
		m2 = m1-m3+MatrixSize2(N)-1;

	for(j = n1; j <= n2; j++) for(i = m1; i <= m2; i++)
		N[j+n1-n3][i+m1-m3] = M[j][i];
}

Vector Matrixrow(const Matrix M, int row)
{
	Vector V;
	int m, n, i;
	MatrixTypePointer Mrow;

	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	if(row < 0) row = n-row;
	g_assert((row < n) && (row >= 0));
	Mrow = M[row];

	V = newVector(m);
	
	for(i = 0; i < m; i++) V[i] = Mrow[i];

	return V;
}

Vector Matrixcolumn(const Matrix M, int column)
{
	Vector V;
	int m, n, j;
	
	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	if(column < 0) column = m-column;
	g_assert((column < m) && (column >= 0));

	V = newVector(n);

	for(j = 0; j < n; j++) V[j] = M[j][column];

	return V;
}

Vector sumMatrixrows(const Matrix M)
{
	Vector V;
	int m, n, i, j;

	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	V = newVector(m);
	zeroVector(V);
	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		V[i] += M[j][i];

	return V;
}

Vector sumMatrixcolumns(const Matrix M)
{
	Vector V;
	int m, n, i, j;

	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	V = newVector(n);
	zeroVector(V);
	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		V[j] += M[j][i];

	return V;
}

void MatrixmeanRMS(const Matrix M, double *rms, double *mean)
{
	double s=0.0, ss = 0.0;
	int m, n, i, j;

	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	
	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
	{
		s  += M[j][i];
		ss += M[j][i]*M[j][i];
	}
	s  /= (float)(m*n);
	ss /= (float)(m*n);
	*rms = s;
	*mean = sqrt(ss-s*s);
}

double Matrixmean(const Matrix M)
{
	double mean, rms;
	
	MatrixmeanRMS(M, &mean, &rms);

	return mean;
}

double MatrixRMS(const Matrix M)
{
	double mean, rms;
	
	MatrixmeanRMS(M, &mean, &rms);

	return rms;
}

/* C[i][k] = sum_j A[i][j] B[j][k] */
Matrix Matrixmultiply(const Matrix A, const Matrix B)
{
	int i, j, k, l, m, n;
	Matrix C;
	
	g_assert(MatrixSize2(A) == MatrixSize1(B));

	l = MatrixSize1(A);
	m = MatrixSize2(A);
	n = MatrixSize2(B);

	C = newMatrix(l, n);
	zeroMatrix(C);
	
	for(k = 0; k < n; k++) for(j = 0; j < m; j++) for(i = 0; i < l; i++)
		C[i][k] += A[i][j]*B[j][k];

	return C;
}

/* C[i][k] = sum_j A[i][j] B[j][k] */
void copyMatrixmultiply(Matrix C, const Matrix A, const Matrix B)
{
	int i, j, k, l, m, n;
	
	g_assert(MatrixSize2(A) == MatrixSize1(B));
	g_assert(MatrixSize1(C) == MatrixSize1(A));
	g_assert(MatrixSize2(C) == MatrixSize2(B));

	l = MatrixSize1(A);
	m = MatrixSize2(A);
	n = MatrixSize2(B);

	zeroMatrix(C);
	
	for(k = 0; k < n; k++) for(j = 0; j < m; j++) for(i = 0; i < l; i++)
		C[i][k] += A[i][j]*B[j][k];
}

Vector MatrixVectormultiply(const Matrix M, const Vector V)
{
	int m, n;
	int i, j;
	Vector v;
	
	g_assert(M);
	g_assert(V);
	
	m = MatrixSize1(M);
	n = MatrixSize2(M);
	g_assert(VectorSize(V) == n);

	v = newVector(m);
	zeroVector(v);

	for(i = 0; i < m; i++) for(j = 0; j < n; j++) v[i] += M[i][j]*V[j];

	return v;
}

void Matrixmultiplyelements(Matrix M, const Matrix N)
{
	int m, n, i, j;
	
	g_assert(M);
	g_assert(N);

	if(!MatrixSizeSame(M, N))
	{
		fprintf(stderr, "Matrixmultiply: dissimilar matrices\n");
		return;
	}

	n = MatrixSize1(M);
        m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] *= N[j][i];
}

void Matrixdivideelements(Matrix M, const Matrix N)
{
	int m, n, i, j;
	
	g_assert(M);
	g_assert(N);

	if(!MatrixSizeSame(M, N))
	{
		fprintf(stderr, "Matrixdivide: dissimilar matrices\n");
		return;
	}

	n = MatrixSize1(M);
        m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] /= N[j][i];
}

Matrix Matrixcomplexamplitudes(const Matrix M)
{
	int m, n, i, j, p;
	Matrix N;

	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	p = Matrixrowpad(M);
	g_assert( (!(p%2)) && (!(m%2)) );
	
	N = newpaddedMatrix(n, m/2, p/2);
	for(j = 0; j < n; j++) for(i = 0; i < (m+p)/2; i++)
		N[j][i] = sqrt(M[j][2*i]*M[j][2*i] + M[j][2*i+1]*M[j][2*i+1]);
	
	return N;
}

Matrix Matrixcomplexphases(const Matrix M)
{
	int m, n, i, j, p;
	Matrix N;

	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);
	p = Matrixrowpad(M);
	g_assert( (!(p%2)) && (!(m%2)) );
	
	N = newpaddedMatrix(n, m/2, p/2);
	for(j = 0; j < n; j++) for(i = 0; i < (m+p)/2; i++)
		N[j][i] = atan2(M[j][2*i+1], M[j][2*i]);
	
	return N;
}


/* M[i][j] = M[i][j] * N[i][j] */
void Matrixcomplexmultiply(Matrix M, const Matrix N, int pad)
{
	int m, n, i, j;
	double re, im;
	
	g_assert(M);
	g_assert(N);

	if(!MatrixSizeSame(M, N))
	{
		fprintf(stderr, "Matrixcomplexmultiply: dissimilar matrices\n");
		return;
	}

	n = MatrixSize1(M);
        m = MatrixSize2(M);

	if(m%2)
	{
		fprintf(stderr, "Matrixcomplexmultiply: odd matrix size\n");
		return;
	}

	for(j = 0; j < n; j++) for(i = 0; i < m+pad; i+=2)
	{
		re = M[j][i]*N[j][i] - M[j][i+1]*N[j][i+1];
		im = M[j][i]*N[j][i+1] + M[j][i+1]*N[j][i];
		M[j][i] = re;
		M[j][i+1] = im;
	}
}

/* M[i][j] = M[i][j] * N*[i][j]  (M[i][j] times complex conjugate of N[i][j]) */
void Matrixcomplexconjugatemultiply(Matrix M, const Matrix N, int pad)
{
	int m, n, i, j;
	double re, im;
	
	g_assert(M);
	g_assert(N);

	if(!MatrixSizeSame(M, N))
	{
		fprintf(stderr, "Matrixcomplexmultiply: dissimilar matrices\n");
		return;
	}

	n = MatrixSize1(M);
        m = MatrixSize2(M);

	if(m%2)
	{
		fprintf(stderr, "Matrixcomplexmultiply: odd matrix size\n");
		return;
	}

	for(j = 0; j < n; j++) for(i = 0; i < m+pad; i+=2)
	{
		re = M[j][i]*N[j][i] + M[j][i+1]*N[j][i+1];
		im = -M[j][i]*N[j][i+1] + M[j][i+1]*N[j][i];
		M[j][i] = re;
		M[j][i+1] = im;
	}
}

/* M[i][j] = M[i][j] / N[i][j] */
void Matrixcomplexdivide(Matrix M, const Matrix N, int pad)
{
	int m, n, i, j;
	double re, im, mag2;
	
	g_assert(M);
	g_assert(N);

	if(!MatrixSizeSame(M, N))
	{
		fprintf(stderr, "Matrixcomplexmultiply: dissimilar matrices\n");
		return;
	}

	n = MatrixSize1(M);
        m = MatrixSize2(M);

	if(m%2)
	{
		fprintf(stderr, "Matrixcomplexmultiply: odd matrix size\n");
		return;
	}

	for(j = 0; j < n; j++) for(i = 0; i < m+pad; i+=2)
	{
		re = M[j][i]*N[j][i] + M[j][i+1]*N[j][i+1];
		im = -M[j][i]*N[j][i+1] + M[j][i+1]*N[j][i];
		mag2 = N[j][i]*N[j][i] + N[j][i+1]*N[j][i+1];
		if(mag2 == 0.0) M[j][i] = M[j][i+1] = 0.0;
		else
		{
			M[j][i] = re/mag2;
			M[j][i+1] = im/mag2;
		}
	}
}

Vector Matrixhistogram(const Matrix M, int bins)
{
	Vector V;
	MatrixType min, max;
	double delta;
	int m, n, i, j;

	g_assert(M);
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	Matrixminmax(M, &min, &max);
	V = newVector(bins);
	zeroVector(V);
	
	delta = ((float)bins-0.01)/(max-min);
	
	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		V[((int)((M[j][i]-min)*delta))] += 1.0;
	
	return V;
}

Matrix quarterbinMatrix(const Matrix M)
{
	Matrix N;
	int i, j, m, n;
	
	g_assert(M);

	n = MatrixSize1(M)/2;
        m = MatrixSize2(M)/2;

	N = newMatrix(n, m);

        for(j = 0; j < n; j++) for(i = 0; i < m; i++)
	{
		N[j][i] = 0.25*(M[j*2][i*2] +
				M[j*2+1][i*2] +
				M[j*2][i*2+1] +
				M[j*2+1][i*2+1]);
	}
	
	return N;
}

Matrix rollMatrix(const Matrix M, int delta_n, int delta_m)
{
	Matrix N;
	int i, j, m, n;
		

	g_assert(M);

        n = MatrixSize1(M);
        m = MatrixSize2(M);

        N = newMatrix(n, m);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
                N[j][i] = M[(j + delta_n) % n][(i + delta_m) % m];

	return N;
}

void rollMatrixinplace(Matrix M, int delta_n, int delta_m)
{
	Vector v;
	int rowcount;
	int startrow = 0, row;

        int i, m, n;                

        g_assert(M);

        n = MatrixSize1(M);
        m = MatrixSize2(M);

	delta_n = (n+delta_n) % n;
	delta_m = (m+delta_m) % m;

	v = Matrixrow(M, 0);
	row = (n-delta_n) % n;

	for(rowcount = n-1; rowcount >= 0; rowcount--)
	{
		if(row == startrow)
		{
			for(i = 0; i < m; i++) 
				M[(row + delta_n) % n][(i+delta_m) % m] = v[i];
			startrow = row = (row + 1) % n;
			if(rowcount) 
				for(i = 0; i < m; i++) v[i] = M[row][i];
		}
		else for(i = 0; i < m; i++) 
			M[(row + delta_n) % n][(i + delta_m) % m] = M[row][i];

		row = (row - delta_n + n) % n;
	}

	deleteVector(v);
}

Matrix transposeMatrix(const Matrix M)
{
	Matrix N;
	int i, j, m, n;
	
	g_assert(M);

	n = MatrixSize1(M);
        m = MatrixSize2(M);

	N = newMatrix(m, n);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		N[i][j] = M[j][i];

	return N;
}

void transposeMatrixinplace(Matrix M)
{
	int i, j, n;
	MatrixType t;

	g_assert(M);

	n = MatrixSize1(M);

	if(n != MatrixSize2(M))
	{
		fprintf(stderr, "transposeMatrixinplace : not square\n");
		return;
	}

	for(j = 1; j < n; j++) for(i = 0; i < j; i++)
	{
		t = M[i][j];
		M[i][j] = M[j][i];
		M[j][i] = t;
	}
}


/* reflect matrix such that M[j][i] <--> M[n-j-1][i] */
void reflectMatrix1inplace(Matrix M)
{
	int i, j, m, n, n1, nmax;
	MatrixType t;

	g_assert(M);

	n = MatrixSize1(M);
	m = MatrixSize2(M);

	if(n < 2) return;

	n1 = n-1;
	nmax = n/2;

	for(j = 0; j < nmax; j++)
	{
		for(i = 0; i < m; i++)
		{
			t = M[j][i];
			M[j][i] = M[n1-j][i];
			M[n1-j][i] = t;
		}
	}
}

/* reflect matrix such that M[j][i] <--> M[j][m-i-1] */
void reflectMatrix2inplace(Matrix M)
{
	int i, j, m, n, m1, mmax;
	MatrixType t;

	g_assert(M);

	n = MatrixSize1(M);
	m = MatrixSize2(M);

	if(m < 2) return;

	m1 = m-1;
	mmax = m/2;

	for(j = 0; j < n; j++)
	{
		for(i = 0; i < mmax; i++)
		{
			t = M[j][i];
			M[j][i] = M[j][m1-i];
			M[j][m1-i] = t;
		}
	}
}

void symmetrizeMatrixinplace(Matrix M)
{
	int i, j, n;
	
	g_assert(M);
	g_assert(Matrixissquare(M));

	n = MatrixSize1(M);

	for(i = 1; i < n; i++) for(j = 0; j < i; j++)
		M[i][j] = M[j][i] = 0.5*(M[i][j] + M[j][i]);
}

void antisymmetrizeMatrixinplace(Matrix M)
{
	int i, j, n;
	
	g_assert(M);
	g_assert(Matrixissquare(M));

	n = MatrixSize1(M);

	for(i = 0; i < n; i++) M[i][i] = 0.0;
	for(i = 1; i < n; i++) for(j = 0; j < i; j++)
	{
		M[j][i] = 0.5*(M[j][i] - M[i][j]);
		M[i][j] = -M[j][i];
	}
}

/* Make rotation matrix.  angles in degrees */
Matrix newrotationMatrix(double rx, double ry, double rz)
{
	Matrix M;
	double mx[3][3], my[3][3], mz[3][3], m1[3][3];
	int i, j, k;

	rx = rx*M_PI/180.0;
	ry = ry*M_PI/180.0;
	rz = rz*M_PI/180.0;

	M = newMatrix(3, 3);
	zeroMatrix(M);

	mx[0][0] = 1.0;      mx[0][1] = 0.0;       mx[0][2] = 0.0;
	mx[1][0] = 0.0;      mx[1][1] = cos(rx);   mx[1][2] = sin(rx);
	mx[2][0] = 0.0;      mx[2][1] = -sin(rx);  mx[2][2] = cos(rx);

	my[0][0] = cos(ry);  my[0][1] = 0.0;       my[0][2] = -sin(ry);
	my[1][0] = 0.0;      my[1][1] = 1.0;       my[1][2] = 0.0;
	my[2][0] = sin(ry);  my[2][1] = 0.0;       my[2][2] = cos(ry);

	mz[0][0] = cos(rz);  mz[0][1] = sin(rz);   mz[0][2] = 0.0;
	mz[1][0] = -sin(rz); mz[1][1] = cos(rz);   mz[1][2] = 0.0;
	mz[2][0] = 0.0;      mz[2][1] = 0.0;       mz[2][2] = 1.0;

	for(j = 0; j < 3; j++) for(i = 0; i < 3; i++) m1[j][i] = 0.0;

	for(j = 0; j < 3; j++) for(i = 0; i < 3; i++) for(k = 0; k < 3; k++)
		m1[j][i] += mx[j][k]*my[k][i];
	for(j = 0; j < 3; j++) for(i = 0; i < 3; i++) for(k = 0; k < 3; k++)
		M[j][i] += m1[j][k]*mz[k][i];

	return M;
}

void clipMatrix(Matrix M, double min, double max)
{
	int i, j, m, n;
	
	g_assert(M);

	n = MatrixSize1(M);
        m = MatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
	{	
		if(M[j][i] > max) M[j][i] = max;
		if(M[j][i] < min) M[j][i] = min;
	}
}

/* inverts a lower triangular matrix.  The output is always lower triangular */
Matrix LMatrixinvert(const Matrix L)
{
	Matrix a;
	Vector p;
	int i, j, k, n;
	double sum;
	
	g_assert(L);
	g_assert(Matrixislower(L));
	a = dupMatrix(L);
	n = MatrixSize1(L);
	p = newVector(n);

	for(i = 0; i < n; i++) p[i] = a[i][i];
	
	for(i = 0; i < n; i++)
	{
		a[i][i] = 1.0/p[i];
		for(j = i+1; j < n; j++)
		{
			sum = 0.0;
			for(k = i; k < j; k++) sum -= a[j][k]*a[k][i];
			a[j][i] = sum/p[j];
		}
	}

	deleteVector(p);

	for(i = 1; i < n; i++) for(j = 0; j < i; j++) a[j][i] = 0.0;

	return a;
}

/* returned matrix is lower triangular, ie a[j][i] == 0 for j < i */

Matrix choleskydecompose(const Matrix A)
{
	int N;
	Matrix a;
	Vector p;
	double sum;
	int i, j, k;
	
	g_assert(A);
	g_assert(Matrixissymmetric(A));
	N = MatrixSize1(A);
	g_assert(N > 0);
	g_assert(MatrixSize2(A) == N);

	/* first decompose A a la Num Rec */

	a = dupMatrix(A);
	p = newVector(N);

	for(i = 0; i < N; i++) for(j = i; j < N; j++)
	{
		sum = a[i][j];
		if(i > 0) for(k = i-1; k >= 0; k--) sum -= a[i][k]*a[j][k];
		if(i == j)
		{
			if(sum < 0.0)
			{
				deleteVector(p);
				deleteMatrix(a);
				return 0;
			}
			p[i] = sqrt(sum);
		}
		else a[j][i] = sum/p[i];
	}

	for(i = 0; i < N; i++) a[i][i] = p[i];
	for(i = 1; i < N; i++) for(j = 0; j < i; j++) a[j][i] = 0.0;
	
	deleteVector(p);

	return a;
}

/* same as above, but use long doubles internally */
Matrix choleskydecomposeld(const Matrix A)
{
	int N;
	long double **a;
	long double *p;
	long double sum;
	Matrix M;
	int i, j, k;
	
	g_assert(A);
	g_assert(Matrixissymmetric(A));
	N = MatrixSize1(A);
	g_assert(N > 0);
	g_assert(MatrixSize2(A) == N);

	/* first decompose A a la Num Rec */

	a = g_new(long double *, N);
	for(i = 0; i < N; i++)
	{
		a[i] = g_new(long double, N);
		for(j = 0; j < N; j++) a[i][j] = A[i][j];
	}
	p = g_new(long double, N);

	for(i = 0; i < N; i++) for(j = i; j < N; j++)
	{
		sum = a[i][j];
		if(i > 0) for(k = i-1; k >= 0; k--) sum -= a[i][k]*a[j][k];
		if(i == j)
		{
			if(sum < 0.0)
			{
				g_free(p);
				for(i = 0; i < N; i++) g_free(a[i]);
				g_free(a);
				return 0;
			}
			p[i] = sqrtl(sum);
		}
		else a[j][i] = sum/p[i];
	}

	for(i = 0; i < N; i++) a[i][i] = p[i];
	for(i = 1; i < N; i++) for(j = 0; j < i; j++) a[j][i] = 0.0;

	g_free(p);
	M = newMatrix(N, N);
	for(i = 0; i < N; i++) 
	{
		for(j = 0; j < N; j++) M[i][j] = a[i][j];
		g_free(a[i]);
	}
	g_free(a);
	
	return M;
}

/* Solved A*x = B for symmetric, pos. def. A 
 *
 * call with B = 0 to just factor A into L L^t 
 *
 * returns 0 if decomposition is not possible.  If B != 0, return value
 * is the solution og A x = B.  Else return  non-zero pointer.
 */
Vector cholesky(const Matrix A, const Vector B)
{
	int N;
	Matrix a;
	Vector p, x;
	double sum;
	int i, j, k;
	
	g_assert(A);
	g_assert(Matrixissymmetric(A));
	N = MatrixSize1(A);
	g_assert(N > 0);
	g_assert(MatrixSize2(A) == N);

	/* first decompose A a la Num Rec */

	a = dupMatrix(A);
	p = newVector(N);

	for(i = 0; i < N; i++) for(j = i; j < N; j++)
	{
		sum = a[i][j];
		if(i > 0) for(k = i-1; k >= 0; k--) sum -= a[i][k]*a[j][k];
		if(i == j)
		{
			if(sum <= 0.0)
			{
				deleteVector(p);
				deleteMatrix(a);
				return 0;
			}
			p[i] = sqrt(sum);
		}
		else a[j][i] = sum/p[i];
	}
	
	if(B == 0) return A[0];
	
	g_assert(VectorSize(B) == N);
	
	/* now solve */
	
	x = newVector(N);

	x[0] = B[0]/p[0];
	for(i = 1; i < N; i++)
	{
		sum = B[i];
		for(k = i-1; k >= 0; k--) sum -= a[i][k]*x[k];
		x[i] = sum/p[i];
	}
	
	x[N-1] /= p[N-1];
	for(i = N-2; i >= 0; i--)
	{
		sum = x[i];
		for(k = i+1; k < N; k++) sum -= a[k][i]*x[k];
		x[i] = sum/p[i];
	}

	deleteVector(p);
	deleteMatrix(a);

	return x;
}

int choleskyinvert(Matrix A)
{
	int i, j, k, N;
	Matrix L, Li;
	
	g_assert(Matrixissquare(A));
	N = MatrixSize1(A);
	
	L = choleskydecompose(A);
	if(!L) 
	{
		printf("Inversion failed\n");
		return 0;
	}
	Li = LMatrixinvert(L);
	deleteMatrix(L);

	zeroMatrix(A);
	for(j = 0; j < N; j++) for(k = 0; k <= j; k++) for(i = 0; i <= j; i++)
		A[i][k] += Li[j][i]*Li[j][k];
	
	deleteMatrix(Li);

	return 1;
}

int Matrixisposdef(const Matrix M)
{
	Matrix N;
	int i, j, m;

	g_assert(Matrixissquare(M));
	N = dupMatrix(M);
	m = MatrixSize1(M);
	
	/* symmetrize */
	for(j = 0; j < m; j++) for(i = 0; i < j; i++)
		N[j][i] = N[i][j] = M[j][i] + M[i][j];

	/* if cholesky decomposition is possible, then it is pos def. */
	if(cholesky(N, 0) != 0) return 1;
	else return 0;
}

int Matrixisposdef2(const Matrix M, int niter)
{
	Vector V;
	int i, j, k, m, neg=0;
	double sum;

	g_assert(Matrixissquare(M));
	m = MatrixSize1(M);
	V = newVector(m);

	for(i = 0; i < niter; i++)
	{
		for(j = 0; j < m; j++) V[j] = rand_gauss();
		sum = 0;
		for(j = 0; j < m; j++) for(k = 0; k < m; k++)
			sum += V[j]*M[j][k]*V[k];
		printf("%f\n", sum);
		if(sum < 0) neg++;
	}
	printf("%d of %d were negative\n", neg, niter);
	deleteVector(V);
	if(neg == 0) return 1;
	else return 0;
}

int Matrixisfinite(const Matrix M)
{
	int i, j, m, n;

	n = MatrixSize1(M);
	m = MatrixSize2(M);
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) 
		if(!finite(M[j][i])) return 0;

	return 1;
}

int Matrixislower(const Matrix M)
{
	int i, j, n;

	if(!Matrixissquare(M)) return 0;
	
	n = MatrixSize1(M);

	for(j = 1; j < n; j++) for(i = 0; i < j; i++) 
		if(M[i][j] != 0.0) return 0;

	return 1;
}

int Matrixisupper(const Matrix M)
{
	int i, j, n;

	if(!Matrixissquare(M)) return 0;
	
	n = MatrixSize1(M);

	for(j = 1; j < n; j++) for(i = 0; i < j; i++) 
		if(M[j][i] != 0.0) return 0;

	return 1;
}
	
int Matrixistridiag(const Matrix M)
{
	int i, j, n;

	if(!Matrixissquare(M)) return 0;
	
	n = MatrixSize1(M);
	if(n < 3) return 1;
	
	for(j = 2; j < n; j++) for(i = 0; i < j-1; i++)
		if(M[j][i] != 0.0 || M[i][j] != 0.0) return 0;

	return 1;
}

int Matrixissymmetric(const Matrix M)
{
	int i, j, n;

	if(!Matrixissquare(M)) return 0;
	
	n = MatrixSize1(M);

	for(j = 1; j < n; j++) for(i = 0; i < j; i++) 
		if(M[i][j] != M[j][i]) return 0;

	return 1;
}
