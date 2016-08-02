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
#include <stdarg.h>
#include <glib.h>
#include <stdio.h>
#include <math.h>
#include "vecarray.h"

VecArray newVecArray(int n)
{
	VecArray a;
	int i;

	a = g_new(VectorTypePointer, n+1);
	a++;
	VecArraySize(a) = n;

	for(i = 0; i < n; i++) a[i] = 0;

	return a;
}

VecArray newpopulatedVecArray(int n, int m)
{
	VecArray a;
	int i;

	a = g_new(VectorTypePointer, n+1);
	a++;
	VecArraySize(a) = n;

	for(i = 0; i < n; i++) a[i] = newVector(m);

	return a;
}

void zeroVecArray(VecArray a)
{
	int i, n;

	g_assert(a);

	n = VecArraySize(a);
	for(i = 0; i < n; i++) if(a[i]) zeroVector(a[i]);
}

void deleteVecArray(VecArray a)
{
	g_assert(a);

	VecArraySize(a) = 0;

	g_free(a-1);
}

void deleteVecArrayandVectors(VecArray a)
{
	int i, n;
	
	g_assert(a);

	n = VecArraySize(a);
	for(i = 0; i < n; i++) if(a[i]) deleteVector(a[i]);
	VecArraySize(a) = 0;

	g_free(a-1);
}

VecArray dupVecArray(const VecArray a)
{
	VecArray b;
	int i, n;

	g_assert(a);

	n = VecArraySize(a);
	b = newVecArray(n);
	for(i = 0; i < n; i++) b[i] = a[i];

	return b;
}

VecArray dupVecArrayandVectors(const VecArray a)
{
	VecArray b;
	int i, n;

	g_assert(a);

	n = VecArraySize(a);
	b = newVecArray(n);
	for(i = 0; i < n; i++) if(a[i]) b[i] = dupVector(a[i]);

	return b;
}

void copytoVecArray(VecArray dest, const VecArray src)
{
	int i, n;

	g_assert(dest);
	g_assert(src);

	n = VecArraySize(src);
	
	g_assert(n == VecArraySize(dest));
	
	for(i = 0; i < n; i++) copytoVector(dest[i], src[i]);
}

VecArray VecArraysubset(const VecArray a, int col, int ncol)
{
	VecArray b;
	int i;

	g_assert(a);
	g_assert(col > 0);
	g_assert(ncol > 0);
	g_assert(col+ncol <= VecArraySize(a));

	b = newVecArray(ncol);

	for(i = 0; i < ncol; i++) b[i] = a[i+col];

	return b;
}

VecArray dupVecArraysubset(const VecArray a, int col, int ncol)
{
	VecArray b;
	int i;

	g_assert(a);
	g_assert(col >= 0);
	g_assert(ncol > 0);
	g_assert(col+ncol <= VecArraySize(a));

	b = newVecArray(ncol);

	for(i = 0; i < ncol; i++) if(a[i]) b[i] = dupVector(a[i+col]);

	return b;
}

VecArray dupVecArrayrowsubset(const VecArray a, int row, int nrow)
{
	VecArray b;
	int N, i, j;

	g_assert(a);
	N = VecArraySize(a);
	g_assert(row >= 0);
	g_assert(nrow > 0);
	g_assert(row+nrow <= VecArrayVectorSize(a));

	b = newpopulatedVecArray(N, nrow);
	for(j = 0; j < nrow; j++) for(i = 0; i < N; i++)
		b[i][j] = a[i][j+row];

	return b;
}

int VecArrayVectorSize(const VecArray a)
{
	int i, n, s=-1;

	g_assert(a);

	n = VecArraySize(a);
	for(i = 0; i < n; i++) if(a[i])
	{
		if(s < 0) s = VectorSize(a[i]);
		else if(s != VectorSize(a[i])) return -1;
	}

	return s;
}

int VecArrayVectorcolumns(const VecArray a)
{
	int i, n, c=0;

	g_assert(a);

	n = VecArraySize(a);
	for(i = 0; i < n; i++) if(a[i]) c++;
	
	return c;
}

VecArray VecArraydeletecolumn(VecArray a, int col)
{
	VecArray b;
	int i, n;

	g_assert(a);
	
	n = VecArraySize(a);
	g_assert(col >= 0 && col < n);
	
	if(n <= 1)
	{
		deleteVecArray(a);
		return 0;
	}

	b = newVecArray(n-1);
	
	if(col) for(i = 0; i < col; i++) b[i] = a[i];
	if(col < n-1) for(i = col+1; i < n; i++) b[i-1] = a[i];

	deleteVecArray(a);

	return b;
}

VecArray VecArrayinsertcolumn(VecArray a, int col)
{
	VecArray b;
	int i, n;
	
	g_assert(col >= 0);
	
	if(a == 0) return newVecArray(1);

	n = VecArraySize(a);
	g_assert(col <= n);

	b = newVecArray(n+1);
	
	if(col > 0) for(i = 0; i < col; i++) b[i] = a[i];
	if(col < n) for(i = col; i < n; i++) b[i+1] = a[i];

	deleteVecArray(a);

	return b;
}

VecArray VecArrayappendcolumn(VecArray a)
{
	VecArray b;
	int i, n;

	if(!a) return newVecArray(1);

	n = VecArraySize(a);
	b = newVecArray(n+1);
	for(i = 0; i < n; i++) b[i] = a[i];

	deleteVecArray(a);

	return b;
}

VecArray VecArrayinsertVector(VecArray a, Vector v, int col)
{
	VecArray b;

	b = VecArrayinsertcolumn(a, col);
	b[col] = v;

	return b;
}

VecArray VecArrayappendVector(VecArray a, Vector v)
{
	VecArray b;
	int i, n;

	if(!a)
	{
		b = newVecArray(1);
		b[0] = v;
		return b;
	}

	n = VecArraySize(a);
	b = newVecArray(n+1);
	for(i = 0; i < n; i++) b[i] = a[i];

	deleteVecArray(a);

	b[n] = v;

	return b;
}

VecArray VecArrayappendVectors(VecArray a, ...) 
{
	va_list ap;
	Vector v;

	va_start(ap, a);
	while((v = va_arg(ap, Vector)))
		a = VecArrayappendVector(a, v);
	va_end(ap);

	return a;
}

/*
VecArray newVecArrayfromVectors(...)
{
	va_list ap;
	VecArray a = 0;
	Vector v;

	va_start(ap, a);
	while((v = va_arg(ap, Vector)))
		a = VecArrayappendVector(a, v);
	va_end(ap);

	return a;
}
*/

void VecArraysetcolumn(VecArray a, Vector v, int col)
{
	g_assert(a);
	g_assert(col >= 0 && col < VecArraySize(a));
	
	a[col] = v;
}

Matrix newMatrixfromVecArray(const VecArray a)
{
	Matrix M;
	int m, n, i, j;

	g_assert(a);

	m = VecArrayVectorSize(a);
	if(m <= 0) return 0;
	n = VecArraySize(a);
	
	M = newMatrix(n, m);
	for(j = 0; j < n; j++) for(i = 0; i < m; i++) M[j][i] = a[j][i];

	return M;
}

VecArray newVecArrayfromMatrix(const Matrix M)
{
	int i, j, m, n;
	VecArray a;

	g_assert(M);
	
	n = MatrixSize1(M);
	m = MatrixSize2(M);

	a = newVecArray(n);
	for(j = 0; j < n; j++) 
	{
		a[j] = newVector(m);
		for(i = 0; i < m; i++) a[j][i] = M[j][i];
	}

	return a;
}

VecArray VecArrayconcat(VecArray a, VecArray b)
{
	VecArray c;
	int i, na, nb;

	g_assert(a);
	g_assert(b);
	
	na = VecArraySize(a);
	nb = VecArraySize(b);

	c = newVecArray(na+nb);
	for(i = 0; i < na; i++) c[i] = a[i];
	for(i = 0; i < nb; i++) c[i+na] = b[i];

	deleteVecArray(a);
	deleteVecArray(b);

	return c;
}

VecArray VecArrayappendVecArray(VecArray a, const VecArray b)
{
	VecArray c;
	int i, na, nb;

	g_assert(a);
	g_assert(b);
	
	na = VecArraySize(a);
	nb = VecArraySize(b);

	c = newVecArray(na+nb);
	for(i = 0; i < na; i++) c[i] = a[i];
	for(i = 0; i < nb; i++) c[i+na] = b[i];

	deleteVecArray(a);

	return c;
}

VecArray VecArrayappendVecArrays(VecArray a, ...) 
{
	va_list ap;
	VecArray b;

	va_start(ap, a);
	while((b = va_arg(ap, VecArray)))
		a = VecArrayappendVecArray(a, b);
	va_end(ap);

	return a;
}

VecArray reduceVecArraybymask(VecArray a, const intVector mask)
{
	VecArray b;
	int i, m, j, n, c=0;
	
	g_assert(a);
	g_assert(mask);
	
	m = VecArrayVectorSize(a);
	g_assert(m == VectorSize(mask));

	n = VecArraySize(a);
	
	for(i = 0; i < m; i++) if(mask[i] > 0) c++;
	
	b = newpopulatedVecArray(n, c);

	c = 0;

	for(i = 0; i < m; i++) if(mask[i])
	{
		for(j = 0; j < n; j++) b[j][c] = a[j][i];
		c++;
	}

	deleteVecArrayandVectors(a);

	return b;
}

VecArray VecArrayfromfile(const char *filename, int mincols)
{
	VecArray a;
	Vector v;
	FILE *in;
	int nc = 0, nr = 0, n, i, row;
	char line[1024];
	
	in = fopen(filename, "r");
	if(!in) return 0;

	for(;;)
	{
		fgets(line, 1023, in);
		if(feof(in)) break;
		/* look for "standard" comment characters */
		if(line[0] == '#' || line[0] > 57 || line[0] == '*' 
				  || line[0] == '!' || line[0] == ';') continue;
		n = stringtoVector(line, 0);
		if(n >= mincols)
		{
			if(nc == 0) nc = n;
			if(n >= nc) nr++;
		}
	}
	fclose(in);

	a = newpopulatedVecArray(nc, nr);
	v = newVector(nc);

	in = fopen(filename, "r");
	for(row = 0; row < nr;)
	{
		fgets(line, 1023, in);
		if(feof(in)) break;
		if(line[0] == '#' || line[0] > 57 || line[0] == '*') continue;
		n = stringtoVector(line, v);
		if(n >= nc)
		{
			if(n >= nc)
			{
				for(i = 0; i < nc; i++) a[i][row] = v[i];
				row++;
			}
		}
	}
	fclose(in);

	deleteVector(v);

	return a;
}

int VecArraytofile(const VecArray a, const char *filename, 
	const char *format)
{
	int m, n, i, j, needspace = 0;
	const char *f;
	FILE *out;

	g_assert(a);

	if(format == 0) f = "%10.5f";
	else f = format;

	m = VecArrayVectorSize(a);
	if(m <= 0) 
	{
		fprintf(stderr, "VecArraycolumnstofile : "
			"VecArray is too complicated\n");
		return -1;
	}
	n = VecArraySize(a);

	out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "VecArraycolumnstofile : "
			"Cannot open file %s for write\n", filename);
		return -1;
	}

	for(i = 0; i < m; i++) 
	{
		needspace = 0;
		for(j = 0; j < n; j++) if(a[j]) 
		{
			if(needspace++) fprintf(out, " ");
			fprintf(out, f, a[j][i]);
		}
		fprintf(out, "\n");
	}

	fclose(out);

	return needspace;
}

void VecArrayresize(const VecArray a, int newsize)
{
	int i, j, m, nr, nc;
	Vector V;
	
	g_assert(a);
	
	nc = VecArraySize(a);
	nr = VecArrayVectorSize(a);

	g_assert(nr);
	g_assert(nc);

	if(newsize < nr) m = newsize;
	else m = nr;

	for(i = 0; i < nc; i++)
	{
		V = newVector(newsize);
		for(j = 0; j < m; j++) V[j] = a[i][j];
		if(newsize > nr)
			for(j = nr; j < newsize; j++) a[i][j] = 0.0;
		deleteVector(a[i]);
		a[i] = V;
	}
}

/* Arithmetic operations */

void scaleVecArray(VecArray a, double f)
{
	int N, j;

	g_assert(a);

	N = VecArraySize(a);
	for(j = 0; j < N; j++) if(a[j]) scaleVector(a[j], f);
}

void addtoVecArray(VecArray a, const VecArray b)
{
	int i, M;
	
	g_assert(a);
	g_assert(b);
	
	M = VecArraySize(a);
	g_assert(M == VecArraySize(a));

	for(i = 0; i < M; i++) addtoVector(a[i], b[i]);
}

Vector VecArrayVectoraverage(const VecArray a)
{
	int N, M, j;
	Vector v;

	g_assert(a);

	N = VecArrayVectorSize(a);
	if(N <= 0) return 0;

	M = VecArraySize(a);

	v = newVector(N);
	zeroVector(v);
	for(j = 0; j < M; j++) addtoVector(v, a[j]);
	scaleVector(v, 1.0/M);

	return v;
}


/* reduce symmetric matrix to tridiagonal form 
 *
 * From numerical recipes tred2
 */
VecArray tridiagfromMatrix(Matrix a)
{
	int n, i, j, k, l;
	double scale, hh, h, g, f;
	VecArray out;
	Vector d, e;

	g_assert(Matrixissquare(a));
	
	n = MatrixSize1(a);

	out = newpopulatedVecArray(2, n);
	d = out[0];
	e = out[1];

	for(i = n-1; i > 0; i--)
	{
		l = i-1;
		h = scale = 0.0;
		if(l > 0)
		{
			for(k = 0; k <= l; k++)
				scale += fabs(a[i][k]);
			if(scale == 0.0)
				e[i] = a[i][i];
			else
			{
				for(k = 0; k <= l; k++)
				{
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f = a[i][i];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale*g;
				h -= f*g;
				a[i][l] = f - g;
				f = 0.0;
				for(j = 0; j <= l; j++)
				{
					a[j][i] = a[i][j]/h;
					g = 0.0;
					for(k = 0; k <= j; k++)
						g += a[j][k]*a[i][k];
					for(k = j+1; k <= l; k++)
						g += a[k][j]*a[i][k];
					e[j] = g/h;
					f += e[j]*a[i][j];
				}
				hh = 0.5*f/h;
				for(j = 0; j <= l; j++)
				{
					f=a[i][j];
					e[j] = g = e[j] - hh*f;
					for(k = 0; k <= j; k++)
						a[j][k] -= (f*e[k] + g*a[i][k]);
				}
			}
		}
		else e[i] = a[i][l];
		d[i] = h;
	}
	d[0] = e[0] = 0.0;

	for(i = 0; i < n; i++)
	{
		l = i - 1;
		if(i)
		{
			for(j = 0; j <= l; j++)
			{
				g = 0.0;
				for(k = 0; k <= l; k++)
					g += a[i][k]*a[k][j];
				for(k = 1; k <= l; k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for(j = 0; j <= l;j++) a[j][i] = a[i][j] = 0.0;
	}
	
	return out;
}

void multVectorbytridiag(Vector v, const VecArray a)
{
	int i, n;
	Vector w, d, e;

	g_assert(VecArraySize(a) == 2);
	n = VectorSize(v);
	g_assert(VecArrayVectorSize(a) == n);
	d = a[0];
	e = a[1];
	
	w = dupVector(v);

	v[0] = w[0]*d[0] + w[1]*e[1];
	v[n-1] = w[n-1]*d[n-1] + w[n-2]*e[n-1];
	for(i = 1; i < n-1; i++)
		v[i] = w[i]*d[i] + w[i-1]*e[i-1] + w[i+1]*e[i];
	
	deleteVector(w);
}
