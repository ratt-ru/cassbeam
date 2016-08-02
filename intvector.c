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
#include <glib.h>
#include <stdio.h>
#include "intvector.h"
#include "randdist.h"



/* intVector routines */

intVector newintVector(int n)
{
	intVector v;

	v = g_new(intVectorType, n+1);
	v++;
	intVectorSize(v) = n;

	return v;
}

void zerointVector(intVector v)
{
	int i, m;

	m = intVectorSize(v);
	for(i = 0; i < m; i++) v[i] = 0;
}

void fillintVector(intVector v, int r)
{
	int i, m;

	m = intVectorSize(v);
	for(i = 0; i < m; i++) v[i] = r;
}

void deleteintVector(intVector v)
{
	g_free(v-1);
}

intVector dupintVector(const intVector v)
{
	int i, m;
	intVector v2;

	m = intVectorSize(v);
	v2 = newintVector(m);
	for(i = 0; i < m; i++) v2[i] = v[i];

	return v2;
}

void printintVector(const intVector v)
{
	int i, m;

	m = intVectorSize(v);

	for(i = 0; i < m; i++)
	{
		if(i != 0) printf(" ");
		printf("%d", v[i]);
	}
	printf("\n");
}

void saveintVectorasascii(const intVector v, const char *filename)
{
	FILE *out;
	int i, m;

	g_assert(v);
	out = fopen(filename, "w");
	g_assert(out);
	m = intVectorSize(v);

	for(i = 0; i < m; i++) fprintf(out, "%d %d\n", i, v[i]);

	fclose(out);
}

int intVectorrandomelement(const intVector v)
{
	return v[rand_int(intVectorSize(v))];
}

int intVectormax(const intVector v)
{
	int i, n, m;

	g_assert(v);
	n = intVectorSize(v);
	m = v[0];
	for(i = 1; i < n; i++) if(v[i] > m) m = v[i];

	return m;
}

int intVectormin(const intVector v)
{
	int i, n, m;

	g_assert(v);
	n = intVectorSize(v);
	m = v[0];
	for(i = 1; i < n; i++) if(v[i] < m) m = v[i];

	return m;
}

intMatrix newpaddedintMatrix(int n, int m, int rowpad)
{
	intMatrix M;
	intMatrixTypePointer D;
	int j, w;

	w = m+rowpad;
	D = g_new(intMatrixType, w*n+INTMATRIXVALUESPERBLOCK)
		+ INTMATRIXVALUESPERBLOCK;
	M = g_new(intMatrixTypePointer, ((n>1) ? (n+INTMATRIXSTUBINDICES) :
					         (2+INTMATRIXSTUBINDICES)));
	M += INTMATRIXSTUBINDICES;
	intMatrixStub(M)->datastart = D;
	intMatrixStub(M)->n = n;
	intMatrixStub(M)->m = m;
	intMatrixStub(M)->rowpad = rowpad;
	intMatrixrefcount(M) = 1;
	for(j = 0; j < n; j++) M[j] = D + (j*w);  /* Index array for [][] */
	if(n == 1) M[1] = D + w;

	return M;
}

void deleteintMatrix(intMatrix M)
{
	intMatrixStubType *ms;

	g_assert(M);
	ms = intMatrixStub(M);
	intMatrixrefcount(M)--;
	if(intMatrixrefcount(M) <= 0)
		g_free(ms->datastart-INTMATRIXVALUESPERBLOCK);
	g_free(ms);
}

intMatrix refsubintMatrix(const intMatrix M, int n1, int m1, int n2, int m2)
{
	intMatrix N;
	intMatrixTypePointer D;
	int j, m, n, w;

	g_assert(M);

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	if(n1 < 0) n1 = 0;
	if(n2 < 0) n2 = n-1;
	if(m1 < 0) m1 = 0;
	if(m2 < 0) m2 = m-1;

	if(n2 < n1 || m2 < m1 || m2 >= m || n2 >= n)
	{
		fprintf(stderr, "subintMatrix: limit error\n");
		return 0;
	}

	w = intMatrixstride(M);
	N = g_new(intMatrixTypePointer, ((n>1) ? (n+INTMATRIXSTUBINDICES) :
					         (2+INTMATRIXSTUBINDICES)));
	N += INTMATRIXSTUBINDICES;
	intMatrixStub(N)->datastart = D = intMatrixStub(M)->datastart;
	intMatrixStub(N)->n = n2-n1+1;
	intMatrixStub(N)->m = m2-m1+1;
	intMatrixStub(N)->rowpad = w-(m2-m1+1);
	intMatrixrefcount(M)++;
	for(j = 0; j <= n2-n1; j++) N[j] = M[n1] + m1 + (j*w);
	if(n2 == n1) N[1] = N[0]+w;

	return N;
}

void zerointMatrix(intMatrix M)
{
	int i, j, m, n;
	intMatrixTypePointer Mrow;

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);
	for(j = 0; j < n; j++)
	{
		Mrow = M[j];
		for(i = 0; i < m; i++) Mrow[i] = 0;
	}
}

intMatrix dupintMatrix(const intMatrix M)
{
	int i, j, m, n;
	intMatrix M2;

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);
	M2 = newpaddedintMatrix(n, m, intMatrixrowpad(M));
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) M2[j][i] = M[j][i];

	return M2;
}

intMatrixType intMatrixmax(const intMatrix M)
{
	intMatrixType max;
	int i, j, m, n;

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	max = M[0][0];
	for(i = 0; i < m; i++) for(j = 0; j < n; j++)
		if(M[j][i] > max) max = M[j][i];

	return max;
}

intMatrixType intMatrixmin(const intMatrix M)
{
	intMatrixType min;
	int i, j, m, n;

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	min = M[0][0];
	for(i = 0; i < m; i++) for(j = 0; j < n; j++)
		if(M[j][i] < min) min = M[j][i];

	return min;
}

void intMatrixminmax(const intMatrix M, intMatrixType *min, intMatrixType *max)
{
	int i, j, m, n;

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	*min = *max = M[0][0];

	for(i = 0; i < m; i++) for(j = 0; j < n; j++)
	{
		if(M[j][i] < *min) *min = M[j][i];
		if(M[j][i] > *max) *max = M[j][i];
	}
}

void scaleintMatrix(intMatrix M, double f)
{
	int i, j, m, n;

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++) M[j][i] *= f;
}

void biasintMatrix(intMatrix M, int b)
{
	int i, j, m, n;

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++) M[j][i] += b;
}

void addtointMatrix(intMatrix M, const intMatrix N)
{
	int i, j, m, n;

	g_assert(intMatrixSizeSame(M, N));

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] += N[j][i];
}

intMatrix addintMatrices(const intMatrix M, const intMatrix N)
{
	int i, j, m, n;
	intMatrix Sum;

	g_assert(intMatrixSizeSame(M, N));

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	Sum = newintMatrix(n, m);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		Sum[j][i] = M[j][i] + N[j][i];

	return Sum;
}

void copytointMatrix(intMatrix M, const intMatrix N)
{
	int i, j, m, n;

	g_assert(intMatrixSizeSame(M, N));

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] = N[j][i];
}

void applyfunctointMatrix(intMatrix M, int (*func)(int x))
{
	int i, j, m, n;

	g_assert(M);

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] = func(M[j][i]);
}

intMatrix intMatrixproduct(intVector v2, intVector v1)
{
	intMatrix M;
	int i, j, m, n;

	g_assert(v2);
	g_assert(v1);

	M = newintMatrix(intVectorSize(v2), intVectorSize(v1));

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	for(j = 0; j < n; j++) for(i = 0; i < m; i++)
		M[j][i] = v2[j]*v1[i];

	return M;
}

intMatrix dupsubintMatrix(const intMatrix M, int n1, int m1, int n2, int m2)
{
	intMatrix N;
	int i, j, m, n; 

	g_assert(M);

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	if(n1 < 0) n1 = 0;
	if(n2 < 0) n2 = n-1;
	if(m1 < 0) m1 = 0;
	if(m2 < 0) m2 = m-1;

	if(n2 < n1 || m2 < m1 || m2 >= m || n2 >= n)
	{
		fprintf(stderr, "dupsubintMatrix: limit error\n");
		return 0;
	}

	N = newintMatrix(n2-n1+1, m2-m1+1);

	for(j = n1; j <= n2; j++) for(i = m1; i <= m2; i++)
		N[j-n1][i-m1] = M[j][i];

	return N;
}

void copysubintMatrix(intMatrix N, const intMatrix M, int n1, int m1, 
        int n2, int m2, int n3, int m3)
{
	int i, j, m, n;

	g_assert(M);
	g_assert(N);

	n = intMatrixSize1(M);
	m = intMatrixSize2(M);

	if(n1 < 0) n1 = 0;
	if(n2 < 0) n2 = n-1;
	if(m1 < 0) m1 = 0;
	if(m2 < 0) m2 = m-1;

	if(n2 < n1 || m2 < m1 || m2 >= m || n2 >= n ||
		n3 >= intMatrixSize1(N) || m3 >= intMatrixSize2(N))
	{
		fprintf(stderr, "copysubMatrix: limit error\n");
		return;
	}

	if(n3 < 0) n3 = intMatrixSize1(N)-(n2-n1);
	if(m3 < 0) m3 = intMatrixSize2(N)-(m2-m1);

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

	if(n3+n2-n1 >= intMatrixSize1(N))
		n2 = n1-n3+intMatrixSize1(N)-1;
	if(m3+m2-m1 >= intMatrixSize2(N))
		m2 = m1-m3+intMatrixSize2(N)-1;

	for(j = n1; j <= n2; j++) for(i = m1; i <= m2; i++)
		N[j+n1-n3][i+m1-m3] = M[j][i];
}

	
