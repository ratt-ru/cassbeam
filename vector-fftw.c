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
#include "vector-fftw.h"

inline double sinc(double x)
{
	if(x == 0.0) return 1.0;
	else return sin(x)/x;
}

int Matrixcomplexgrid(Matrix M, double u, double v, double re, double im, int p)
{
	int m, n;
	int i, j;
	int a, b, r, s;
	double du, dv, f;

	n = MatrixSize1(M);
	m = MatrixSize2(M)/2;
	
	while(u < 0.0) u += m;
	while(u >= m) u -= m;

	while(v < 0.0) v += n;
	while(v >= n) v -= n;

	a = (int)(u+0.5);
	b = (int)(v+0.5);

	du = u-a;
	dv = v-b;

	for(j = -p; j <= p; j++) for(i = -p; i <= p; i++)
	{
		f = sinc(M_PI*(i-du))*sinc(M_PI*(j-dv));
		r = (b+j+n)%n;
		s = (a+i+m)%m;
		M[r][2*s] += f*re;
		M[r][2*s+1] += f*im;
	}

	return 1;
}

Vector VectorFFT(Vector V, fftw_direction dir)
{
	Vector out;
	fftw_plan p;

	p = fftw_create_plan(Vectorlength(V)/2, dir, FFTW_ESTIMATE);
	out = newVector(Vectorlength(V));
	fftw_one(p, (fftw_complex *)V, (fftw_complex *)out);
	fftw_destroy_plan(p);

	return out;
}

/* output vector is in FFTW half-complex format */
Vector VectorRFFT(Vector V, fftw_direction dir)
{
	Vector out;
	rfftw_plan p;
	
	p = rfftw_create_plan(Vectorlength(V), dir, FFTW_ESTIMATE);
	out = newVector(Vectorlength(V));
	rfftw_one(p, (fftw_real *)V, (fftw_real *)out);
	rfftw_destroy_plan(p);

	return out;
}

void smoothVector(Vector V, double width, int type)
{
	int pad, trans, m, n, i;
	double fac, re, im;
	Vector K, ftK, W, ftW;
	
	if(width <= 0.0) return;

	g_assert(V);
	n = VectorSize(V);
	g_assert(n > 2.0*width);

	if(type == SMOOTH_GAUSS) pad = width*3;
	else pad = width;
	if(n%2 == 1) trans = width+1;
	else trans = width;

	m = n + pad*2 + trans;
	i = m%256;
	if(i) 
	{
		trans += (256-i);
		m += (256-i);
	}
	
	W = newVector(m);
	for(i = 0; i < n; i++) W[i] = V[i];
	for(i = n; i < n+pad; i++) W[i] = V[n-1];
	for(i = n+pad+trans; i < m; i++) W[i] = V[0];
	for(i = 0; i < trans; i++) W[n+pad+i] = V[n-1] 
		+ (V[0]-V[n-1])*(float)i/(float)trans;

	ftW = VectorRFFT(W, FFTW_REAL_TO_COMPLEX);
	deleteVector(W);

	K = newVector(m);
	zeroVector(K);
	K[0] = 1.0;
	switch(type)
	{
		case SMOOTH_GAUSS:
			for(i = 1; i <= m/2; i++)
				K[i] = K[m-i] = exp(-0.5*i*i/(width*width));
			break;
		case SMOOTH_COS:
			for(i = 1; i < width; i++)
				K[i] = K[m-i] = cos(0.5*M_PI*i/width);
			break;
		case SMOOTH_COS1:
			for(i = 1; i < width; i++)
				K[i] = K[m-i] = 0.5+0.5*cos(M_PI*i/width);
			break;
		case SMOOTH_BOX:
			for(i = 1; i < width; i++)
				K[i] = K[m-i] = 1.0;
			break;
		case SMOOTH_TRIANGLE:
			for(i = 1; i < width; i++)
				K[i] = K[m-i] = 1.0 - i/width;
			break;
	}

	scaleVector(K, 1.0/Vectorsum(K));

	ftK = VectorRFFT(K, FFTW_REAL_TO_COMPLEX);
	deleteVector(K);

	ftW[0]   *= ftK[0];
	ftW[m/2] *= ftK[m/2];
	for(i = 1; i < m/2; i++)
	{
		re       = ftW[i];
		im       = ftW[m-i];
		ftW[i]   = re*ftK[i] - im*ftK[m-i];
		ftW[m-i] = im*ftK[i] + re*ftK[m-i];
	}
	
	fac = 1.0/m;
	
	deleteVector(ftK);
	
	W = VectorRFFT(ftW, FFTW_COMPLEX_TO_REAL);
	deleteVector(ftW);

	for(i = 0; i < n; i++) V[i] = W[i]*fac;

	deleteVector(W);
}
	
/* This will not work on a padded Matrix */
void MatrixFFT(Matrix M, fftw_direction dir)
{
	fftwnd_plan p;

	g_assert(M);
	g_assert(Matrixrowpad(M) == 0);

	p = fftw2d_create_plan(MatrixSize1(M), MatrixSize2(M)/2, dir, 
		FFTW_ESTIMATE | FFTW_IN_PLACE);
	fftwnd_one(p, (fftw_complex *)(*M), 0);
	fftwnd_destroy_plan(p);
}

/* This will work on a padded or unpadded Matrix */
void MatrixFFTrows(Matrix M, fftw_direction dir)
{
	fftw_plan p;

	g_assert(M);

	p = fftw_create_plan(MatrixSize2(M)/2, dir, 
		FFTW_ESTIMATE | FFTW_IN_PLACE);
	fftw(p, MatrixSize1(M), (fftw_complex *)(*M), 1, Matrixstride(M)/2, 
		0, 0, 0);
	fftw_destroy_plan(p);
}

void MatrixFFTcolumns(Matrix M, fftw_direction dir)
{
	fftw_plan p;

	g_assert(M);

	p = fftw_create_plan(MatrixSize1(M), dir,
		FFTW_ESTIMATE | FFTW_IN_PLACE);
	fftw(p, MatrixSize2(M)/2, (fftw_complex *)(*M), Matrixstride(M)/2, 
		1, 0, 0, 0);
	fftw_destroy_plan(p);
}

void MatrixRFFT(Matrix M, fftw_direction dir)
{
	rfftwnd_plan p;

	g_assert(M);

	if(Matrixrowpad(M) != 2 - (MatrixSize2(M) % 2))
	{
		fprintf(stderr, "MatrixRFFT: padding is: %d, should be %d\n",
			Matrixrowpad(M), 2 - (MatrixSize2(M) % 2));
		g_assert(0);
	}

	p = rfftw2d_create_plan(MatrixSize1(M), MatrixSize2(M), dir,
		FFTW_ESTIMATE | FFTW_IN_PLACE);
	if(dir == FFTW_REAL_TO_COMPLEX)
		rfftwnd_one_real_to_complex(p, (fftw_real *)(*M), 0);
	else
		rfftwnd_one_complex_to_real(p, (fftw_complex *)(*M), 0);

	rfftwnd_destroy_plan(p);
}

void MatrixRFFTrows(Matrix M, fftw_direction dir)
{
	rfftw_plan p;

	g_assert(M);
	
	p = rfftw_create_plan(MatrixSize2(M), dir, 
		FFTW_ESTIMATE | FFTW_IN_PLACE);
	
	rfftw(p, MatrixSize1(M), (fftw_real *)(*M), 1, Matrixstride(M),
		0, 0, 0);

	rfftw_destroy_plan(p);
}

void Matrixconvolve(Matrix M, const Matrix N)
{
	int m, n, i, j;
	Matrix R, T;

	g_assert(M);
	g_assert(N);

	n = MatrixSize1(M);
	m = MatrixSize2(M);


	if(MatrixSize1(N) > n || MatrixSize2(N) > m)
	{
		fprintf(stderr, "Matrixconvolve: Size Problem\n");
		return;
	}

	R = newpaddedMatrix(n, m, 2);
	copytoMatrix(R, M);
	
	T = newpaddedMatrix(n, m, 2);
	zeroMatrix(T);
	
	for(j = 0; j < MatrixSize1(N); j++) 
		for(i = 0; i < MatrixSize2(N); i++)
			T[j][i] = N[j][i];

	MatrixRFFT(T, FFTW_REAL_TO_COMPLEX);
	MatrixRFFT(R, FFTW_REAL_TO_COMPLEX);
	Matrixcomplexmultiply(R, T, 2);
	MatrixRFFT(R, FFTW_COMPLEX_TO_REAL);

	copytoMatrix(M, R);
	deleteMatrix(T);
	deleteMatrix(R);
}

Matrix Matrixautocorrelate(const Matrix M)
{
	int m, n;
	Matrix R;

	g_assert(M);

	n = MatrixSize1(M);
	m = MatrixSize2(M);

	R = newpaddedMatrix(n, m, 2);
	copytoMatrix(R, M);
	
	MatrixRFFT(R, FFTW_REAL_TO_COMPLEX);
	Matrixcomplexconjugatemultiply(R, R, 2);
	MatrixRFFT(R, FFTW_COMPLEX_TO_REAL);

	return R;
}

Matrix Matrixcrosscorrelate(const Matrix M, const Matrix N)
{
	int m, n;
	Matrix R, T;

	g_assert(M);
	g_assert(N);

	n = MatrixSize1(M);
	m = MatrixSize2(M);


	if(MatrixSize1(N) != n || MatrixSize2(N) != m)
	{
		fprintf(stderr, "Matrixcrosscorrelate: Size Problem\n");
		return 0;
	}

	R = newpaddedMatrix(n, m, 2);
	copytoMatrix(R, M);
	
	T = newpaddedMatrix(n, m, 2);
	copytoMatrix(T, N);

	MatrixRFFT(T, FFTW_REAL_TO_COMPLEX);
	MatrixRFFT(R, FFTW_REAL_TO_COMPLEX);
	Matrixcomplexconjugatemultiply(R, T, 2);
	MatrixRFFT(R, FFTW_COMPLEX_TO_REAL);

	deleteMatrix(T);

	return R;
}

void Matrixdeconvolve(Matrix M, const Matrix N)
{
	int m, n, i, j;
	Matrix R, T;

	g_assert(M);
	g_assert(N);

	n = MatrixSize1(M);
	m = MatrixSize2(M);


	if(MatrixSize1(N) > n || MatrixSize2(N) > m)
	{
		fprintf(stderr, "Matrixconvolve: SizeProblem\n");
		return;
	}

	R = newpaddedMatrix(n, m, 2);
	copytoMatrix(R, M);
	
	T = newpaddedMatrix(n, m, 2);
	zeroMatrix(T);
	
	for(j = 0; j < MatrixSize1(N); j++) 
		for(i = 0; i < MatrixSize2(N); i++)
			T[j][i] = N[j][i];

	MatrixRFFT(T, FFTW_REAL_TO_COMPLEX);
	MatrixRFFT(R, FFTW_REAL_TO_COMPLEX);
	Matrixcomplexdivide(R, T, 2);
	MatrixRFFT(R, FFTW_COMPLEX_TO_REAL);

	copytoMatrix(M, R);
	deleteMatrix(T);
	deleteMatrix(R);
}
		
