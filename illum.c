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
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include "vector.h"
#include "image-vector.h"
#include "vector-fftw.h"
#include "keyvalue.h"
#include "illum.h"

Illum *newIllum(int N)
{
	Illum *I;
	int i;

	I = g_new(Illum, 1);

	I->A = newMatrix(N, N);		
	I->P = newMatrix(N, N);		
	I->B = newMatrix(N, N);		
	zeroMatrix(I->A);
	zeroMatrix(I->P);
	zeroMatrix(I->B);
	I->E = g_new(Matrix, 12);
	for(i = 0; i < 12; i++) 
	{
		I->E[i] = newMatrix(N, N);
		zeroMatrix(I->E[i]);
	}
	I->groundfrac = 0.0;
	/* Fill in reasonable temperatures.  should be updated */
	I->Tsky = 3.0;
	I->Tground = 290.0;
	I->Trec = 20.0;
	I->pointx = I->pointy = 0.0;
	I->spilleff = 0.0;
	I->prispilleff = 0.0;
	I->subspilleff = 0.0;
	I->blockeff = 0.0;
	I->surfeff = 0.0;
	I->illumeff = 0.0;
	I->ampeff = 0.0;
	I->phaseeff = 0.0;
	I->efficiency = 0.0;
	I->diffeff = 0.0;
	I->misceff = 1.0;
	I->aeff = 0.0;
	I->fwhm_x = I->fwhm_y = 0.0;
	I->peaksidelobe = 0.0;
	I->lambda = 0.0;
	I->Xangle = 0.0;
	I->prefix = g_strdup("illum");
	I->pixelsperbeam = 32;
	I->beampixelscale = 0.0;

	I->savejonesmatrices = 1;
	I->savestokesimages = 1;
	I->saveapertureimages = 1;
	I->saveparams = 1;

	return I;
}

void setIllumfromKeyValue(Illum *I, const struct KeyValue *kv)
{
	int i;
	double v;
	const char *str;

	v = getKeyValuedouble(kv, "misceff");
	if(v != KV_FLOATERR) I->misceff = v;
	
	v = getKeyValuedouble(kv, "diffeff");
	if(v != KV_FLOATERR) I->diffeff = v;

	str = getKeyValuestring(kv, "out");
	if(str) 
	{
		if(I->prefix) g_free(I->prefix);
		I->prefix = g_strdup(str);
	}
	else if(!I->prefix) I->prefix = g_strdup("illum");

	v = getKeyValuedouble(kv, "pixelsperbeam");
	if(v != KV_FLOATERR) I->pixelsperbeam = v;

	str = getKeyValuestring(kv, "compute");
	if(str)
	{
		if(strcmp(str, "all") == 0 || strcmp(str, "ALL") == 0)
		{
			I->savejonesmatrices = 1;
			I->savestokesimages = 1;
			I->saveapertureimages = 1;
			I->saveparams = 1;
		}
		else if(strcmp(str, "none") == 0 || strcmp(str, "NONE") == 0)
		{
			I->savejonesmatrices = 0;
			I->savestokesimages = 0;
			I->saveapertureimages = 0;
			I->saveparams = 0;
		}
		else
		{
			I->savejonesmatrices = 0;
			I->savestokesimages = 0;
			I->saveapertureimages = 0;
			I->saveparams = 0;
			for(i = 0; str[i]; i++)
			{
				if(str[i] == 'j' || str[i] == 'J')
					I->savejonesmatrices = 1;
				if(str[i] == 's' || str[i] == 'S')
					I->savestokesimages = 1;
				if(str[i] == 'a' || str[i] == 'A')
					I->saveapertureimages = 1;
				if(str[i] == 'p' || str[i] == 'P')
					I->saveparams = 1;
			}
		}
		
	}
}

void deleteIllum(Illum *I)
{
	int i;
	if(!I) return;

	if(I->A) deleteMatrix(I->A);
	if(I->P) deleteMatrix(I->P);
	if(I->B) deleteMatrix(I->B);
	if(I->E)
	{
		for(i = 0; i < 12; i++) if(I->E[i]) 
			deleteMatrix(I->E[i]);
		g_free(I->E);
	}
}

void calcIllumparams(Illum *I)
{
	double er, ei, e, e2;
	/* sums over unblocked portion of aperture */
	double aEr=0.0, aEi=0.0, aE=0.0, aE2=0.0, aA=0.0;
	/* sums over blocked portion of aperture */
	double bEr=0.0, bEi=0.0;
	double a, b, abEr, abEi;
	int i, j, m, n;

	g_assert(MatrixSizeSame(I->A, I->P));
	
	n = MatrixSize1(I->A);
	m = MatrixSize2(I->A);
	
	for(i = 0; i < m; i++) for(j = 0; j < n; j++) if(I->A[j][i] > 0.0)
	{
		b = I->B[j][i];
		a = 1.0-b;
		er = I->A[j][i]*cos(I->P[j][i]);
		ei = I->A[j][i]*sin(I->P[j][i]);
		e  = I->A[j][i];
		e2 = I->A[j][i]*I->A[j][i];
		
		aEr += a*er;
		aEi += a*ei;
		aE  += a*e;
		aE2 += a*e2;
		aA  += a;
		bEr += b*er;
		bEi += b*ei;
	}
	aEr *= I->dA;
	aEi *= I->dA;
	aE  *= I->dA;
	aE2 *= I->dA;
	aA  *= I->dA;
	bEr *= I->dA;
	bEi *= I->dA;

	abEr = aEr + bEr;
	abEi = aEi + bEi;
	
	I->ampeff   = (aE*aE)/(aA*aE2);
	I->phaseeff = (aEr*aEr + aEi*aEi)/(aE*aE);
	I->blockeff = (aEr*aEr + aEi*aEi)/(abEr*abEr + abEi*abEi);
	I->illumeff = I->ampeff*I->phaseeff;
	I->efficiency = I->spilleff * I->blockeff * I->surfeff * I->illumeff
		* I->diffeff * I->misceff;
	I->gain = I->efficiency*4.0*M_PI*I->a0/(I->lambda*I->lambda);
	I->aeff = I->lambda*I->lambda*I->gain/(4.0*M_PI);

	I->pTground = I->Tground*I->groundfrac;
	I->pTsky = I->Tsky*(1.0-I->groundfrac);
	I->pTrec = I->Trec;
	I->Tsys = I->pTground + I->pTsky + I->pTrec;
}

/* returns the matrix that rotates the XY Efield into RL.  Element order is:
 * Input  : ReX, ImX, ReY, ImY
 * Output : ReR, ImR, ReL, ImL
 */
static Matrix calcSmatrix(const Illum *I)
{
	Matrix S;
	double c, s;
	
	S = newMatrix(4, 4);
	c = cos(I->Xangle);
	s = sin(I->Xangle);
	S[0][0] =  c; S[0][1] = -s; S[0][2] =  s; S[0][3] =  c;
	S[1][0] =  s; S[1][1] =  c; S[1][2] = -c; S[1][3] =  s;
	S[2][0] =  c; S[2][1] =  s; S[2][2] =  s; S[2][3] = -c;
	S[3][0] = -s; S[3][1] =  c; S[3][2] =  c; S[3][3] =  s;

//	printf("S = "); printMatrix(S);

	return S;
}

/* This function takes the electric fields on the aperture and computes the
 * far field pattern via Fraunhofer diffraction.  See Jackson equation
 * 9.156
 */
void calcIllumpolparams(Illum *I)
{
	Matrix XRbeam, XLbeam, YRbeam, YLbeam;
	Matrix RRbeam, RLbeam, LRbeam, LLbeam;
	Matrix *beam, *stokes;
	Matrix T;
	const char stokename[][3] = {"I","Q","U","V","QI","UI","VI"};
	char filename[1000];
	int i, j, ii, jj, N;
	double kx, ky, kz, Rx, Ix, Ry, Iy;
	Matrix S;
	FILE *out;
	int start, stop;
	double f, factor, a, b, c;
	int F;
	
	if(I->savejonesmatrices == 0 && I->savestokesimages == 0)
		return;

	F = I->pixelsperbeam;

	S = calcSmatrix(I);

	N = MatrixSize1(I->A);

	beam = g_new(Matrix, 4);
	stokes = g_new(Matrix, 7);

	for(i = 0; i < 4; i++)
	{
		beam[i] = newMatrix(N, 2*N);
		zeroMatrix(beam[i]);
	}
	for(i = 0; i < 7; i++)
		stokes[i] = newMatrix(N, N);

	RRbeam = XRbeam = beam[0];
	RLbeam = XLbeam = beam[1];
	LRbeam = YRbeam = beam[2];
	LLbeam = YLbeam = beam[3];

	/* compute \hat{n} \cross \vec{E} on aperture surface 
	 *
	 * \hat{n} = \hat{z} = (0, 0, 1), so cross product does:
	 *	x -> y   y -> -x
	 */

	for(j = 0; j < N; j++) 
	{
		if(j >= N/2) jj = (j-N/2)/F + N/2;
		else jj = N/2 - 1 - (N/2-j)/F;
		for(i = 0; i < N; i++)
		{
			if(i >= N/2) ii = (i-N/2)/F + N/2;
			else ii = N/2 - 1 - (N/2-i)/F;
			b = 1.0-I->B[j][i];
			if(b <= 0.0) continue;

			/* components of beam excited with R pol */
			XRbeam[jj][2*ii]   += -b*I->E[RYR][j][i];
			XRbeam[jj][2*ii+1] += -b*I->E[IYR][j][i];
			YRbeam[jj][2*ii]   +=  b*I->E[RXR][j][i];
			YRbeam[jj][2*ii+1] +=  b*I->E[IXR][j][i];

			/* components of beam excited with L pol */
			XLbeam[jj][2*ii]   += -b*I->E[RYL][j][i];
			XLbeam[jj][2*ii+1] += -b*I->E[IYL][j][i];
			YLbeam[jj][2*ii]   +=  b*I->E[RXL][j][i];
			YLbeam[jj][2*ii+1] +=  b*I->E[IXL][j][i];
		}
	}

	for(i = 0; i < 4; i++)
	{
		rollMatrixinplace(beam[i], -N/2, -N);
		MatrixFFT(beam[i], FFTW_FORWARD);
		rollMatrixinplace(beam[i], N/2, N);
		/* convert from x, y to l, m by reflecting about y=m axis */
		reflectMatrix2inplace(beam[i]);
	}

	/* compute k cross the above FFT. 
	 *
	 * To recover stokes R and L, dot this vector with the conjugate of
	 * the circular polarization vectors.  p_r and p_l.
	 * These vectors need to be rotated
	 * at each point to be orthogonal to k.  The final product is then
	 * p_r^* . (k x I), where I is the Integral performed above in the
	 * FFT.  This is rearranged as -(k x p_r) . I.  In the circular basis
	 * k x p_r = -i p_r, and k x p_l = i p_l.  
	 */

	factor = I->lambda/(F*N*sqrt(I->dA));

	for(j = 0; j < N; j++) for(i = 0; i < 2*N; i+=2)
	{
		kx = (i/2-N/2)*(factor)-I->pointx;
		ky = (j-N/2)*(factor)-I->pointy;
		kz = sqrt(1.0-kx*kx-ky*ky);
		
		/* rotated p_r and p_l have the following components:
		 *   [note that the z component is not used]
		 * (a + b i) X + (b + c i) Y  and  (a - b i) X + (b - c i) Y
		 * k cross these are:
		 * (b - a i) X + (c - b i) Y  and  (b + a i) X + (c + b i) Y
		 * thus the conjugates used in the dotting are:
		 * (b + a i) X + (c + b i) Y  and  (b - a i) X + (c - b i) Y
		 */
		f = 1.0/(M_SQRT2*(kx*kx + ky*ky));
		a = f*(kx*kx*kz + ky*ky);
		b = f*(kx*ky - kx*ky*kz);
		c = f*(kx*kx + ky*ky*kz);

		/* compute RRbeam and LRbeam.  note that the beam matrices
		 * are being recycled
		 *
		 * The x <--> y here is done to reflect that radio astronomers
		 * have e_1 along the y = m axis and e_2 along the x=-l axis.
		 */
		Ry = XRbeam[j][i];
		Iy = XRbeam[j][i+1];
		Rx = YRbeam[j][i];
		Ix = YRbeam[j][i+1];

		RRbeam[j][i]   =  b*Rx - a*Ix + c*Ry - b*Iy;
		RRbeam[j][i+1] =  a*Rx + b*Ix + b*Ry + c*Iy;
		LRbeam[j][i]   =  b*Rx + a*Ix + c*Ry + b*Iy;
		LRbeam[j][i+1] = -a*Rx + b*Ix - b*Ry + c*Iy;

		/* compute RLbeam and LLbeam.  note that the beam matrices
		 * are being recycled
		 *
		 * The x <--> y here is done to reflect that radio astronomers
		 * have e_1 along the y = m axis and e_2 along the x=-l axis.
		 */
		Ry = XLbeam[j][i];
		Iy = XLbeam[j][i+1];
		Rx = YLbeam[j][i];
		Ix = YLbeam[j][i+1];

		RLbeam[j][i]   =  b*Rx - a*Ix + c*Ry - b*Iy;
		RLbeam[j][i+1] =  a*Rx + b*Ix + b*Ry + c*Iy;
		LLbeam[j][i]   =  b*Rx + a*Ix + c*Ry + b*Iy;
		LLbeam[j][i+1] = -a*Rx + b*Ix - b*Ry + c*Iy;
	}

	start = N/4 ;
	stop  = N - start + 1;
		
	if(I->savejonesmatrices)
	{
		sprintf(filename, "%s.jones.dat", I->prefix);
		out = fopen(filename, "w");
		for(jj = start; jj < stop; jj++) 
			for(ii = start; ii < stop; ii++)
		{
			fprintf(out, "%f %f %f %f %f %f %f %f\n",
				RRbeam[jj][2*ii],
				RRbeam[jj][2*ii+1],
				LRbeam[jj][2*ii],
				LRbeam[jj][2*ii+1],
				RLbeam[jj][2*ii],
				RLbeam[jj][2*ii+1],
				LLbeam[jj][2*ii],
				LLbeam[jj][2*ii+1]);
		}
		fclose(out);
	}

	for(j = 0; j < N; j++) for(i = 0; i < N; i++)
	{
		double RR, LL, RL;

		RR = RRbeam[j][2*i]   * RRbeam[j][2*i]   + 
		     RRbeam[j][2*i+1] * RRbeam[j][2*i+1] +
		     RLbeam[j][2*i]   * RLbeam[j][2*i]   +
		     RLbeam[j][2*i+1] * RLbeam[j][2*i+1] ;
		LL = LRbeam[j][2*i]   * LRbeam[j][2*i]   + 
		     LRbeam[j][2*i+1] * LRbeam[j][2*i+1] +
		     LLbeam[j][2*i]   * LLbeam[j][2*i]   +
		     LLbeam[j][2*i+1] * LLbeam[j][2*i+1] ;
		stokes[0][j][i] = 0.5*(RR + LL);
		stokes[3][j][i] = 0.5*(RR - LL);
		RL = RLbeam[j][2*i]   * LLbeam[j][2*i]   +
		     RLbeam[j][2*i+1] * LLbeam[j][2*i+1] +
		     RRbeam[j][2*i]   * LRbeam[j][2*i]   +
		     RRbeam[j][2*i+1] * LRbeam[j][2*i+1] ;
		stokes[1][j][i] = RL;
		RL = RLbeam[j][2*i+1] * LLbeam[j][2*i]   -
		     RLbeam[j][2*i]   * LLbeam[j][2*i+1] +
		     RRbeam[j][2*i+1] * LRbeam[j][2*i]   -
		     RRbeam[j][2*i]   * LRbeam[j][2*i+1];
		stokes[2][j][i] = RL;

		stokes[4][j][i] = stokes[1][j][i]/stokes[0][j][i];
		stokes[5][j][i] = stokes[2][j][i]/stokes[0][j][i];
		stokes[6][j][i] = stokes[3][j][i]/stokes[0][j][i];
	}

	if(I->savestokesimages)
	{
		printf("\n");
		printf("Output image scale is %f deg/pixel\n", 
			180.0*factor/M_PI);

		I->beampixelscale = factor;
	
		for(i = 0; i < 7; i++)
		{
			sprintf(filename, "%s.%s.pgm", I->prefix, stokename[i]);
			/* flip again, because pgm stores the image upside down
			 */
			reflectMatrix1inplace(stokes[i]);
			T = refsubMatrix(stokes[i], start, start, 
				stop-1, stop-1);
			saveMatrixaspgm(T, filename);
			deleteMatrix(T);
		}
	}

	/* cleanup */
	for(i = 0; i < 7; i++) deleteMatrix(stokes[i]);
	g_free(stokes);
	for(i = 0; i < 4; i++) deleteMatrix(beam[i]);
	g_free(beam);

	deleteMatrix(S);
}

static double peaksidelobe(const Matrix M, int j0, int i0)
{
	Matrix A;
	int i, j, imax, n, o;
	int rx=0, ry=0, rx1=0, rx2=0, ry1=0, ry2=0;
	double peak, side;

	/* 0. build centered amplitude image */
	A = Matrixcomplexamplitudes(M);
	n = MatrixSize1(A);
	o = n/2;
	if(j0 > o) j0 -= n;
	if(i0 > o) i0 -= n;
	rollMatrixinplace(A, o-j0, o-i0);
	peak = A[o][o];
	
	/* 1. find first null in x direction */
	for(i = 2; i < o-2; i++)
	{
		if(rx1 == 0 && A[o][o+i] < A[o][o+i-1] && 
			       A[o][o+i] < A[o][o+i+1]) rx1 = i;
		if(rx2 == 0 && A[o][o-i] < A[o][o-i-1] && 
			       A[o][o-i] < A[o][o-i+1]) rx2 = i;
		if(rx1 != 0 && rx2 != 0) break;
	}
	if(rx1 > rx2) rx = rx1; else rx = rx2;
	
	/* 2. find first null in y direction */
	for(j = 2; j < o-2; j++)
	{
		if(ry1 == 0 && A[o+j][o] < A[o+j-1][o] && 
			       A[o+j][o] < A[o+j+1][o]) ry1 = j;
		if(ry2 == 0 && A[o-j][o] < A[o-j-1][o] && 
			       A[o-j][o] < A[o-j+1][o]) ry2 = j;
		if(ry1 != 0 && ry2 != 0) break;
	}
	if(ry1 > ry2) ry = ry1; else ry = ry2;

	/* 3. blank ellipse in middle */
	for(j = -ry; j <= ry; j++)
	{
		imax = sqrt(ry*ry - j*j)*(double)rx/(double)ry;
		for(i = -imax; i <= imax; i++) A[o+j][o+i] = 0.0;
	}

	side = Matrixmax(A)/peak;

	deleteMatrix(A);

	return side*side;
}

/* Todo : take F as a parameter, automatically peak up on area that scales
 * with F. */
void dephaseIllum(Illum *I)
{
	Matrix M, A;
	int F=6;
	double max=0.0, mag=0.0;
	int maxi=0, maxj=0, i, j, k, n, ii, jj;
	double dx, xfac, yfac;
	double di, dj;
	double pointx, pointy;
	double c, s, re, im;
	double Axx, Axy, Ayy;
	double sidelobe;

	/* Make an F times oversampled beam... */
	
	n = MatrixSize1(I->A);

	M = newMatrix(n*F, 2*n*F);
	zeroMatrix(M);
	for(j = 0; j < n; j++) for(i = 0; i < n; i++)
	{
		M[j][2*i]   = (1.0-I->B[j][i])*I->A[j][i]*cos(I->P[j][i]);
		M[j][2*i+1] = (1.0-I->B[j][i])*I->A[j][i]*sin(I->P[j][i]);
	}

	rollMatrixinplace(M, -n/2, -n);
	MatrixFFT(M, FFTW_FORWARD);
	
	for(j = 0; j < n*F; j++) for(i = 0; i < n*F; i++)
	{
		mag = M[j][2*i]*M[j][2*i]+M[j][2*i+1]*M[j][2*i+1];
		if(mag > max)
		{
			max = mag;
			maxi = i;
			maxj = j;
		}
	}

	dx = sqrt(I->dA);

	if(maxi > F*n/2) maxi = maxi-F*n;
	if(maxj > F*n/2) maxj = maxj-F*n;

	/* Peak up the maximum */
	A = newMatrix(5, 5);
	for(j = 0; j < 5; j++) for(i = 0; i < 5; i++)
	{
		ii = maxi-2+i;
		if(ii < 0) ii += n*F;
		jj = maxj-2+j;
		if(jj < 0) jj += n*F;
		A[j][i] = M[jj][2*ii]*M[jj][2*ii]+M[jj][2*ii+1]*M[jj][2*ii+1];
	}

	scaleMatrix(A, 1.0/A[2][2]);
	applyfunctoMatrix(A, log);

	Matrixpeakup(A, 2, 2, 2, 2, &dj, &di, &Ayy, &Axy, &Axx);
	I->fwhm_x = 2.0*sqrt(-log(2.0)/Axx)*I->lambda/(F*n*dx);
	I->fwhm_y = 2.0*sqrt(-log(2.0)/Ayy)*I->lambda/(F*n*dx);

	di -= 2.0;
	dj -= 2.0;

	if(fabs(di) > 0.8 || fabs(dj) > 0.8)
	{
		fprintf(stderr, "Warning -- peak up interpolation fails\n");
		di = dj = 0.0;
	}

	/* incremental pointing offsets in radians */
	pointx = (maxi+di)*I->lambda/(F*n*dx);
	pointy = (maxj+dj)*I->lambda/(F*n*dx);

	if(fabs(pointx) < 0.000001 && fabs(pointy) < 0.000001)
	{
		pointx = 0.000001;
	}

	/* should there be a sin() below around point? */
	xfac = 2.0*M_PI*pointx*dx/I->lambda;
	yfac = 2.0*M_PI*pointy*dx/I->lambda;

	/* subtract phase wedge from scalar phase */
	for(j = 0; j < n; j++) for(i = 0; i < n; i++)
		I->P[j][i] -= (xfac*i + yfac*j);

	/* subtract phase wedge from vector data */
	for(j = 0; j < n; j++) for(i = 0; i < n; i++)
	{
		c = cos(xfac*i + yfac*j);
		s = -sin(xfac*i + yfac*j);
		for(k = 0; k < 12; k+=2)
		{
			re = I->E[k  ][j][i];
			im = I->E[k+1][j][i];
			I->E[k  ][j][i] = re*c - im*s;
			I->E[k+1][j][i] = re*s + im*c;
		}
	}

	/* accumulated pointing offsets in radians */
	I->pointx += pointx;
	I->pointy += pointy;

	sidelobe = peaksidelobe(M, maxj, maxi);
	if(sidelobe > I->peaksidelobe) I->peaksidelobe = sidelobe;

	deleteMatrix(M);
}

void printIllum(const Illum *I)
{
	printf("Output: %p\n", I);
	if(I->spilleff > 0.0)
	{
		printf("  Spillover eff = %f\n", I->spilleff);
		if(I->subspilleff != 0.0 > I->prispilleff > 0.0)
		{
			printf("    primary     = %f\n", I->prispilleff);
			printf("    subreflector= %f\n", I->subspilleff);
		}
	}
	if(I->blockeff > 0.0)
		printf("  Blockage eff  = %f\n", I->blockeff);
	if(I->surfeff > 0.0)
		printf("  Surface eff   = %f\n", I->surfeff);
	if(I->illumeff > 0.0)
	{
		printf("  Illum eff     = %f\n", I->illumeff);
		if(I->phaseeff > 0.0)
			printf("    phase eff   = %f\n", I->phaseeff);
		if(I->ampeff > 0.0)
			printf("    amp eff     = %f\n", I->ampeff);
	}
	if(I->diffeff > 0.0)
		printf("  Diffract eff  = %f\n", I->diffeff);
	if(I->misceff > 0.0)
		printf("  Misc eff      = %f\n", I->misceff);
	if(I->efficiency > 0.0)
		printf("  Total eff     = %f\n", I->efficiency);
	if(I->gain > 0.0)
		printf("  Gain          = %6.2f  =  %5.2f dBi\n",
			I->gain, 10.0*log10(I->gain));
	if(I->Tsys > 0.0)
	{
		printf("  Tsys          = %6.3f K\n",  I->Tsys);
		if(I->pTground > 0.0)
			printf("    ground      = %6.3f K\n",  I->pTground);
		if(I->pTsky > 0.0)
			printf("    sky         = %6.3f K\n",  I->pTsky);
		if(I->pTrec > 0.0)
			printf("    rec         = %6.3f K\n",  I->pTrec);
	}
	if(I->aeff > 0.0)
	{
		printf("  Aeff          = %f m^2\n",   I->aeff);
		if(I->Tsys > 0.0)
			printf("  Aeff/Tsys     = %f m^2/K\n", I->aeff/I->Tsys);
	}
	if(I->pointx != 0.0 || I->pointy != 0.0)
	{
		printf("  l beamshift   = %f deg\n", -I->pointx*180.0/M_PI);
		printf("  m beamshift   = %f deg\n", I->pointy*180.0/M_PI);
	}
	if(I->fwhm_x != 0.0 || I->fwhm_y != 0.0)
	{
		printf("  l beam FWHM   = %f deg\n", I->fwhm_x*180.0/M_PI);
		printf("  m beam FWHM   = %f deg\n", I->fwhm_y*180.0/M_PI);
	}
	if(I->peaksidelobe > 0.0)
		printf("  Peak sidelobe = %f = %f dB\n", I->peaksidelobe,
			10.0*log10(I->peaksidelobe));
}

void KeyValueupdateIllum(struct KeyValue *kv, const Illum *I)
{
	if(I->spilleff > 0.0)
	{
		KeyValueupdateparmdouble(kv, "spilleff",  I->spilleff);
		if(I->subspilleff != 0.0 > I->prispilleff > 0.0)
		{
		   KeyValueupdateparmdouble(kv, "prispilleff",  I->prispilleff);
		   KeyValueupdateparmdouble(kv, "subspilleff",  I->subspilleff);
		}
	}
	if(I->blockeff > 0.0)
		KeyValueupdateparmdouble(kv, "blockeff",  I->blockeff);
	if(I->surfeff > 0.0)
		KeyValueupdateparmdouble(kv, "surfeff",   I->surfeff);
	if(I->illumeff > 0.0)
	{
		KeyValueupdateparmdouble(kv, "illumeff",  I->illumeff);
		if(I->phaseeff > 0.0)
			KeyValueupdateparmdouble(kv, "phaseeff",  I->phaseeff);
		if(I->ampeff > 0.0)
			KeyValueupdateparmdouble(kv, "ampeff",    I->ampeff);
	}
	if(I->diffeff > 0.0)
		KeyValueupdateparmdouble(kv, "diffeff", I->diffeff);
	if(I->misceff > 0.0)
		KeyValueupdateparmdouble(kv, "misceff", I->misceff);
	if(I->efficiency > 0.0)
		KeyValueupdateparmdouble(kv, "totaleff",  I->efficiency);
	if(I->gain > 0.0)
		KeyValueupdateparmdouble(kv, "gain", I->gain);
	if(I->Tsys > 0.0)
		KeyValueupdateparmdouble(kv, "Tsys", I->Tsys);
	if(I->aeff > 0.0)
	{
		KeyValueupdateparmdouble(kv, "Aeff",      I->aeff);
		if(I->Tsys > 0.0)
			KeyValueupdateparmdouble(kv, "Aeff_Tsys", 
				I->aeff/I->Tsys);
	}
	if(I->pointx != 0.0 || I->pointy != 0.0)
	{
		KeyValueupdateparmdouble(kv, "point_l", -I->pointx*180.0/M_PI);
		KeyValueupdateparmdouble(kv, "point_m", I->pointy*180.0/M_PI);
	}
	if(I->fwhm_x != 0.0 || I->fwhm_y != 0.0)
	{
		KeyValueupdateparmdouble(kv, "fwhm_l", I->fwhm_x*180.0/M_PI);
		KeyValueupdateparmdouble(kv, "fwhm_m", I->fwhm_y*180.0/M_PI);
	}
	if(I->peaksidelobe > 0.0)
		KeyValueupdateparmdouble(kv, "peaksidelobe", I->peaksidelobe);
	if(I->beampixelscale > 0.0)
		KeyValueupdateparmdouble(kv, "beampixelscale", 
			180.0*I->beampixelscale/M_PI);
}
