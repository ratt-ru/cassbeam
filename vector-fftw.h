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
#ifndef __VECTOR_FFTW_H__
#define __VECTOR_FFTW_H__

#include <math.h>
#include <fftw.h>
#include <rfftw.h>

#include "vector.h"

#define SMOOTH_GAUSS	1
#define SMOOTH_COS	2
#define SMOOTH_COS1	3
#define SMOOTH_BOX	4
#define SMOOTH_TRIANGLE	5

/* This function grids a non-evenly sampled datum onto the FT plane */
int Matrixcomplexgrid(Matrix M, double u, double v, double re, double im,int p);

/* These use FFTW for the transform.  Don't use these if many fast transforms
 * are needed -- use plans and the fftw functions in those cases 
 *
 * dir is one of: FFTW_FORWARD or FFTW_BACKWARD for complex (FFT) transforms
 *            or: FFTW_REAL_TO_COMPLEX or FFTW_REAL_TO_COMPLEX for real (RFFT)
 */
Vector VectorFFT(Vector V, fftw_direction dir);
Vector VectorRFFT(Vector V, fftw_direction dir);

void smoothVector(Vector V, double width, int type);

/* Notes for Matrix transforms: 
 *	1. all real side of RFFT transforms are represented in FFTW's 
 *		native order
 *	2. input arrays for all RFFT operations are assumed to have 2 extra
 *		floating point values padding each row.  These matrices can be
 *		created with newpaddedMatrix(Ny, Nx, pad = 2)
 */
void MatrixFFT(Matrix M, fftw_direction dir);
void MatrixFFTrows(Matrix M, fftw_direction dir);
void MatrixFFTcolumns(Matrix M, fftw_direction dir);
void MatrixRFFT(Matrix M, fftw_direction dir);
void MatrixRFFTrows(Matrix M, fftw_direction dir);
Matrix Matrixcrosscorrelate(const Matrix M, const Matrix N);
Matrix Matrixautocorrelate(const Matrix M);
void Matrixconvolve(Matrix M, const Matrix N);
void Matrixdeconvolve(Matrix M, const Matrix N);

#endif
