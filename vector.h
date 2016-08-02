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
#ifndef __VECTOR__
#define __VECTOR__

#include "intvector.h"

/* TODO: Need to speed up matrix functions by defining Mrow = M[j] */


/* Uncomment next line to use floats instead of doubles */
// #define FLOATVECTOR */

#ifdef FLOATVECTOR
typedef  float*  Vector;
typedef  float   VectorType;
typedef  float*  VectorTypePointer;
typedef  float** Matrix;
typedef  float   MatrixType;
typedef  float*  MatrixTypePointer;
#else
typedef double*  Vector;
typedef double   VectorType;
typedef double*  VectorTypePointer;
typedef double** Matrix;
typedef double   MatrixType;
typedef double*  MatrixTypePointer;
#endif

typedef struct 
{
	MatrixTypePointer datastart;
	int n, m, rowpad;
} MatrixStubType;

#define MATRIXALIGN		32

#define MATRIXSTUBBLOCKS	(((sizeof(MatrixStubType)-1)/MATRIXALIGN)+1)
#define MATRIXSTUBINDICES	(((MATRIXALIGN*MATRIXSTUBBLOCKS-1)/ \
					sizeof(MatrixTypePointer))+1)
#define MATRIXVALUESPERBLOCK	(((MATRIXALIGN-1)/sizeof(MatrixType))+1)
#define MatrixStub(M)		((MatrixStubType *)(M-MATRIXSTUBINDICES))
#define constMatrixStub(M)	((const MatrixStubType *)(M-MATRIXSTUBINDICES))

#define VectorSize(v)   	(((int *)(v))[-1])
#define Vectorlength(v) 	(((const int *)(v))[-1])

#define MatrixSize1(M)  	(constMatrixStub((M))->n)
#define MatrixSize2(M)  	(constMatrixStub((M))->m)
#define Matrixrowpad(M)		(constMatrixStub((M))->rowpad)
#define Matrixstride(M)		(((M)[1]-(M)[0]))
#define Matrixrefcount(M)	(*((int *)((MatrixStub((M))->datastart)- \
					MATRIXVALUESPERBLOCK)))

#define newMatrix(n, m)	newpaddedMatrix((n), (m), 0)
#define refMatrix(M)	refsubMatrix((M), 0, 0, -1, -1)

#define VectorSizeSame(U, V)	(Vectorlength((U)) == Vectorlength((V)))

#define MatrixSizeSame(M, N)	((MatrixSize1((M)) == MatrixSize1((N))) && \
				 (MatrixSize2((M)) == MatrixSize2((N))))
#define Matrixissquare(M)	 (MatrixSize1((M)) == MatrixSize2((M)))

int gettotalVectordata();

int gaussjordan(Matrix a, Vector b);
Vector polyfit(const Vector X, const Vector Y, int order);
Vector weightedpolyfit(const Vector X, const Vector Y, const Vector W, 
	int order);
Vector polyfitwithlag(const Vector X, const Vector Y, int order, int lag);


Vector newVector(int n);
Vector newVectorfromarray(int n, const double *arr);
Vector newVectorfromintVector(intVector v);
void zeroVector(Vector v);
void fillVector(Vector v, double value);
Vector dupVector(const Vector v);
void deleteVector(Vector v);
int Vectorisfinite(const Vector v);
double dotVectors(const Vector v1, const Vector v2);
void subfromVector(Vector a, const Vector b);
void addtoVector(Vector a, const Vector b);
void addscaledVectortoVector(Vector a, const Vector b, double scale);
Vector addVectors(const Vector a, const Vector b);
Vector subVectors(const Vector a, const Vector b);
void copytoVector(Vector a, const Vector b);
void scaleVector(Vector v, double f);
void biasVector(Vector v, double f);
void squareVector(Vector v);
void normalizeVector(Vector v);
void printVector(const Vector v);
void applyfunctoVector(Vector v, double (*func)(double x));
VectorType Vectormax(const Vector v);
VectorType Vectormin(const Vector v);
int Vectorfindmax(const Vector v);
int Vectorfindmin(const Vector v);
void Vectorminmax(const Vector v, VectorType *min, VectorType *max);
VectorType interpolateVector(const Vector V, double index);
void saveVectorasascii(const Vector v, const char *filename);
void saveVectorasformattedascii(const Vector v, const char *filename,
	const char *format);
void Vectorsavebinary(const Vector v, const char *filename);
Vector Vectorloadbinary(const char *filename);
void VectormeanRMS(const Vector v, double *mean, double *rms);
double Vectorsum(const Vector v);
double Vectorsumsquare(const Vector v);
double Vectormean(const Vector v);
double Vectormeansquare(const Vector v);
double VectorRMS(const Vector v);
double Vectormoment(const Vector v, int moment);
void rollVector(Vector v, int deltan);
Vector rebinVector(const Vector src, int factor);
void Vectorboxcar(const Vector src, Vector dest, int s);
double VectorAllanStandardDeviation(const Vector v, int lag, int pad, int step);
intVector ComputeAllanLags(int n, double power);
Vector VectorAllanStandardDeviations(const Vector v, const intVector lags);
int stringtoVector(const char *str, Vector v);
char *Vectortostring(const Vector v);
Vector newVectorfromstring(const char *str);
Vector *Vectorcolumnsfromfile(const char *filename, int mincols, int *ncol);
int Vectorcolumnstofile(const Vector *data, const char *filename, int ncol);
int Vectorcolumnstofilewithformat(const Vector *data, const char *filename, int ncol,
        const char *format);
int Vectorcolumnstofilewithmask(const Vector *data, const Vector mask,
	const char *filename, int ncol);
int Vectorcolumnstofilewithgaps(const Vector *data, const char *filename, 
	int ncol, int refcol, double gap);
Vector subVector(const Vector v, int n1, int n2);
Vector newunitVector(const Vector v);

int Vectorisincreasing(const Vector v);
int Vectorisdecreasing(const Vector v);

Matrix newpaddedMatrix(int n, int m, int rowpad);
void deleteMatrix(Matrix M);
Matrix refsubMatrix(const Matrix M, int n1, int m1, int n2, int m2);
void zeroMatrix(Matrix M);
Matrix dupMatrix(const Matrix M);
void printMatrix(const Matrix M);
Matrix repadMatrix(const Matrix M, int pad);
void Matrixsavebinary(const Matrix M, const char *filename);
Matrix Matrixloadbinary(const char *filename);
MatrixType Matrixmax(const Matrix M);
MatrixType Matrixmaxvalue(const Matrix M, int *bestj, int *besti);
MatrixType Matrixmin(const Matrix M);
MatrixType Matrixminvalue(const Matrix M, int *bestj, int *besti);
void Matrixminmax(const Matrix M, MatrixType *min, MatrixType *max);
double Matrixpeakup(const Matrix D, int y, int x, int py, int px,
        double *y0, double *x0, double *Ayy, double *Axy, double *Axx);
void scaleMatrix(Matrix M, double f);
void addscaledMatrixtoMatrix(Matrix a, const Matrix b, double scale);
void biasMatrix(Matrix M, double f);
void addtoMatrix(Matrix M, const Matrix N);
void subfromMatrix(Matrix M, const Matrix N);
Matrix addMatrices(const Matrix M, const Matrix N);
void copytoMatrix(Matrix M, const Matrix N);
void applyfunctoMatrix(Matrix M, double (*func)(double x));
Matrix Matrixproduct(Vector v2, Vector v1);
Matrix dupsubMatrix(const Matrix M, int n1, int m1, int n2, int m2);
void copysubMatrix(Matrix N, const Matrix M, int n1, int m1, int n2, int m2,
	int n3, int m3);
Vector Matrixrow(const Matrix M, int row);
Vector Matrixcolumn(const Matrix M, int column);
Vector sumMatrixrows(const Matrix M);
Vector sumMatrixcolumns(const Matrix M);
void MatrixmeanRMS(const Matrix M, double *rms, double *mean);
double Matrixmean(const Matrix M);
double MatrixRMS(const Matrix M);
Matrix Matrixmultiply(const Matrix A, const Matrix B);
void copyMatrixmultiply(Matrix C, const Matrix A, const Matrix B);
Vector MatrixVectormultiply(const Matrix M, const Vector V);
void Matrixmultiplyelements(Matrix M, const Matrix N);
void Matrixdivideelements(Matrix M, const Matrix N);
Matrix Matrixcomplexamplitudes(const Matrix M);
Matrix Matrixcomplexphases(const Matrix M);
void Matrixcomplexmultiply(Matrix M, const Matrix N, int pad);
void Matrixcomplexconjugatemultiply(Matrix M, const Matrix N, int pad);
void Matrixcomplexdivide(Matrix M, const Matrix N, int pad);
Vector Matrixhistogram(const Matrix M, int bins);
Matrix quarterbinMatrix(const Matrix M);
Matrix rollMatrix(const Matrix M, int delta_n, int delta_m);
void rollMatrixinplace(Matrix M, int delta_n, int delta_m);
Matrix transposeMatrix(const Matrix M);
void transposeMatrixinplace(Matrix M);
void reflectMatrix1inplace(Matrix M);
void reflectMatrix2inplace(Matrix M);
void symmetrizeMatrixinplace(Matrix M);
void antisymmetrizeMatrixinplace(Matrix M);
Matrix newrotationMatrix(double rx, double ry, double rz);
void clipMatrix(Matrix M, double min, double max);
Matrix LMatrixinvert(const Matrix L);
Matrix choleskydecompose(const Matrix A);
Matrix choleskydecomposeld(const Matrix A);
Vector cholesky(const Matrix A, const Vector B);
int choleskyinvert(Matrix A);

int Matrixisposdef(const Matrix M);
int Matrixisposdef2(const Matrix M, int niter);
int Matrixisfinite(const Matrix M);
int Matrixislower(const Matrix M);
int Matrixisupper(const Matrix M);
int Matrixistridiag(const Matrix M);
int Matrixissymmetric(const Matrix M);

#endif
