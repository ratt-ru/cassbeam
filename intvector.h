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
#ifndef __INT_VECTOR__
#define __INT_VECTOR__


typedef int*     intVector;
typedef int      intVectorType;
typedef int*     intVectorTypePointer;
typedef int**	 intMatrix;
typedef int      intMatrixType;
typedef int*     intMatrixTypePointer;


typedef struct
{
	intMatrixTypePointer datastart;
	int n, m, rowpad;
} intMatrixStubType;

#define INTMATRIXALIGN             32

#define INTMATRIXSTUBBLOCKS     (((sizeof(intMatrixStubType)-1)/INTMATRIXALIGN)+1)
#define INTMATRIXSTUBINDICES    (((INTMATRIXALIGN*INTMATRIXSTUBBLOCKS-1)/ \
                                     sizeof(intMatrixTypePointer))+1)
#define INTMATRIXVALUESPERBLOCK (((INTMATRIXALIGN-1)/sizeof(intMatrixType))+1)
#define intMatrixStub(M)        ((intMatrixStubType *)(M-INTMATRIXSTUBINDICES))
#define constintMatrixStub(M)   ((const intMatrixStubType *)(M-INTMATRIXSTUBINDICES))


#define intVectorSize(v)	((v)[-1])
#define intVectorSizeSame(U, V)	(intVectorSize((U)) == intVectorSize((V)))

#define intMatrixSize1(M)       (constintMatrixStub((M))->n)
#define intMatrixSize2(M)       (constintMatrixStub((M))->m)
#define intMatrixrowpad(M)      (constintMatrixStub((M))->rowpad)
#define intMatrixstride(M)      (((M)[1]-(M)[0]))
#define intMatrixrefcount(M)    (*((int *)((intMatrixStub((M))->datastart)- \
                                     INTMATRIXVALUESPERBLOCK)))


#define newintMatrix(n, m)	newpaddedintMatrix((n), (m), 0)
#define refintMatrix(M)		refsubintMatrix((M), 0, 0, -1, -1)

#define intMatrixSizeSame(M, N) ((intMatrixSize1((M)) == intMatrixSize1((N))) && \
                                 (intMatrixSize2((M)) == intMatrixSize2((N))))


intVector newintVector(int n);
void zerointVector(intVector v);
void fillintVector(intVector v, int r);
void deleteintVector(intVector v);
intVector dupintVector(const intVector v);
void printintVector(const intVector v);
void saveintVectorasascii(const intVector v, const char *filename);
int intVectorrandomelement(const intVector v);
int intVectormax(const intVector v);
int intVectormin(const intVector v);

intMatrix newpaddedintMatrix(int n, int m, int rowpad);
void deleteintMatrix(intMatrix M);
intMatrix refsubintMatrix(const intMatrix M, int n1, int m1, int n2, int m2);
void zerointMatrix(intMatrix M);
intMatrix dupintMatrix(const intMatrix M);
intMatrixType intMatrixmax(const intMatrix M);
intMatrixType intMatrixmin(const intMatrix M);
void intMatrixminmax(const intMatrix M, intMatrixType *min, intMatrixType *max);
void scaleintMatrix(intMatrix M, double f);
void biasintMatrix(intMatrix M, int b);
void addtointMatrix(intMatrix M, const intMatrix N);
intMatrix addintMatrices(const intMatrix M, const intMatrix N);
void copytointMatrix(intMatrix M, const intMatrix N);
void applyfunctointMatrix(intMatrix M, int (*func)(int x));
intMatrix intMatrixproduct(intVector v2, intVector v1);
intMatrix dupsubintMatrix(const intMatrix M, int n1, int m1, int n2, int m2);
void copysubintMatrix(intMatrix N, const intMatrix M, int n1, int m1,
	int n2, int m2, int n3, int m3);
intVector intMatrixrow(const intMatrix M, int row);
intVector intMatrixcolumn(const intMatrix M, int column);
intVector sumintMatrixrows(const intMatrix M);
intVector sumintMatrixcolumns(const intMatrix M);
void intMatrixmeanRMS(const intMatrix M, double *rms, double *mean);
double intMatrixmean(const intMatrix M);
double intMatrixRMS(const intMatrix M);
intMatrix intMatrixcomplexamplitudes(const intMatrix M);
intMatrix intMatrixcomplexphases(const intMatrix M);
intVector intMatrixhistogram(const intMatrix M, int bin0, int bin1);




#endif  /* __INT_VECTOR__ */
