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
#ifndef __VECARRAY__
#define __VECARRAY__

#include "vector.h"

#ifdef FLOATVEC
typedef  float** VecArray;
#else
typedef double** VecArray;
#endif

#define VecArraySize(a) (((int *)(a))[-1])

VecArray newVecArray(int n);
VecArray newpopulatedVecArray(int n, int m);
//VecArray newVecArrayfromVectors(...);
void zeroVecArray(VecArray a);
void deleteVecArray(VecArray a);
void deleteVecArrayandVectors(VecArray a);
VecArray dupVecArray(const VecArray a);
VecArray dupVecArrayandVectors(const VecArray a);
void copytoVecArray(VecArray dest, const VecArray src);
VecArray VecArraysubset(const VecArray a, int col, int ncol);
VecArray dupVecArraysubset(const VecArray a, int col, int ncol);
VecArray dupVecArrayrowsubset(const VecArray a, int row, int nrow);
int VecArrayVectorSize(const VecArray a);
int VecArrayVectorcolumns(const VecArray a);
VecArray VecArraydeletecolumn(VecArray a, int col);
VecArray VecArrayinsertcolumn(VecArray a, int col);
VecArray VecArrayappendcolumn(VecArray a);
VecArray VecArrayinsertVector(VecArray a, Vector v, int col);
VecArray VecArrayappendVector(VecArray a, Vector v);
VecArray VecArrayappendVectors(VecArray a, ...);
VecArray appendtoVecArray(VecArray a, VecArray b);
void VecArraysetcolumn(VecArray a, Vector v, int col);
Matrix newMatrixfromVecArray(const VecArray a);
VecArray newVecArrayfromMatrix(const Matrix M);
VecArray VecArrayconcat(VecArray a, VecArray b);
VecArray VecArrayappendVecArray(VecArray a, const VecArray b);
VecArray VecArrayappendVecArrays(VecArray a, ...);
VecArray reduceVecArraybymask(VecArray a, const intVector mask);
VecArray VecArrayfromfile(const char *filename, int mincols);
int VecArraytofile(const VecArray a, const char *filename, 
	const char *format);
void VecArrayresize(const VecArray a, int newsize);

/* arithmetic */
void scaleVecArray(VecArray a, double f);
void addtoVecArray(VecArray a, const VecArray b);
Vector VecArrayVectoraverage(const VecArray a);

/* tridiagonal matrix stuff */
VecArray tridiagfromMatrix(Matrix a);

#endif
