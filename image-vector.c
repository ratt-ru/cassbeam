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
#include "image-vector.h"

Image newImagefromMatrix(const Matrix M)
{
	Image im;
	int m, n, i, j;
	MatrixType min, max;

	n = MatrixSize1(M);
	m = MatrixSize2(M);

	im = newImage(m, n);

	Matrixminmax(M, &min, &max);

	for(i = 0; i < m; i++) for(j = 0; j < n; j++)
		im->data[j][i] = 255.0*(M[j][i]-min)/(max-min);
	
	return im;
}

Matrix newMatrixfromImage(const Image im)
{
	Matrix M;
	int m, n, i, j;
	m = im->xres;
	n = im->yres;

	M = newMatrix(n, m);
	for(i = 0; i < m; i++) for(j = 0; j < n ; j++) 
		M[j][i] = im->data[j][i]/255.0;
	
	return M;
}

void saveMatrixaspgm(const Matrix M, const char *filename)
{
	Image I;

	g_assert(M);
	I = newImagefromMatrix(M);
	saveImageaspgm(I, filename);
	deleteImage(I);
}

/* Note -- this doesn't blank the pad */
Matrix newpaddedMatrixfromImage(const Image im, int rowpad)
{
	Matrix M;
	int m, n, i, j;
	m = im->xres;
	n = im->yres;

	M = newpaddedMatrix(n, m, rowpad);
	for(i = 0; i < m; i++) for(j = 0; j < n ; j++) 
		M[j][i] = im->data[j][i]/255.0;
	
	return M;
}
