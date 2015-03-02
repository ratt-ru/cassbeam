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
