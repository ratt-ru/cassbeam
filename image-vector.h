#ifndef __IMAGE_VECTOR__
#define __IMAGE_VECTOR__

#include "image.h"
#include "vector.h"

Image newImagefromMatrix(const Matrix M);
void saveMatrixaspgm(const Matrix M, const char *filename);
Matrix newMatrixfromImage(const Image im);
Matrix newpaddedMatrixfromImage(const Image im, int rowpad);

#endif
