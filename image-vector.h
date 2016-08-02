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
#ifndef __IMAGE_VECTOR__
#define __IMAGE_VECTOR__

#include "image.h"
#include "vector.h"

Image newImagefromMatrix(const Matrix M);
void saveMatrixaspgm(const Matrix M, const char *filename);
Matrix newMatrixfromImage(const Image im);
Matrix newpaddedMatrixfromImage(const Image im, int rowpad);

#endif
