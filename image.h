#ifndef __IMAGE_H__
#define __IMAGE_H__

/* data[Y][X] */

struct _Image
{
	unsigned char **data;
	int xres, yres;
};

typedef struct _Image* Image;

#define zeroImage(p)  blankImage(p, 0)

Image newImage(int xres, int yres);
void deleteImage(Image p);

Image loadImageaspgm(const char *filename);
Image loadImageasweightedpgm(const char *filename, 
	float rw, float gw, float bw);
int saveImageaspgm(const Image p, const char *filename);

Image dupImage(const Image p);
void blankImage(Image p, unsigned char value);
void Imagedrawhline(Image p, int x0, int x1, int y, unsigned char value);
void stretchImage(Image p);
Image scaleImage(const Image p, int newxres, int newyres);
void applyfunctoImage(Image p, unsigned char (* func)(unsigned char v));
Image halveImage(Image p);

#endif
