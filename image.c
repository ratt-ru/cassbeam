#include <glib.h>
#include <string.h>
#include <stdio.h>
#include "image.h"

/* Warning!  Does not blank Image */
Image newImage(int xres, int yres)
{
	int i;
	Image img;
	unsigned char *raw;
	
	img = g_new(struct _Image, 1);
	img->xres = xres;
	img->yres = yres;
	img->data = g_new(unsigned char *, yres);
	raw = g_new(unsigned char, xres*yres);
	for(i = 0; i < yres; i++) img->data[i] = raw+(i*xres);
	return img;
}

void deleteImage(Image p)
{
	if(!p) 
	{
		fprintf(stderr, "Warning: Deleteing nullImage\n");
		return;
	}
	if(p->data[0]) g_free(p->data[0]);
	if(p->data) g_free(p->data);
	g_free(p);
}

/* 0 on success */
int readpgmline(FILE *in, char line[100], const char *filename)
{
	for(;;)
	{
		fgets(line, 99, in);
		if(feof(in))
		{
			fprintf(stderr, "Short file: %s\n", filename);
			fclose(in);
			return -1;
		}
		if(strlen(line) > 98)
		{
			fprintf(stderr, "loadImageaspgm: Bad file: %s\n",
				filename);
			fclose(in);
			return -1;
		}
		if(line[0] != '#') 
			return 0;
	}
}

Image loadImageaspgm(const char *filename)
{
	Image img;
	FILE *in;
	int i, j, k, xres, yres, maxgrey;
	char line[100];
	int mode, c;
	guchar *row;

	in = fopen(filename, "r");
	if(!in)
	{
		fprintf(stderr, "loadImageaspgm: File not found: %s\n",
			filename);
		return 0;
	}
	if(readpgmline(in, line, filename)) return 0;
	if(line[0] != 'P' || (line[1] != '5' && line[1] != '6'))
	{
		fprintf(stderr, "Unsupported file: %s\n", filename);
		fclose(in);
		return 0;
	}
	if(line[1] == '5') mode = 0;
	else mode = 1;
	if(readpgmline(in, line, filename)) return 0;
	sscanf(line, "%d%d", &xres, &yres);
	if(xres <= 0 || yres <= 0)
	{
		fprintf(stderr, "loadImageaspgm: Bad file: %s\n", filename);
		fclose(in);
		return 0;
	}
	if(readpgmline(in, line, filename)) return 0;
	sscanf(line, "%d", &maxgrey);

	img = newImage(xres, yres);
	if(mode == 0) 
	{
		for(j = 0; j < yres; j++) 
		{
			if(fread(img->data[j], sizeof(unsigned char), xres, in)
				!= xres)
			{
				fprintf(stderr, "Warning: loadImageaspgm: "
					"File Short: %s\n", filename);
				break;
			}
		}
	}
	else 
	{
		row = g_new(guchar, 3*xres);
		for(j = 0; j < yres; j++)
		{
			if(fread(row, sizeof(unsigned char), 3*xres,in) 
				!= 3*xres)
			{
				fprintf(stderr, "Warning: loadImageaspgm:"
					" File Short: %s\n", filename);
				break;
			}
			k = 0;
			for(i = 0; i < xres; i++)
			{
				c = 71*row[k++];
				c += 150*row[k++];
				c += 28*row[k++];
				img->data[j][i] = c>>8;
			}
		}
		g_free(row);
	}

	fclose(in);

	return img;
}

Image loadImageasweightedpgm(const char *filename, float rw, float gw, float bw)
{
	Image img;
	FILE *in;
	int i, j, k, xres, yres, maxgrey;
	char line[100];
	int mode;
	float wc;
	guchar *row;

	in = fopen(filename, "r");
	if(!in)
	{
		fprintf(stderr, "loadImageaspgm: File not found: %s\n",
			filename);
		return 0;
	}
	if(readpgmline(in, line, filename)) return 0;
	if(line[0] != 'P' || (line[1] != '5' && line[1] != '6'))
	{
		fprintf(stderr, "Unsupported file: %s\n", filename);
		fclose(in);
		return 0;
	}
	if(line[1] == '5') mode = 0;
	else mode = 1;
	if(readpgmline(in, line, filename)) return 0;
	sscanf(line, "%d%d", &xres, &yres);
	if(xres <= 0 || yres <= 0)
	{
		fprintf(stderr, "loadImageaspgm: Bad file: %s\n", filename);
		fclose(in);
		return 0;
	}
	if(readpgmline(in, line, filename)) return 0;
	sscanf(line, "%d", &maxgrey);

	img = newImage(xres, yres);
	if(mode == 0) 
	{
		for(j = 0; j < yres; j++) 
		{
			if(fread(img->data[j], sizeof(unsigned char), xres, in)
				!= xres)
			{
				fprintf(stderr, "Warning: loadImageaspgm: "
					"File Short: %s\n", filename);
				break;
			}
		}
	}
	else 
	{
		row = g_new(guchar, 3*xres);
		for(j = 0; j < yres; j++)
		{
			if(fread(row, sizeof(unsigned char), 3*xres,in) 
				!= 3*xres)
			{
				fprintf(stderr, "Warning: loadImageaspgm:"
					" File Short: %s\n", filename);
				break;
			}
			k = 0;
			for(i = 0; i < xres; i++)
			{
				wc = rw*row[k++];
				wc += gw*row[k++];
				wc += bw*row[k++];
				img->data[j][i] = wc/(rw+gw+bw);
			}
		}
		g_free(row);
	}

	fclose(in);

	return img;
}

/* 0 on success, non-0 otherwise */
int saveImageaspgm(const Image p, const char *filename)
{
	int i;
	FILE *out;

	if(!p) 
	{
		fprintf(stderr, "saveImageaspgm: Null Image.  Dest: %s\n",
			filename);
		return -1;
	}

	out = fopen(filename, "w");
	if(!out)
	{
		fprintf(stderr, "saveImageaspgm: Cannot open %s for write\n",
			filename);
		return -1;
	}
	fprintf(out, "P5\n%d %d\n255\n", p->xres, p->yres);
	for(i = 0; i < p->yres; i++) 
		fwrite(p->data[i], sizeof(unsigned char), p->xres, out);
	fclose(out);
	return 0;
}

Image dupImage(const Image p)
{
	Image img;
	int i;
	
	if(!p) return 0;
	img = newImage(p->xres, p->yres);
	for(i = 0; i < p->yres; i++) 
		img->data[i] = g_memdup(p->data[i], p->xres);
	return img;
}

void blankImage(Image p, unsigned char value)
{
	int i;
	if(!p) return;
	for(i = 0; i < p->yres; i++)
		memset(p->data[i], value, p->xres);
}

void Imagedrawhline(Image p, int x0, int x1, int y, unsigned char value)
{
	if(!p) return;
	memset(p->data[y]+x0, value, x1-x0+1);
}

void stretchImage(Image p)
{
}

Image scaleImage(const Image p, int newxres, int newyres)
{
	return 0;
}

void applyfunctoImage(Image p, unsigned char (* func)(unsigned char v))
{
	int i, j;
	for(j = 0; j < p->yres; j++) for(i = 0; i < p->xres; i++)
		p->data[j][i] = func(p->data[j][i]);
}

Image halveImage(Image p)
{
	Image q;
	int i, j;
	
	q = newImage(p->xres/2, p->yres/2);
	for(j = 0; j < q->yres; j++) for(i = 0; i < q->xres; i++)
		q->data[j][i] = (p->data[2*j][2*i] + p->data[2*j][2*i+1] +
			p->data[2*j+1][2*i] + p->data[2*j+1][2*i+1]) >> 2;
	return q;
}
