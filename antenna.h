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
#ifndef __ANTENNA_H__
#define __ANTENNA_H__

#include "vector.h"
#include "vecarray.h"
#include "keyvalue.h"

/*
 * Antenna parameters
 */
typedef struct
{
	double sub_h;		/* height of subreflector (on axis) */
	double feed[3];		/* position of the feed */
	double feeddir[3];	/* unit vector pointing along feed */
	double radius;		/* antenna radius (m) */
	double K;		
	double deltar;		
	double zedge;		/* height at the edge of the dish */
	double bestparabola;	/* best fit parabola quadratic coef */
	double ftaper;		/* taper of feed */
	double thmax;		/* maximum angle of feed */
	double fa2pi;
	double legwidth;
	double legfoot, legfootz;
	double legapex;
	double leggroundscatter;
	double hole_radius;
	double freq, lambda;
	double oversamp;
	double roughness;	/* RMS surface roughness in meters */
	double Trec, Tsky, Tground;
	double elev;		/* elevation in radians.  0 = point at horiz */
	Vector dir;
	Vector hhat, vhat;	/* unit vectors orthogonal to dir */
	Vector z;
	Vector m;
	Vector k;
	Vector pol;	/* 4 component electric field -- ReH, ImH, ReV, ImV */
	Vector E;	/* 6 component electric field -- ReX, ImX, ... */
	Vector feedpattern;
	double feedpatterndelta;
	char *name;
	int gridsize;
} Antenna;

typedef struct
{
	Matrix subrot;		/* 3x3 matrix rotating x,y,z or nx,ny,nz */
	Matrix feedrot;		/* 3x3 matrix rotating x,y,z or nx,ny,nz */
	Vector subshift;	/* 3 length vector */
	Vector feedshift;	/* 3 length vector */
	Vector subrotpoint;	/* 3 vector describing point to rotate sub. */
	double az_offset;	/* azimuth pointing offset (radians) */
	double el_offset;	/* elevation pointing offset (radians) */
	double phase_offset;	/* DC offset in phase (radians) */
	double focus;		/* meters out of focus toward subreflector */
} Pathology;

typedef struct
{
	Vector aper;		/* aperture x, y, z, nx, ny, nz */
	Vector dish;		/* dish x, y, z, nx, ny, nz */
	Vector sub;		/* subreflector x, y, z, nx, ny, nz */
	Vector feed;		/* feed x, y, z */
} Ray;

Antenna *newAntenna(double sub_h, double feed_x, double feed_y, double feed_z,
	double ftaper, double thmax, const char *geomfile);
void deleteAntenna(Antenna *a);
void Antennasetfreq(Antenna *a, double freq);
void Antennasetdir(Antenna *a, const Vector dir);
void alignfeed(Antenna *a, const Pathology *p);
VecArray getfeedbasis(const Antenna *a, const Pathology *p);
VecArray calcsubreflectorrim(const Antenna *a, const Pathology *p, int nrim);
double subpower(const VecArray rim, const Antenna *a);
Vector Efield(const Antenna *a, const Pathology *p, const Vector pol);
void Antennasetpol(Antenna *a, const Vector pol);
int Antennasetfeedpattern(Antenna *a, const char *filename, double scale);
Antenna *newAntennafromKeyValue(struct KeyValue *kv, const char *paramfilename);
void printAntenna(const Antenna *a);
void defocusAntenna(Pathology *P, const Antenna *a);
int dishvalue(const Antenna *a, double r, double *z, double *m);
int subfromdish(const Antenna *a, double x, double y, Vector subpoint);
int dishfromsub(const Antenna *a, double x, double y, Vector dishpoint);

/*
void antennaR(const Antenna *a, double r, double *x, double *z);
*/
Ray *newRay(const Vector sub);
void deleteRay(Ray *ray);
Pathology *newPathology();
Pathology *dupPathology(const Pathology *P);
Pathology *newPathologyfromKeyValue(struct KeyValue *kv);
void printPathology(const Pathology *P);
void deletePathology(Pathology *P);
double dAdOmega(const Antenna *a, const Ray *ray1, const Ray *ray2, 
	const Ray *ray3, const Pathology *p);
double dOmega(const Antenna *a, const Ray *ray1, const Ray *ray2, 
	const Ray *ray3, const Pathology *p);
double Raylen(const Ray *ray);
double feedfunc(const Antenna *a, double theta);
double feedgain(const Antenna *a, const Ray *ray, const Pathology *p);
double feedpower(const Antenna *a);
double feedbackpower(const Antenna *a);
void Pathologize(Vector sub, const Pathology *p);
void intersectdish(const Antenna *a, const Vector sub, const Vector unitdir, 
	Vector dish, int niter);
void intersectaperture(const Antenna *a, const Vector dish, 
	const Vector unitdir, Vector aper);
Ray *trace(const Antenna *a, double x, double y, const Pathology *p);

Vector tracepol(const Vector E0, const Ray *ray);

int legplanewaveblock(const Antenna *a, double x, double y);
int legplanewaveblock2(const Antenna *a, const Ray *ray);
int legsphericalwaveblock(const Antenna *a, const Ray *ray);

#endif
