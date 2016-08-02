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
#include <stdio.h>
#include <math.h>
#include <glib.h>
#include <string.h>
#include "vector.h"
#include "image.h"
#include "image-vector.h"
#include "vector-fftw.h"
#include "vecarray.h"
#include "polygon.h"
#include "randdist.h"
#include "antenna.h"
#include "illum.h"

const char program[] = "cassbeam";
const char version[] = "1.0";
const char versiondate[] = "08/18/2003";

/*  TODO
 *
 *
 * POST 1.0
 *   Add opacity for gain and Tsys calculations
 *   Compute phase gradient in primary beam (DFT beam cuts)
 *   Oversample blockage calculations?
 *   More accurate diffraction loss?
 *   Feed radius/azimuth instead of x,y
 */

/* given an array of alternating re and im components, mult by exp(i phase) */
void phaseshift(Vector e, double phase)
{
	double c, s;
	Vector E;
	int i, N;

	N = VectorSize(e);
	g_assert(N % 2 == 0);

	c = cos(phase);
	s = sin(phase);

	E = dupVector(e);

	for(i = 0; i < N; i+=2)
	{
		e[i  ] = c*E[i] - s*E[i+1];
		e[i+1] = s*E[i] + c*E[i+1];
	}

	deleteVector(E);
}

/* Returns 4 matrices : Re[Ex], Im[Ex], Re[Ey], Im[Ey] */
Illum *cassillum(const Antenna *a, const Pathology *p, 
	const struct KeyValue *kv)
{
	int i, j, ii, jj, N, n, iter, niter=7;
	double dx, dA, R2, H2, eps;
	Ray *ray, *rayx, *rayy, *bray;
	double x, y, r2, x1, y1, bx, by;
	double dx1, dy1, dx2, dy2;
	double dO, dP;
	double amp, L, L0;
	Vector E, Er, El, Pr, Pl;
	Illum *I;
	double Pforward, Pbackward, Pspill, Pground, Pleg=0.0, Phole=0.0;
	double Psub=0.0, Ppri=0.0, Ptot;
	double phase, cp, sp;
	double subperimeter;
	const int nrim=120;
	VecArray subrim;  /* columns are : X, Y, Z, theta, rho, dtheta */
	VecArray B;

	/* compute total power coming out of the feed */
	Pforward = feedpower(a);
	Pbackward = feedbackpower(a);
	Ptot = Pforward + Pbackward;

	/* compute rim of subreflector and the total power hitting it */
	subrim = calcsubreflectorrim(a, p, nrim);
	Psub = subpower(subrim, a);
	subperimeter = polygonperimeter(subrim[0], subrim[1]);
	
	/* compute image size */
	if(a->gridsize > 0) n = a->gridsize/2;
	else n = 2.0*a->oversamp*a->radius/a->lambda + 0.999;
	if(n%2 == 1) n++;
	N = 2*n;
	
	dx = a->radius/n;
	dA = dx*dx;
	eps = dx/9.237;
	R2 = a->radius*a->radius;
	H2 = a->hole_radius*a->hole_radius;

	/* compute central ray pathlength */
	ray = trace(a, 0.0, 0.00001, p);
	L0 = Raylen(ray);
	deleteRay(ray);

	/* compute polarization vectors */
	Pr = newVector(4);
	Pl = newVector(4);
	Pr[0] = 1.0/M_SQRT2; Pr[1] = 0.0; Pr[2] = 0.0; Pr[3] = 1.0/M_SQRT2;
	Pl[0] = 1.0/M_SQRT2; Pl[1] = 0.0; Pl[2] = 0.0; Pl[3] = -1.0/M_SQRT2;
	Er = Efield(a, p, Pr);
	El = Efield(a, p, Pl);
	deleteVector(Pr);
	deleteVector(Pl);

	B = getfeedbasis(a, p);

	/* Initialize the illumination structure */
	I = newIllum(N);
	setIllumfromKeyValue(I, kv);
	I->dA = dA;
	I->a0 = a->radius*a->radius*M_PI;
	I->lambda = a->lambda;
	I->Tsky = a->Tsky;
	I->Tground = a->Tground;
	I->Trec = a->Trec;
	I->Xangle = atan2(B[0][1], B[0][0]);
	I->gridsize = a->gridsize;
	
	phaseshift(Er, -I->Xangle);
	phaseshift(El, I->Xangle);
	I->Xangle = 0.0;
	
	/* compute field on primary */
	for(j = -n; j < n; j++)
	{
		y = (j + 0.5)*dx;
		jj = j + n;
		for(i = -n; i < n; i++) 
		{
			x = (i + 0.5)*dx;
			ii = i + n;
			ray = rayx = rayy = bray = 0;
			
			r2 = x*x + y*y;
			if(r2 > R2) continue;

			x1 = x;
			y1 = y;

			/* Iterate to find starting point */
			for(iter = 0; iter < niter; iter++)
			{
				ray = trace(a, x1, y1, p);
				if(!ray) goto nextpoint;
				x1 += (x - ray->aper[0]);
				y1 += (y - ray->aper[1]);
				deleteRay(ray);
				ray = 0;
				if(x1*x1 + y1*y1 > R2) goto nextpoint;
			}
	
			if(y < 0) rayy = trace(a, x1, y1+eps, p);
			else rayy = trace(a, x1, y1-eps, p);
			
			if(x < 0) rayx = trace(a, x1+eps, y1, p);
			else rayx = trace(a, x1-eps, y1, p);
			
			ray = trace(a, x1, y1, p);
			
			if(ray == 0 || rayx == 0 || rayy == 0)
				goto nextpoint;

			dx1 = rayx->aper[0]-ray->aper[0];
			dy1 = rayx->aper[1]-ray->aper[1];
			dx2 = rayy->aper[0]-ray->aper[0];
			dy2 = rayy->aper[1]-ray->aper[1];

			/* To make blockage statistically correct... */
			bx = x + 0.5*rand_pm_one()*(dx1/eps)*dx;
			by = y + 0.5*rand_pm_one()*(dy2/eps)*dx;
			bray = trace(a, bx, by, p);

			dA = 0.5*fabs(dx1*dy2 - dx2*dy1);

			dO = (dOmega(a, rayx, rayy, ray, p)/dA)*dx*dx;
			dP = dO*feedgain(a, ray, p);

			/* Check for blockage */
			if(bx*bx+by*by < H2) /* hole in the middle */
			{
				I->B[jj][ii] = 1.0;
				Phole += dP;
			}
			else if(legplanewaveblock2(a, bray) ||   /* legs */
				legsphericalwaveblock(a, bray))
			{
				I->B[jj][ii] = 1.0;
				Pleg += dP;
			}
			
#if 0
			if(polygonside(subrim[0], subrim[1], 
				ray->sub[0], ray->sub[1]) < 0)
			{
				Po += dP;
				printf(".");
				goto nextpoint;
			}
#endif

			Ppri += dP;

			/* Fill in Scalar field matrices */
			amp = sqrt(dP);
			L = Raylen(ray);
			phase = 2.0*M_PI*(L-L0)/a->lambda;
			sp = sin(phase);
			cp = cos(phase);
			I->A[jj][ii] = amp;
			I->P[jj][ii] = phase;
			
			/* Fill in Vector field matrices */
			E = tracepol(Er, ray);
			I->E[RXR][jj][ii] = amp*(E[0]*cp - E[1]*sp);
			I->E[IXR][jj][ii] = amp*(E[0]*sp + E[1]*cp);
			I->E[RYR][jj][ii] = amp*(E[2]*cp - E[3]*sp);
			I->E[IYR][jj][ii] = amp*(E[2]*sp + E[3]*cp);
			I->E[RZR][jj][ii] = amp*(E[4]*cp - E[5]*sp);
			I->E[IZR][jj][ii] = amp*(E[4]*sp + E[5]*cp);
			deleteVector(E);
			E = tracepol(El, ray);
			I->E[RXL][jj][ii] = amp*(E[0]*cp - E[1]*sp);
			I->E[IXL][jj][ii] = amp*(E[0]*sp + E[1]*cp);
			I->E[RYL][jj][ii] = amp*(E[2]*cp - E[3]*sp);
			I->E[IYL][jj][ii] = amp*(E[2]*sp + E[3]*cp);
			I->E[RZL][jj][ii] = amp*(E[4]*cp - E[5]*sp);
			I->E[IZL][jj][ii] = amp*(E[4]*sp + E[5]*cp);
			deleteVector(E);

		nextpoint:
			if(ray)  deleteRay(ray);
			if(rayx) deleteRay(rayx);
			if(rayy) deleteRay(rayy);
			if(bray) deleteRay(bray);
		}
	}

	deleteVecArrayandVectors(B);
	deleteVecArrayandVectors(subrim);
	deleteVector(Er);
	deleteVector(El);

	/* remove phase gradient, populate pointing fields */
	dephaseIllum(I);
	/* again, with incremental inprovements */
	dephaseIllum(I);

	Pspill = Ptot - Ppri;

	Pground = Pspill - (Ptot-Psub) + a->leggroundscatter*Pleg;

	I->groundfrac = Pground/Ptot;
	I->spilleff = 1.0-Pspill/Ptot;
	I->subspilleff = Psub/Ptot;
	I->prispilleff = I->spilleff/I->subspilleff;
	x = 4.0*M_PI*a->roughness/a->lambda;
	I->surfeff = exp(-x*x);
	/* Emperical -- FIXME */
	if(I->diffeff == 0.0)
		I->diffeff = 1.0-5.5*a->lambda/subperimeter;

	calcIllumparams(I);

	return I;
}

void savecassdata(const Antenna *a, const Illum *I)
{
	char filename[1000];

	if(I->saveapertureimages)
	{
		sprintf(filename, "%s.illumamp.pgm", I->prefix);
		saveMatrixaspgm(I->A, filename);
		sprintf(filename, "%s.illumphase.pgm", I->prefix);
		saveMatrixaspgm(I->P, filename);
		sprintf(filename, "%s.illumblock.pgm", I->prefix);
		saveMatrixaspgm(I->B, filename);
	}
}

int makecassbeam(struct KeyValue *kv, const char *paramfile)
{
	Pathology *p;
	Antenna *a;
	Illum *I;
	char filename[1000];
	
	a = newAntennafromKeyValue(kv, paramfile);
	g_assert(a);
	p = newPathologyfromKeyValue(kv);
	g_assert(p);

	alignfeed(a, p);

	printAntenna(a);
	printPathology(p);

	if(p->focus != 0.0) defocusAntenna(p, a);

	I = cassillum(a, p, kv);
	printIllum(I);
	calcIllumpolparams(I);
	
	savecassdata(a, I);
	
	KeyValueupdateIllum(kv, I);

	if(I->saveparams)
	{
		sprintf(filename, "%s.params", I->prefix);
		saveKeyValue(kv, filename);
	}
	
	deleteIllum(I);

	deletePathology(p);
	deleteAntenna(a);

	return 1;
}

int main(int argc, char **argv)
{
	struct KeyValue *kv;

	int i, j;
	char k[100], v[100];

	if(argc < 2)
	{
		fprintf(stderr, "%s version %s\n\n", program, version);
		fprintf(stderr, "Walter Brisken, NRAO, %s\n\n", versiondate);
		fprintf(stderr, "usage : %s <filename> [options]\n", argv[0]);
		return 0;
	}

	kv = loadKeyValue(argv[1]);
	if(!kv) return 0;
	
	if(argc > 2) for(i = 2; i < argc; i++)
	{
		if(argv[i][0] == '=') fprintf(stderr,
			"Warning -- option begins with =\n");
		for(j = 0; argv[i][j] != 0; j++) if(argv[i][j] == '=')
			argv[i][j] = ' ';
		if(sscanf(argv[i], "%s %s", k, v) != 2) continue;
		
		j = KeyValuekeyindex(kv, k);
		if(j < 0) KeyValueaddparm(kv, k, v);
		else
		{
			g_free(kv->value[j]);
			kv->value[j] = g_strdup(v);
		}
	}

	/* force some default values */
	if(KeyValuekeyindex(kv, "out") < 0)
		KeyValueupdateparm(kv, "out", program);
	if(KeyValuekeyindex(kv, "compute") < 0)
		KeyValueupdateparm(kv, "compute", "all");
	KeyValueupdateparm(kv, "program", program);
	KeyValueupdateparm(kv, "version", version);

	makecassbeam(kv, argv[1]);

	deleteKeyValue(kv);

	printf("\n%s has finished\n", program);

	return 0;
}
				
