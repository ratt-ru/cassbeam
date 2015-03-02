#include <stdio.h>
#include <math.h>
#include <glib.h>
#include <string.h>
#include "constants.h"
#include "antenna.h"
// #include "holog.h"
// #include "vecsurface.h"
// #include "vecdiffract.h"

Antenna *newAntenna(double sub_h, double feed_x, double feed_y, double feed_z,
	double ftaper, double thmax, const char *geomfile)
{
	Antenna *a;
	double z;
	VecArray arr;
	int n;
	double d;

	arr = VecArrayfromfile(geomfile, 3);
	if(!arr) 
	{
		fprintf(stderr, "File %s not found\n", geomfile);
		return 0;
	}

	n = VecArrayVectorSize(arr);
	g_assert(n>2);

	z = sub_h-feed_z;

	a = g_new(Antenna, 1);
	a->sub_h = sub_h;
	a->feed[0] = feed_x;
	a->feed[1] = feed_y;
	a->feed[2] = feed_z;
	d = sqrt(feed_x*feed_x + feed_y*feed_y + z*z);
	if(z > 0.0)
	{
		a->K = sub_h + d;
		a->feeddir[0] = -feed_x/d;
		a->feeddir[1] = -feed_y/d;
		a->feeddir[2] = (sub_h-feed_z)/d;
	}
	else
	{
		a->K = sqrt(feed_x*feed_x + feed_y*feed_y + feed_z*feed_z);
		a->feeddir[0] = -feed_x/d;
		a->feeddir[1] = -feed_y/d;
		a->feeddir[2] = (sub_h-feed_z)/d;
	}
	a->radius = arr[0][n-1];
	a->zedge = arr[1][n-1];
	a->deltar = a->radius/(float)(n-1);
	a->bestparabola = a->zedge/(a->radius*a->radius);
	a->z = arr[1];
	a->m = arr[2];
	g_assert(a->m[0] == 0.0);
	if(a->ftaper < 0.0)
		fprintf(stderr, "Warning : changing sign of feedtaper\n");
	a->ftaper = fabs(ftaper);
	a->thmax = thmax;
	a->fa2pi = 2.0*M_PI*sqrt(ftaper)*0.1874/sin(thmax*M_PI/180.0);
	a->legwidth = 0.0;
	a->legfoot = a->radius/2.0;
	a->legapex = sub_h*1.2;
	a->hole_radius = 0.0; 
	a->roughness = 0.0;
	a->dir = newVector(3);
	a->dir[0] = a->dir[1] = 0.0;
	a->dir[2] = 1.0;
	a->hhat = newVector(3);
	a->vhat = newVector(3);
	a->k = newVector(3);
	a->name = g_strdup("unnamed");
	zeroVector(a->k);
	/* default to no polarization state */
	a->pol = newVector(4);
	zeroVector(a->pol);
	a->E = newVector(6);
	zeroVector(a->E);
	Antennasetfreq(a, 1.0);
	Antennasetdir(a, 0);  /* compute hhat and vhat */
	a->Trec = 50.0;
	a->Tground = 290.0;
	a->oversamp = 1.5;
	a->elev = M_PI/2.0;
	a->feedpattern = 0;
	a->feedpatterndelta = 0.0;
	a->gridsize = 0;
	a->leggroundscatter = 0.2;
	dishvalue(a, a->legfoot, &a->legfootz, 0);

	deleteVector(arr[0]);
	deleteVecArray(arr);

	return a;
}

void deleteAntenna(Antenna *a)
{
	if(!a) return;

	if(a->z) deleteVector(a->z);
	if(a->m) deleteVector(a->m);
	if(a->k) deleteVector(a->k);
	if(a->dir) deleteVector(a->dir);
	if(a->hhat) deleteVector(a->hhat);
	if(a->vhat) deleteVector(a->vhat);
	if(a->pol) deleteVector(a->pol);
	if(a->E) deleteVector(a->E);
	
	g_free(a);
}

void Antennasetfreq(Antenna *a, double freq)
{
	int i;
	
	g_assert(freq > 0.0);

	a->freq = freq;
	a->lambda = NS_METER/freq;
	if(freq < 1.0) a->Tsky = 3.0*pow(freq, -2.5);
	else a->Tsky = 3.0;
	for(i = 0; i < 3; i++) a->k[i] = -2.0*M_PI*a->dir[i]/a->lambda;
}

void Antennasetdir(Antenna *a, const Vector dir)
{
	double hmag;
	int i, l;

	if(dir)
	{
		l = VectorSize(dir);
		g_assert(l == 1 || l == 3);
		
		if(l == 1)
		{
			a->dir[0] = -sin(dir[0]*M_PI/180.0);
			a->dir[1] = 0.0;
			a->dir[2] = cos(dir[0]*M_PI/180.0);
		}
		else
		{
			deleteVector(a->dir);
			a->dir = newunitVector(dir);
		}

		if(a->dir[0] == 0.0 && a->dir[1] == 0.0)
		{
			a->hhat[0] = 1.0;
			a->hhat[1] = a->hhat[2] = 0.0;
			a->vhat[1] = 1.0;
			a->vhat[0] = a->vhat[2] = 0.0;
		}
		else
		{
			a->hhat[0] = a->dir[1];
			a->hhat[1] = -a->dir[0];
			a->hhat[2] = 0.0;
			hmag = sqrt(a->hhat[0]*a->hhat[0]
				  + a->hhat[1]*a->hhat[1]);
			a->hhat[0] /= hmag;
			a->hhat[1] /= hmag;

			a->vhat[0] = a->hhat[1]*a->dir[2] 
				   - a->hhat[2]*a->dir[1];
			a->vhat[1] = a->hhat[2]*a->dir[0] 
				   - a->hhat[0]*a->dir[2];
			a->vhat[2] = a->hhat[0]*a->dir[1] 
				   - a->hhat[1]*a->dir[0];
		}
	}
	for(i = 0; i < 3; i++) a->k[i] = -2.0*M_PI*a->dir[i]/a->lambda;
	Antennasetpol(a, 0); /* recompute polarization */
}

void Antennasetpol(Antenna *a, const Vector pol)
{
	if(pol)
	{
		g_assert(VectorSize(pol) == 4);
		copytoVector(a->pol, pol);
		normalizeVector(a->pol);
	}

	a->E[0] = a->hhat[0]*a->pol[0] + a->vhat[0]*a->pol[2];
	a->E[1] = a->hhat[0]*a->pol[1] + a->vhat[0]*a->pol[3];
	a->E[2] = a->hhat[1]*a->pol[0] + a->vhat[1]*a->pol[2];
	a->E[3] = a->hhat[1]*a->pol[1] + a->vhat[1]*a->pol[3];
	a->E[4] = a->hhat[2]*a->pol[0] + a->vhat[2]*a->pol[2];
	a->E[5] = a->hhat[2]*a->pol[1] + a->vhat[2]*a->pol[3];
}

/* sets feeddir after pathology is considered */
void alignfeed(Antenna *a, const Pathology *p)
{
	int i, j;
	double f[3], s0[3], s[3], d[3], m=0.0;

	for(i = 0; i < 3; i++) f[i] = a->feed[i] + p->feedshift[i];
	for(i = 0; i < 3; i++) s0[i] = -p->subrotpoint[i];
	s0[2] += a->sub_h;
	for(i = 0; i < 3; i++) 
	{
		s[i] = 0.0;
		for(j = 0; j < 3; j++) 
			s[i] += p->subrot[i][j]*s0[j];
		s[i] += p->subrotpoint[i] + p->subshift[i];
		d[i] = s[i]-f[i];
		m += d[i]*d[i];
	}
	m = sqrt(m);
	for(i = 0; i < 3; i++) a->feeddir[i] = d[i]/m;
}

VecArray getfeedbasis(const Antenna *a, const Pathology *p)
{
	int i, j;
	VecArray B;
	Vector dir, vhat, hhat;

	B = newpopulatedVecArray(3, 3);

	hhat = B[0];
	vhat = B[1];
	dir = B[2];

	if(p == 0) for(i = 0; i < 3; i++) dir[i] = a->feeddir[i];
	else 
	{
		zeroVector(dir);
		for(j = 0; j < 3; j++) for(i = 0; i < 3; i++)
			dir[j] += p->feedrot[j][i]*a->feeddir[i];
	}

	if(dir[0] == 0.0 && dir[1] == 0.0)
	{
		vhat[0] = 1.0;
		vhat[1] = vhat[2] = 0.0;
		hhat[1] = 1.0;
		hhat[0] = hhat[2] = 0.0;
	}
	else
	{
		vhat[0] = dir[1];
		vhat[1] = -dir[0];
		vhat[2] = 0.0;
		normalizeVector(vhat);

		hhat[0] = vhat[1]*dir[2] - vhat[2]*dir[1];
		hhat[1] = vhat[2]*dir[0] - vhat[0]*dir[2];
		hhat[2] = vhat[0]*dir[1] - vhat[1]*dir[0];
	}

	return B;
}

/* Returns vecarray.  each vector is nrim long.  The columns are:
 *
 *   [0] x coordinate of rim
 *   [1] y coordinate of rim
 *   [2] z coordinate of rim
 *   [3] theta (azimuth angle) of point as seen by the feed
 *   [4] phi (boresight angle) of point as seen by the feed
 *   [5] dtheta, useful for integrations.
 */
VecArray calcsubreflectorrim(const Antenna *a, const Pathology *p, int nrim)
{
	VecArray subrim, feedbasis;
	Vector s1, subr, h, v, feedaxis;
	double Reff, theta;
	int i, j;
	
	feedbasis = getfeedbasis(a, p);
	h = feedbasis[0];
	v = feedbasis[1];
	feedaxis = feedbasis[2];
	
	subrim = newpopulatedVecArray(6, nrim);
	s1 = newVector(6);
	subr = newVector(3);
	
	/* Reff is such that the polygon circumcsribes the primary */
	Reff = a->radius/cos(M_PI/nrim);

	for(i = 0; i < nrim; i++)
	{
		theta = 2.0*M_PI*i/nrim;
		subfromdish(a, Reff*cos(theta), Reff*sin(theta), s1);
		Pathologize(s1, p);
		subrim[0][i] = s1[0];
		subrim[1][i] = s1[1];
		subrim[2][i] = s1[2];
		for(j = 0; j < 3; j++)
			subr[j] = s1[j] - (a->feed[j]+p->feedshift[j]);
		normalizeVector(subr);

		subrim[3][i] = atan2(dotVectors(subr, v), dotVectors(subr, h));
		subrim[4][i] = acos(dotVectors(subr, feedaxis));
	}

	deleteVecArrayandVectors(feedbasis);
	deleteVector(subr);
	deleteVector(s1);
	
	/* compute dtheta */
	subrim[5][0] = subrim[3][1]-subrim[3][nrim-1];
	subrim[5][nrim-1] = subrim[3][0]-subrim[3][nrim-2];
	for(i = 1; i < nrim-1; i++)
		subrim[5][i] = subrim[3][i+1]-subrim[3][i-1];
	for(i = 0; i < nrim; i++)
	{
		if(subrim[5][i] < -M_PI) subrim[5][i] += 2.0*M_PI;
		if(subrim[5][i] >  M_PI) subrim[5][i] -= 2.0*M_PI;
		subrim[5][i] = 0.5*fabs(subrim[5][i]);
	}

	return subrim;
}

/* rim is derived from calcsubreflectorrim */
double subpower(const VecArray rim, const Antenna *a)
{
	const int ntheta = 2000;
	int i, t, N;
	double theta, dtheta, power=0.0;
	double domega;

	g_assert(rim);
	g_assert(VecArraySize(rim) == 6);

	N = VecArrayVectorSize(rim);

	g_assert(N > 2);

	dtheta = Vectormax(rim[4])/(ntheta-0.5);

	for(t = 0; t < ntheta; t++)
	{
		theta = (t+0.5)*dtheta;
		domega = 0.0;
		for(i = 0; i < N; i++) 
			if(rim[4][i] > theta) domega += rim[5][i];

		domega*=sin(theta)*dtheta;

		power += domega*feedfunc(a, theta);
	}

	return power;
}

Vector Efield(const Antenna *a, const Pathology *p, const Vector pol)
{
	VecArray B;
	Vector hhat, vhat, E;

	B = getfeedbasis(a, p);
	hhat = B[0];
	vhat = B[1];

	E = newVector(6);
	
	E[0] = hhat[0]*pol[0] + vhat[0]*pol[2];
	E[1] = hhat[0]*pol[1] + vhat[0]*pol[3];
	E[2] = hhat[1]*pol[0] + vhat[1]*pol[2];
	E[3] = hhat[1]*pol[1] + vhat[1]*pol[3];
	E[4] = hhat[2]*pol[0] + vhat[2]*pol[2];
	E[5] = hhat[2]*pol[1] + vhat[2]*pol[3];

	deleteVecArrayandVectors(B);

	return E;
}

int Antennasetfeedpattern(Antenna *a, const char *filename, double scale)
{
	int i, N, Nmax;
	double x, delta;
	VecArray pat;
	
	a->feedpatterndelta = 0.0;
	if(a->feedpattern) deleteVector(a->feedpattern);

	if(filename == 0) return 1;

	pat = VecArrayfromfile(filename, 2);

	if(!pat) return 0;
	N = VectorSize(pat[0]);
	g_assert(N > 2);
	g_assert(pat[0][0] == 0.0);
	
	delta = pat[0][1];
	g_assert(delta > 0.0);
	for(i = 2; i < N; i++) 
	{
		x = pat[0][i]-pat[0][i-1]-delta;
		g_assert(fabs(x) < delta/10000.0);
	}

	/* convert to radians */
	delta *= M_PI/180.0;

	/* and scale it */
	if(scale > 0.0) delta *= scale;

	/* Do we need to truncate the pattern? */
	Nmax = M_PI/delta;
	if(N > Nmax)
	{
		a->feedpattern = newVector(Nmax);
		for(i = 0; i < Nmax; i++) 
			a->feedpattern[i] = fabs(pat[1][i]);
		deleteVector(pat[1]);
	}
	else a->feedpattern = pat[1];

	a->feedpatterndelta = delta;
	deleteVector(pat[0]);
	deleteVecArray(pat);

	return 1;
}

Antenna *newAntennafromKeyValue(struct KeyValue *kv, const char *paramfilename)
{
	Antenna *a;
	Vector V;
	double v;
	char paramfile[100];
	double sub_h, feed_x, feed_y, feed_z, thmax, ftaper;
	const char *geomfile, *feedfile;
	const char *name;
	int i;

	if(paramfilename) strncpy(paramfile, paramfilename, 99);
	else strncpy(paramfile, "(input file)", 99);
	paramfile[99] = 0;
	
	sub_h = getKeyValuedouble(kv, "sub_h");
	if(sub_h == KV_FLOATERR)
	{
		fprintf(stderr, "sub_h not defined in %s\n", paramfile);
		return 0;
	}
	
	feed_x = getKeyValuedouble(kv, "feed_x");
	if(feed_x == KV_FLOATERR)
	{
		fprintf(stderr, "feed_x not defined in %s\n", paramfile);
		return 0;
	}
	
	feed_y = getKeyValuedouble(kv, "feed_y");
	if(feed_y == KV_FLOATERR)
	{
		fprintf(stderr, "feed_y not defined in %s\n", paramfile);
		return 0;
	}
	
	feed_z = getKeyValuedouble(kv, "feed_z");
	if(feed_z == KV_FLOATERR)
	{
		fprintf(stderr, "feed_z not defined in %s\n", paramfile);
		return 0;
	}

	feedfile = getKeyValuestring(kv, "feedpattern");
	
	thmax = getKeyValuedouble(kv, "feedthetamax");
	if(thmax == KV_FLOATERR && !feedfile)
	{
		fprintf(stderr, "feedthetamax not defined in %s\n", paramfile);
		return 0;
	}
	
	ftaper = getKeyValuedouble(kv, "feedtaper");
	if(ftaper == KV_FLOATERR && !feedfile)
	{
		fprintf(stderr, "feedthetamax not defined in %s\n", paramfile);
		return 0;
	}
	
	geomfile = getKeyValuestring(kv, "geom");
	if(geomfile == 0)
	{
		fprintf(stderr, "geom not defined in %s\n", paramfile);
		return 0;
	}

	a = newAntenna(sub_h, feed_x, feed_y, feed_z, ftaper, thmax, geomfile);
	g_assert(a);
	
	name = getKeyValuestring(kv, "name");
	if(name)
	{
		g_free(a->name);
		a->name = g_strdup(name);
	}
	
	/* feed pattern file is two column text file containing 
	 * angle (in degrees) and power (in dBi) 
	 */
	if(feedfile != 0)
	{
		double scale;
		scale = getKeyValuedouble(kv, "feedpatternscale");
		if(!Antennasetfeedpattern(a, feedfile, scale)) 
		{
			deleteAntenna(a);
			fprintf(stderr, "Problem with feed file <%s>\n",
				feedfile);
			return 0;
		}
	}
	
	i = getKeyValueint(kv, "gridsize");
	if(i > 0) 
	{
		if(i < 32) a->gridsize = 32;
		else a->gridsize = i;
	}
	
	v = getKeyValuedouble(kv, "freq");
	if(v > 0.0) Antennasetfreq(a, v);
	
	v = getKeyValuedouble(kv, "OS");
	if(v > 0.0) a->oversamp = v;
	
	v = getKeyValuedouble(kv, "oversamp");
	if(v > 0.0) a->oversamp = v;
	
	v = getKeyValuedouble(kv, "legwidth");
	if(v != KV_FLOATERR) a->legwidth = v;

	v = getKeyValuedouble(kv, "legfoot");
	if(v != KV_FLOATERR) 
	{
		a->legfoot = v;
		dishvalue(a, a->legfoot, &a->legfootz, 0);
	}

	v = getKeyValuedouble(kv, "legapex");
	if(v != KV_FLOATERR) a->legapex = v;

	v = getKeyValuedouble(kv, "hole_radius");
	if(v != KV_FLOATERR) a->hole_radius = v;

	v = getKeyValuedouble(kv, "Tground");
	if(v != KV_FLOATERR) a->Tground = v;

	v = getKeyValuedouble(kv, "Trec");
	if(v != KV_FLOATERR) a->Trec = v;

	v = getKeyValuedouble(kv, "elev");
	if(v != KV_FLOATERR) a->elev = v*M_PI/180.0;

	v = getKeyValuedouble(kv, "roughness");
	if(v != KV_FLOATERR) a->roughness = v;

	v = getKeyValuedouble(kv, "leggroundscatter");
	if(v != KV_FLOATERR) a->leggroundscatter = v;

	V = getKeyValueVector(kv, "point");
	if(V)
	{
		Antennasetdir(a, V);
		deleteVector(V);
	}
	else
	{
		V = newVector(1);
		V[0] = 0.0;
		Antennasetdir(a, V);
		deleteVector(V);
	}

	V = getKeyValueVector(kv, "pol");
	if(V)
	{
		Antennasetpol(a, V);
		deleteVector(V);
	}
	else
	{
		V = newVector(4);
		zeroVector(V);
		V[0] = 1.0;	/* liner pol */
		Antennasetpol(a, V);
		deleteVector(V);
	}

	return a;
}

int dishvalue(const Antenna *a, double r, double *z, double *m)
{
	double ma, mb, mc, zav, A, B, C, D;
	double x, d, dd;
	double s = 1.0;
	int n;

	g_assert(a);

	if(r == 0)
	{
		*z = a->z[0];
		*m = 0.0;
		return 1;
	}
	
	if(r < 0) 
	{
		s = -1.0;
		r = -r;
	}
	d = a->deltar;
	dd = d*d;
	
	n = floor(r/d + 0.5);	/* the middle point */
	if(n > VectorSize(a->m)-2) n = VectorSize(a->m)-2;

	x = r - n*d;

	if(n == 0)
	{
		mc = a->m[1];
		ma = -mc;
		mb = 0.0;
		zav = 2.0*a->z[1] + a->z[0];
	}
	else
	{
		ma = a->m[n-1];
		mb = a->m[n];
		mc = a->m[n+1];
		zav = a->z[n-1] + a->z[n] + a->z[n+1];
	}
	
	A = mb;
	B = 0.5*(mc - ma)/d;
	C = 0.5*(mc - 2.0*mb + ma)/dd;

	D = (zav - B*dd)/3.0;

	if(m) *m = s*(A + B*x + C*x*x);
	if(z) *z = s*(D + A*x + B*x*x/2.0 + C*x*x*x/3.0);
	
	return 1;
}

/* normalizes a "vector" of 3 doubles in the vector sense */
static inline void norm3(double *v)
{
	double s;
	s = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	v[0] /= s;
	v[1] /= s;
	v[2] /= s;
}

/* Returns position of subreflector piece (x, y, z) and
 * its normal (u, v, w)
 */
int subfromdish(const Antenna *a, double x, double y, Vector subpoint)
{
	double r, z, m, u, v, w;
	double dx, dy, dz, dl, t;
	double n[3], sf[3], sd[3];
	int i;

	r = sqrt(x*x + y*y);

	if(r == 0)
	{
		subpoint[0] = 0.0;
		subpoint[1] = 0.0;
		subpoint[2] = a->sub_h;
	}
	else
	{
		dishvalue(a, r, &z, &m);

		/* Compute reflected unit vector direction */
		m = tan(2.0*atan(m));
		w = 1.0/sqrt(1.0+m*m);
		u = -m*(x/r)*w;
		v = -m*(y/r)*w;
		
		dx = a->feed[0]-x;
		dy = a->feed[1]-y;
		dz = a->feed[2]-z;
		dl = a->K + z;

		t = 0.5*(dx*dx + dy*dy + dz*dz - dl*dl)
		      / (-dl + u*dx + v*dy + w*dz);

		subpoint[0] = x + u*t;
		subpoint[1] = y + v*t;
		subpoint[2] = z + w*t;
	}

	for(i = 0; i < 3; i++) sf[i] = a->feed[i] - subpoint[i];
	sd[0] = x - subpoint[0];
	sd[1] = y - subpoint[1];
	sd[2] = z - subpoint[2];

	norm3(sf);
	norm3(sd);

	for(i = 0; i < 3; i++) n[i] = sd[i]+sf[i];

	norm3(n);

	for(i = 0; i < 3; i++) subpoint[3+i] = n[i];

	return 1;
}

int dishfromsub(const Antenna *a, double x, double y, Vector dishpoint)
{
	VecArray sub;
	double x1, y1, dx, dy, mx, my, r, d;
	const double eps = 0.001;
	int iter, niter=500;

	sub = newpopulatedVecArray(5, 6);

	x1 = x;
	y1 = y;

	for(iter = 0; iter < niter; iter++)
	{
		subfromdish(a, x1, y1, sub[0]);
		subfromdish(a, x1-eps, y1, sub[1]);
		subfromdish(a, x1+eps, y1, sub[2]);
		subfromdish(a, x1, y1-eps, sub[3]);
		subfromdish(a, x1, y1+eps, sub[4]);
		mx = 0.5*(sub[2][0]-sub[1][0])/eps;
		my = 0.5*(sub[4][1]-sub[3][1])/eps;
		dx = (x-sub[0][0])/mx;
		dy = (y-sub[0][1])/my;
		if(fabs(dx) > a->radius/7.0) 
		{
			if(dx < 0) dx = -a->radius/7.0;
			else dx = a->radius/7.0;
		}
		if(fabs(dy) > a->radius/7.0) 
		{
			if(dy < 0) dy = -a->radius/7.0;
			else dy = a->radius/7.0;
		}
		r = sqrt(x1*x1 + y1*y1);
		if(r >= a->radius)
			if(x1*dx + y1*dy > 0.0) return 0;
		x1 += 0.5*dx;
		y1 += 0.5*dy;
		if(fabs(dx) < 0.005*eps && fabs(dy) < 0.005*eps) break;
	}
	if(iter == niter) return 0;

	r = sqrt(x1*x1 + y1*y1);
	dishpoint[0] = x1;
	dishpoint[1] = y1;
//	dishpoint[2] = polyvalue(a->shape, r);
	dishpoint[3] = sub[0][0];
	dishpoint[4] = sub[0][1];
	dishpoint[5] = sub[0][2];
	d = sqrt(1.0+mx*mx+my*my);
	dishpoint[6] = mx/d;
	dishpoint[7] = my/d;
	dishpoint[8] = 1.0/d;
	dishpoint[9] = sub[0][3];
	dishpoint[10] = sub[0][4];
	dishpoint[11] = sub[0][5];

	deleteVecArrayandVectors(sub);

	if(r > a->radius) return 0;
	else return 1;
}

void printAntenna(const Antenna *a)
{
	printf("Antenna: %s  %p\n", a->name, a);
	printf("  freq    = %f GHz  lambda = %f m\n", a->freq, a->lambda);
	printf("  Tsky    = %f K Tground = %f K Trec = %f K\n",
		a->Tsky, a->Tground, a->Trec);
	printf("  dir     = "); printVector(a->dir);
	printf("  feeddir = %f, %f, %f\n", 
		a->feeddir[0], a->feeddir[1], a->feeddir[2]); 
	if(!a->feedpattern) 
	{
		printf("  ftaper  = %f\n", a->ftaper);
		printf("  thmax   = %f\n", a->thmax);
	}
	printf("\n");
}

Ray *newRay(const Vector sub)
{
	Ray *ray;

	ray = g_new(Ray, 1);
	ray->aper = newVector(6);
	ray->dish = newVector(6);
	ray->sub = dupVector(sub);
	ray->feed = newVector(3);

	return ray;
}

void printRay(const Ray *ray)
{
	printf("Ray : ");
	printVector(ray->sub);
	printf("      ");
	printVector(ray->dish);
	printf("      ");
	printVector(ray->aper);
}

void deleteRay(Ray *ray)
{
	g_assert(ray);

	if(ray->aper) deleteVector(ray->aper);
	if(ray->dish) deleteVector(ray->dish);
	if(ray->sub ) deleteVector(ray->sub);
	if(ray->feed) deleteVector(ray->feed);

	g_free(ray);
}

Pathology *newPathology()
{
	Pathology *P;
	int i;

	P = g_new(Pathology, 1);

	P->subrot = newMatrix(3, 3);
	zeroMatrix(P->subrot);
	for(i = 0; i < 3; i++) P->subrot[i][i] = 1.0;

	P->subrotpoint = newVector(3);
	zeroVector(P->subrotpoint);

	P->feedrot = newMatrix(3, 3);
	zeroMatrix(P->feedrot);
	for(i = 0; i < 3; i++) P->feedrot[i][i] = 1.0;

	P->subshift = newVector(3);
	zeroVector(P->subshift);

	P->feedshift = newVector(3);
	zeroVector(P->feedshift);

	P->az_offset = 0.0;
	P->el_offset = 0.0;
	P->phase_offset = 0.0;
	P->focus = 0.0;

	return P;
}

Pathology *dupPathology(const Pathology *P)
{
	Pathology *Q;

	Q = g_new(Pathology, 1);
	Q->subrot = dupMatrix(P->subrot);
	Q->subrotpoint = dupVector(P->subrotpoint);
	Q->feedrot = dupMatrix(P->feedrot);
	Q->subshift = dupVector(P->subshift);
	Q->feedshift = dupVector(P->feedshift);
	Q->az_offset = P->az_offset;
	Q->el_offset = P->el_offset;
	Q->phase_offset = P->phase_offset;
	Q->focus = P->focus;

	return Q;
}

Pathology *newPathologyfromKeyValue(struct KeyValue *kv)
{
	Pathology *P;
	double v;
	double srx=0.0, sry=0.0, srz=0.0, frx=0.0, fry=0.0, frz=0.0;
	Vector V;

	P = newPathology();
	
	v = getKeyValuedouble(kv, "dfeed_x");
	if(v != KV_FLOATERR) P->feedshift[0] = v;
	v = getKeyValuedouble(kv, "dfeed_y");
	if(v != KV_FLOATERR) P->feedshift[1] = v;
	v = getKeyValuedouble(kv, "dfeed_z");
	if(v != KV_FLOATERR) P->feedshift[2] = v;
	v = getKeyValuedouble(kv, "dsub_x");
	if(v != KV_FLOATERR) P->subshift[0] = v;
	v = getKeyValuedouble(kv, "dsub_y");
	if(v != KV_FLOATERR) P->subshift[1] = v;
	v = getKeyValuedouble(kv, "dsub_z");
	if(v != KV_FLOATERR) P->subshift[2] = v;
	v = getKeyValuedouble(kv, "rsub_x");
	if(v != KV_FLOATERR) srx = v;
	v = getKeyValuedouble(kv, "rsub_y");
	if(v != KV_FLOATERR) sry = v;
	v = getKeyValuedouble(kv, "rsub_z");
	if(v != KV_FLOATERR) srz = v;
	v = getKeyValuedouble(kv, "rfeed_x");
	if(v != KV_FLOATERR) frx = v;
	v = getKeyValuedouble(kv, "rfeed_y");
	if(v != KV_FLOATERR) fry = v;
	v = getKeyValuedouble(kv, "rfeed_z");
	if(v != KV_FLOATERR) frz = v;
	v = getKeyValuedouble(kv, "focus");
	if(v != KV_FLOATERR) P->focus = v;
	v = getKeyValuedouble(kv, "az_offset");
	if(v != KV_FLOATERR) P->az_offset = v*M_PI/180.0;
	v = getKeyValuedouble(kv, "el_offset");
	if(v != KV_FLOATERR) P->el_offset = v*M_PI/180.0;
	v = getKeyValuedouble(kv, "phase_offset");
	if(v != KV_FLOATERR) P->phase_offset = v*M_PI/180.0;

	v = getKeyValuedouble(kv, "sub_h");
	if(v != KV_FLOATERR) 
	{
		P->subrotpoint[0] = 0.0;
		P->subrotpoint[1] = 0.0;
		P->subrotpoint[2] = v;	
	}

	V = getKeyValueVector(kv, "subrotpoint");
	if(V)
	{
		switch(VectorSize(V))
		{
		case 1:
			P->subrotpoint[0] = 0.0;
			P->subrotpoint[1] = 0.0;
			P->subrotpoint[2] = V[0];	
			break;
		case 2:
			fprintf(stderr, 
				"Warning : assuming a height for sub rot\n");
			P->subrotpoint[0] = V[0];
			P->subrotpoint[1] = V[1];
			break;
		case 3:
			P->subrotpoint[0] = V[0];
			P->subrotpoint[1] = V[1];
			P->subrotpoint[2] = V[2];	
			break;
		default:
			deletePathology(P);
			fprintf(stderr, "Pathology error : subrotpoint\n");
			return 0;
		}
	}

	deleteMatrix(P->subrot);
	P->subrot = newrotationMatrix(srx, sry, srz);
	deleteMatrix(P->feedrot);
	P->feedrot = newrotationMatrix(frx, fry, frz);

	return P;
}

void printPathology(const Pathology *P)
{
	printf("Pathology: %p\n", P);
	printf("  subrot = "); printMatrix(P->subrot);
	printf("  feedrot = "); printMatrix(P->feedrot);
	printf("  subshift = "); printVector(P->subshift);
	printf("  subrotpoint = "); printVector(P->subrotpoint);
	printf("  feedshift = "); printVector(P->feedshift);
	printf("\n");
}

void deletePathology(Pathology *P)
{
	g_assert(P);

	if(P->subrot) deleteMatrix(P->subrot);
	if(P->feedrot) deleteMatrix(P->feedrot);
	if(P->subshift) deleteVector(P->subshift);
	if(P->feedshift) deleteVector(P->feedshift);
	if(P->subrotpoint) deleteVector(P->subrotpoint);

	g_free(P);
}

static void normvec(const double *a, const double *b, double *c)
{
	int i;
	double r;
	for(i = 0; i < 3; i++) c[i] = a[i] - b[i];
	r = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
	for(i = 0; i < 3; i++) c[i] /= r;
}

double dAdOmega(const Antenna *a, const Ray *ray1, const Ray *ray2,
	const Ray *ray3, const Pathology *p)
{
	double A, Omega;
	double n1[3], n2[3], n3[3], f[3], ci, cj, ck;
	int i;

	/* Area in aperture is in a plane z = const */
	A = 0.5*fabs(
		(ray1->aper[0]-ray2->aper[0])*(ray1->aper[1]-ray3->aper[1]) -
		(ray1->aper[0]-ray3->aper[0])*(ray1->aper[1]-ray2->aper[1]) );

	for(i = 0; i < 3; i++) f[i] = a->feed[i] + p->feedshift[i];

	normvec(ray1->sub, f, n1);
	normvec(ray2->sub, f, n2);
	normvec(ray3->sub, f, n3);

	for(i = 0; i < 3; i++)
	{
		n1[i] -= n3[i];
		n2[i] -= n3[i];
	}
	
	ci = n1[1]*n2[2] - n1[2]*n2[1];
	cj = n1[2]*n2[0] - n1[0]*n2[2];
	ck = n1[0]*n2[1] - n1[1]*n2[0];
	
	Omega = 0.5*sqrt(ci*ci + cj*cj + ck*ck);
	
	return A/Omega;
}

double dOmega(const Antenna *a, const Ray *ray1, const Ray *ray2,
	const Ray *ray3, const Pathology *p)
{
	double Omega;
	double n1[3], n2[3], n3[3], f[3], ci, cj, ck;
	int i;

	for(i = 0; i < 3; i++) f[i] = a->feed[i] + p->feedshift[i];

	normvec(ray1->sub, f, n1);
	normvec(ray2->sub, f, n2);
	normvec(ray3->sub, f, n3);

	for(i = 0; i < 3; i++)
	{
		n1[i] -= n3[i];
		n2[i] -= n3[i];
	}
	
	ci = n1[1]*n2[2] - n1[2]*n2[1];
	cj = n1[2]*n2[0] - n1[0]*n2[2];
	ck = n1[0]*n2[1] - n1[1]*n2[0];
	
	Omega = 0.5*sqrt(ci*ci + cj*cj + ck*ck);
	
	return Omega;
}

double Raylen(const Ray *ray)
{
	double len = 0, d[3];
	int i;

	/* feed to subreflector */
	for(i = 0; i < 3; i++) 
		d[i] = ray->feed[i] - ray->sub[i];
	len += sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

	/* subreflector to dish */
	for(i = 0; i < 3; i++) 
		d[i] = ray->sub[i] - ray->dish[i];
	len += sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

	/* dish to aperture */
	for(i = 0; i < 3; i++) 
		d[i] = ray->dish[i] - ray->aper[i];
	len += sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
 
	return len;
}

void Pathologize(Vector sub, const Pathology *p)
{
	int i;
	int j;
	double tmp[6];
	
	for(i = 0; i < 3; i++) sub[i] -= p->subrotpoint[i];
	
	for(i = 0; i < 3; i++) 
	{
		tmp[i] = 0.0;
		tmp[i+3] = 0.0;
		for(j = 0; j < 3; j++) tmp[i] += p->subrot[i][j]*sub[j];
		for(j = 0; j < 3; j++) tmp[i+3] += p->subrot[i][j]*sub[j+3];
	}
	
	for(i = 0; i < 3; i++) 
		sub[i] = tmp[i] + p->subrotpoint[i] + p->subshift[i];
	for(i = 4; i < 6; i++) 
		sub[i] = tmp[i];
}

void defocusAntenna(Pathology *P, const Antenna *a)
{
	double dx[3];
	int i;
	
	dx[0] = -a->feed[0];
	dx[1] = -a->feed[1];
	dx[2] = a->sub_h-a->feed[2];
	norm3(dx);
	for(i = 0; i < 3; i++) P->feedshift[i] += P->focus*dx[i];

	P->focus = 0.0;
}


void intersectdish(const Antenna *a, const Vector sub, const Vector unitdir,
	Vector dish, int niter)
{
	double A, B, C, t, m, r;
	double x[3], n[3];
	int i, iter;
	
	/* First intersect with ideal paraboloid */
	A = a->bestparabola*(unitdir[0]*unitdir[0] + unitdir[1]*unitdir[1]);
	B = 2.0*a->bestparabola*(unitdir[0]*sub[0] + unitdir[1]*sub[1])
		-unitdir[2];
	C = a->bestparabola*(sub[0]*sub[0] + sub[1]*sub[1]) - sub[2];
	t = 0.5*(sqrt(B*B-4.0*A*C) - B)/A; /* take greater root */
	
	for(iter = 0; ; iter++)
	{
		/* get position (x) and normal (n) on the real dish */
		for(i = 0; i < 2; i++) x[i] = sub[i] + t*unitdir[i];
		r = sqrt(x[0]*x[0] + x[1]*x[1]);
		dishvalue(a, r, &(x[2]), &m);
		n[2] = 1.0/sqrt(1.0+m*m);
		n[0] = -m*(x[0]/r)*n[2];
		n[1] = -m*(x[1]/r)*n[2];

		if(iter >= niter) break;

		A = B = 0;
		for(i = 0; i < 3; i++)
		{
			A += n[i]*(x[i]-sub[i]);	/* n dot (x-sub) */
			B += n[i]*unitdir[i];		/* n dot unitdir */
		}
		t = A/B;
	}

	for(i = 0; i < 3; i++)
	{
		dish[i] = x[i];
		dish[i+3] = n[i];
	}
}

void intersectaperture(const Antenna *a, const Vector dish, 
	const Vector unitdir, Vector aper)
{
	double t;
	int i;
	
	t = (a->zedge-dish[2])/unitdir[2];
	for(i = 0; i < 3; i++) aper[i] = dish[i] + t*unitdir[i];
	
	aper[3] = aper[4] = 0.0;
	aper[5] = 1.0;
}

/* gain in power */
double feedfunc(const Antenna *a, double theta)
{
	double stheta;

	if(a->feedpattern == 0)
	{
		stheta = sin(theta);
		return exp(2.0*(-0.083)*a->fa2pi*a->fa2pi*stheta*stheta);
	}

	else return pow(10.0, 0.1*interpolateVector(a->feedpattern,
		theta/a->feedpatterndelta));
}

/* gain in power */
double feedgain(const Antenna *a, const Ray *ray, const Pathology *p)
{
	double theta, costheta = 0.0;
	double v[3], f[3]; 
	int i, j;

	for(i = 0; i < 3; i++) 
		v[i] = ray->sub[i] - ray->feed[i];
	norm3(v);

	for(i = 0; i < 3; i++) f[i] = 0.0;
	for(j = 0; j < 3; j++) for(i = 0; i < 3; i++)
		f[j] += p->feedrot[j][i]*a->feeddir[i];

	for(i = 0; i < 3; i++) costheta += f[i]*v[i];

	if(a->feedpattern == 0)
	    return exp(2.0*(-0.083)*a->fa2pi*a->fa2pi*(1.0-costheta*costheta));

	theta = acos(costheta);
	return pow(10.0, 0.1*interpolateVector(a->feedpattern, 
		theta/a->feedpatterndelta));
}

/* returns integral of feedgain(theta)*sin(theta) over theta = 0 to pi/2 
 * this is equivalent to the power over the forward hemisphere
 */
double feedpower(const Antenna *a)
{
	const double dtheta = 0.00005;
	double theta, stheta, thetamax;
	double v, sum = 0.0;

	if(a->feedpattern == 0)
	{
		for(theta = 0.5*dtheta; theta < 0.5*M_PI; theta+=dtheta)
		{
			stheta = sin(theta);
			v = exp(2.0*(-0.083)*a->fa2pi*a->fa2pi*stheta*stheta);
			sum += v*stheta*dtheta;
		}
	}
	else
	{
		thetamax = a->feedpatterndelta*(VectorSize(a->feedpattern)-1.0);
		if(thetamax > M_PI/2.0) thetamax = M_PI/2.0;
		for(theta = 0.5*dtheta; theta < thetamax; theta+=dtheta)
		{
			stheta = sin(theta);
			v = pow(10.0, 0.1*interpolateVector(a->feedpattern,
				theta/a->feedpatterndelta));
			sum += v*stheta*dtheta;
		}
	}
	return 2.0*M_PI*sum;
}

double feedbackpower(const Antenna *a)
{
	const double dtheta = 0.00005;
	double theta, stheta, thetamax;
	double v, sum = 0.0;
	
	if(a->feedpattern == 0) return 0.0;
	else
	{
		thetamax = a->feedpatterndelta*(VectorSize(a->feedpattern)-1.0);
		if(thetamax < M_PI/2.0) return 0.0;
		if(thetamax > M_PI) thetamax = M_PI;
		
		for(theta = M_PI/2.0; theta < thetamax; theta+=dtheta)
		{
			stheta = sin(theta);
			v = pow(10.0, 0.1*interpolateVector(a->feedpattern,
				theta/a->feedpatterndelta));
			sum += v*stheta*dtheta;
		}
		return 2.0*M_PI*sum;
	}
}

Ray *trace(const Antenna *a, double x, double y, const Pathology *p)
{
	Ray *ray;
	Vector idealsub;
	double fu[3], du[3], au[3], ndotf=0.0, ndotd=0.0;
	int i;
	const int niter = 7;

	idealsub = newVector(6);

	subfromdish(a, x, y, idealsub);

	ray = newRay(idealsub);
	deleteVector(idealsub);

	Pathologize(ray->sub, p);

	if(ray->sub[5] < -1.0 || ray->sub[5] > -0.0) 
	{
		deleteRay(ray);
		return 0;
	}

	for(i = 0; i < 3; i++) 
		ray->feed[i] = a->feed[i] + p->feedshift[i];
	
	/* now determine unit vector pointing to dish */

	/* unit toward feed */
	for(i = 0; i < 3; i++) 
		fu[i] = ray->feed[i] - ray->sub[i];

	norm3(fu);

	/* unit toward dish */
	for(i = 0; i < 3; i++) ndotf += ray->sub[i+3]*fu[i];
	for(i = 0; i < 3; i++) du[i] = 2.0*ray->sub[i+3]*ndotf - fu[i];

	/* dish point */
	intersectdish(a, ray->sub, du, ray->dish, niter);

	/* unit toward aperture */
	for(i = 0; i < 3; i++) ndotd += ray->dish[i+3]*du[i];
	for(i = 0; i < 3; i++) au[i] = du[i] - 2.0*ray->dish[i+3]*ndotd;

	intersectaperture(a, ray->dish, au, ray->aper);

	return ray;
}

#if 0
Ray *supertrace(const Antenna *a, const Pathology *p, int niter,
	double x, double y, double *x0, double *y0)
{
	Ray *ray;
	int iter;
	double x1, y1, R2;

	R2 = a->radius*a->radius;

	x1 = x;
	y1 = y;

	for(iter = 0; iter < niter; iter++)
	{
		ray = trace(a, x1, y1, p);
		if(!ray) return 0;
		
		x1 += (x - ray->aper[0]);
		y1 += (y - ray->aper[1]);
		deleteRay(ray);

		if(x1*x1 + y1*y1 > R2) return 0;
	}

	if(x0) *x0 = x1;
	if(y0) *y0 = y1;

	return trace(a, x1, y1, p);
}
#endif
Vector tracepol(const Vector E0, const Ray *ray)
{
	Vector E;
	int i;
	double v1[3], v2[3], v3[3], r[3], C[2];

	for(i = 0; i < 3; i++)
	{
		v1[i] = ray->sub[i]  - ray->feed[i];
		v2[i] = ray->dish[i] - ray->sub[i];
		v3[i] = ray->aper[i] - ray->dish[i];
	}
	norm3(v1);
	norm3(v2);
	norm3(v3);

	E = dupVector(E0);

	for(i = 0; i < 3; i++) r[i] = v1[i] - v2[i];
	norm3(r); 
	C[0] = r[0]*E[0] + r[1]*E[2] + r[2]*E[4];
	C[1] = r[0]*E[1] + r[1]*E[3] + r[2]*E[5];
	for(i = 0; i < 6; i++) E[i] -= 2.0*r[i/2]*C[i%2];
	for(i = 0; i < 6; i++) E[i] = -E[i];

	for(i = 0; i < 3; i++) r[i] = v2[i] - v3[i];
	norm3(r); 
	C[0] = r[0]*E[0] + r[1]*E[2] + r[2]*E[4];
	C[1] = r[0]*E[1] + r[1]*E[3] + r[2]*E[5];
	for(i = 0; i < 6; i++) E[i] -= 2.0*r[i/2]*C[i%2];
	for(i = 0; i < 6; i++) E[i] = -E[i];

	return E;
}

int legplanewaveblock(const Antenna *a, double x, double y)
{
	/* outside the leg foot area, the blockage is spherical wave */
	if(x*x + y*y > a->legfoot*a->legfoot) return 0;
	
	if(a->legwidth == 0.0) return 0;

	if(strcmp(a->name, "VLBA") == 0) 
	{
		const double s=1.457937;
		const double c=1.369094;
		if(fabs(s*x+c*y) < -a->legwidth) return 1;
		if(fabs(s*x-c*y) < -a->legwidth) return 1;
	}
	else if(a->legwidth < 0.0)  /* "x shaped" legs */
	{
		if(fabs(x-y)*M_SQRT2 < -a->legwidth) return 1;
		if(fabs(x+y)*M_SQRT2 < -a->legwidth) return 1;
	}
	else if(a->legwidth > 0.0) /* "+ shaped" legs */
	{
		if(fabs(x)*2.0 < a->legwidth) return 1;
		if(fabs(y)*2.0 < a->legwidth) return 1;
	}

	return 0;
}

int legplanewaveblock2(const Antenna *a, const Ray *ray)
{
	int i, n;
	double dr2;
	double theta, phi;
	double r0[3], dr[3], l0[3], l1[3], dl[3], D[3]; 
	double D2, N[3], ll, rr;
	const double thetaplus[4] = 
		{0, M_PI/2.0, M_PI, 3.0*M_PI/2.0};
	const double thetacross[4] = 
		{0.25*M_PI, 0.75*M_PI, 1.25*M_PI, 1.75*M_PI};
	const double thetavlba[4] =
		{0.816817, 2.3247756, 3.9584096, 5.466368};
	const double *thetalut;

	if(a->legwidth == 0.0) return 0;

	if(strcmp(a->name, "VLBA") == 0) thetalut = thetavlba;
	else if(a->legwidth < 0.0) thetalut = thetacross;
	else thetalut = thetaplus;

	/* inside the leg feet is plane wave blockage */
	dr2 = ray->dish[0]*ray->dish[0] + ray->dish[1]*ray->dish[1];
	if(dr2 >= a->legfoot*a->legfoot) return 0;

	for(i = 0; i < 3; i++)
	{
		r0[i] = ray->dish[i];
		dr[i] = ray->aper[i] - r0[i];
	}
	rr = r0[0]*r0[0] + r0[1]*r0[1];

	l0[2] = a->legfootz;
	l1[0] = l1[1] = 0.0;
	l1[2] = a->legapex;
	phi = atan2(r0[1], r0[0]);

	for(n = 0; n < 4; n++)
	{
		theta = thetalut[n];
		l0[0] = a->legfoot*cos(theta);
		l0[1] = a->legfoot*sin(theta);
		ll = l0[0]*l0[0] + l0[1]*l0[1];
		if((l0[0]*r0[0] + l0[1]*r0[1])/sqrt(ll*rr) < 0.7) continue;
		for(i = 0; i < 3; i++) dl[i] = l1[i] - l0[i];
		for(i = 0; i < 3; i++) D[i] = r0[i] - l0[i];

		N[0] = dr[1]*dl[2] - dr[2]*dl[1];
		N[1] = dr[2]*dl[0] - dr[0]*dl[2];
		N[2] = dr[0]*dl[1] - dr[1]*dl[0];
		norm3(N);

		D2 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];

		if(fabs(D2) <= 0.5*fabs(a->legwidth)) return 1;
	}
	
	return 0;
}

int legsphericalwaveblock(const Antenna *a, const Ray *ray)
{
	int i, n;
	double dr2;
	double theta, phi;
	double r0[3], dr[3], l0[3], l1[3], dl[3], D[3]; 
	double D2, N[3], ll, rr;
	const double thetaplus[4] = 
		{0, M_PI/2.0, M_PI, 3.0*M_PI/2.0};
	const double thetacross[4] = 
		{0.25*M_PI, 0.75*M_PI, 1.25*M_PI, 1.75*M_PI};
	const double thetavlba[4] =
		{0.816817, 2.3247756, 3.9584096, 5.466368};
	const double *thetalut;

	if(a->legwidth == 0.0) return 0;

	if(strcmp(a->name, "VLBA") == 0) thetalut = thetavlba;
	else if(a->legwidth < 0.0) thetalut = thetacross;
	else thetalut = thetaplus;

	/* inside the leg feet is plane wave blockage */
	dr2 = ray->dish[0]*ray->dish[0] + ray->dish[1]*ray->dish[1];
	if(dr2 < a->legfoot*a->legfoot) return 0;

	for(i = 0; i < 3; i++)
	{
		r0[i] = ray->dish[i];
		dr[i] = ray->sub[i] - r0[i];
	}
	rr = r0[0]*r0[0] + r0[1]*r0[1];

	l0[2] = a->legfootz;
	l1[0] = l1[1] = 0.0;
	l1[2] = a->legapex;
	phi = atan2(r0[1], r0[0]);

	for(n = 0; n < 4; n++)
	{
		theta = thetalut[n];
		l0[0] = a->legfoot*cos(theta);
		l0[1] = a->legfoot*sin(theta);
		ll = l0[0]*l0[0] + l0[1]*l0[1];
		if((l0[0]*r0[0] + l0[1]*r0[1])/sqrt(ll*rr) < 0.7) continue;
		for(i = 0; i < 3; i++) dl[i] = l1[i] - l0[i];
		for(i = 0; i < 3; i++) D[i] = r0[i] - l0[i];

		N[0] = dr[1]*dl[2] - dr[2]*dl[1];
		N[1] = dr[2]*dl[0] - dr[0]*dl[2];
		N[2] = dr[0]*dl[1] - dr[1]*dl[0];
		norm3(N);

		D2 = D[0]*N[0] + D[1]*N[1] + D[2]*N[2];

		if(fabs(D2) <= 0.5*fabs(a->legwidth)) return 1;
	}
	
	return 0;
}
