/*
 * lidar_map.c — reference cloud construction and the 3D Euclidean distance
 * field used as the likelihood-field measurement model.
 */

#include "lidar_map.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Reference cloud: scan the world from several viewpoints             */
/* ------------------------------------------------------------------ */

void lidar_map_build_cloud(World *w, int grid, ScanCloud *out)
{
	/*
	 * Scan from a 3D grid of interior viewpoints (grid x grid in the floor
	 * plane, at 3 heights). Near-wall viewpoints sample their nearby wall at
	 * sub-voxel spacing, and the union of all views covers every surface
	 * densely — so the resulting distance field approximates true distance to
	 * surface rather than distance to the nearest sparse sample. Use a dense
	 * LiDAR pattern (world_set_lidar) before calling this.
	 */
	double mn[3], mx[3];
	world_bounds(w, mn, mx);
	if (grid < 2) grid = 2;

	/* interior span: keep ~12% margin off the walls */
	double lox = mn[0] + 0.12*(mx[0]-mn[0]), hix = mx[0] - 0.12*(mx[0]-mn[0]);
	double loy = mn[1] + 0.12*(mx[1]-mn[1]), hiy = mx[1] - 0.12*(mx[1]-mn[1]);
	double heights[3] = {
		mn[2] + 0.25*(mx[2]-mn[2]),
		mn[2] + 0.50*(mx[2]-mn[2]),
		mn[2] + 0.75*(mx[2]-mn[2]),
	};

	ScanCloud tmp;
	scancloud_init(&tmp);

	for (int h = 0; h < 3; h++)
	for (int iy = 0; iy < grid; iy++)
	for (int ix = 0; ix < grid; ix++) {
		Pose p = pose_identity();   /* full-sphere LiDAR: orientation irrelevant */
		p.pos[0] = lox + (hix - lox) * ix / (grid - 1);
		p.pos[1] = loy + (hiy - loy) * iy / (grid - 1);
		p.pos[2] = heights[h];

		world_scan(w, &p, &tmp);
		scancloud_reserve(out, out->n + tmp.n);
		for (int k = 0; k < tmp.n; k++) {
			double local[3] = { tmp.x[k], tmp.y[k], tmp.z[k] };
			double world[3];
			pose_transform_point(&p, local, world);
			scancloud_push(out, (float)world[0], (float)world[1], (float)world[2]);
		}
	}

	scancloud_free(&tmp);
}

/* ------------------------------------------------------------------ */
/* Distance field                                                     */
/* ------------------------------------------------------------------ */

#define DF_INF 1.0e10

/*
 * 1D squared-distance transform of a sampled function (Felzenszwalb &
 * Huttenlocher, 2012). f/d are length n; v (length n) and z (length n+1) are
 * caller-provided scratch. d[q] = min_i ( f[i] + (q-i)^2 ).
 */
static void edt_1d(const double *f, int n, double *d, int *v, double *z)
{
	int k = 0;
	v[0] = 0;
	z[0] = -DF_INF;
	z[1] =  DF_INF;
	for (int q = 1; q < n; q++) {
		double s = ((f[q] + (double)q*q) - (f[v[k]] + (double)v[k]*v[k]))
		         / (2.0*q - 2.0*v[k]);
		while (s <= z[k]) {
			k--;
			s = ((f[q] + (double)q*q) - (f[v[k]] + (double)v[k]*v[k]))
			  / (2.0*q - 2.0*v[k]);
		}
		k++;
		v[k] = q;
		z[k] = s;
		z[k+1] = DF_INF;
	}
	k = 0;
	for (int q = 0; q < n; q++) {
		while (z[k+1] < q) k++;
		double dq = q - v[k];
		d[q] = dq*dq + f[v[k]];
	}
}

static inline long df_idx(const DistanceField *df, int x, int y, int z)
{
	return ((long)z * df->dim[1] + y) * df->dim[0] + x;
}

int df_build(DistanceField *df, const ScanCloud *cloud, double voxel, double pad)
{
	if (cloud->n <= 0 || voxel <= 0.0) return -1;

	/* bounding box of the reference cloud */
	double mn[3] = { cloud->x[0], cloud->y[0], cloud->z[0] };
	double mx[3] = { cloud->x[0], cloud->y[0], cloud->z[0] };
	for (int i = 1; i < cloud->n; i++) {
		double p[3] = { cloud->x[i], cloud->y[i], cloud->z[i] };
		for (int k = 0; k < 3; k++) {
			if (p[k] < mn[k]) mn[k] = p[k];
			if (p[k] > mx[k]) mx[k] = p[k];
		}
	}

	df->voxel = voxel;
	for (int k = 0; k < 3; k++) {
		double lo = mn[k] - pad;
		df->origin[k] = lo;
		df->dim[k] = (int)ceil((mx[k] + pad - lo) / voxel) + 1;
		if (df->dim[k] < 1) df->dim[k] = 1;
	}

	long N = (long)df->dim[0] * df->dim[1] * df->dim[2];
	df->dist = malloc((size_t)N * sizeof(float));
	if (!df->dist) return -1;

	/* occupancy: 0 at surface cells, +INF elsewhere (squared-distance domain) */
	for (long m = 0; m < N; m++) df->dist[m] = (float)DF_INF;
	for (int i = 0; i < cloud->n; i++) {
		int ix = (int)((cloud->x[i] - df->origin[0]) / voxel + 0.5);
		int iy = (int)((cloud->y[i] - df->origin[1]) / voxel + 0.5);
		int iz = (int)((cloud->z[i] - df->origin[2]) / voxel + 0.5);
		if (ix < 0 || ix >= df->dim[0] || iy < 0 || iy >= df->dim[1] ||
		    iz < 0 || iz >= df->dim[2]) continue;
		df->dist[df_idx(df, ix, iy, iz)] = 0.0f;
	}

	/* separable 3D EDT: transform along x, then y, then z */
	int maxd = df->dim[0];
	if (df->dim[1] > maxd) maxd = df->dim[1];
	if (df->dim[2] > maxd) maxd = df->dim[2];
	double *f = malloc((size_t)maxd * sizeof(double));
	double *d = malloc((size_t)maxd * sizeof(double));
	int    *v = malloc((size_t)maxd * sizeof(int));
	double *z = malloc((size_t)(maxd + 1) * sizeof(double));
	if (!f || !d || !v || !z) { free(f); free(d); free(v); free(z); return -1; }

	const int DX = df->dim[0], DY = df->dim[1], DZ = df->dim[2];

	/* along x */
	for (int zz = 0; zz < DZ; zz++)
		for (int yy = 0; yy < DY; yy++) {
			for (int xx = 0; xx < DX; xx++) f[xx] = df->dist[df_idx(df, xx, yy, zz)];
			edt_1d(f, DX, d, v, z);
			for (int xx = 0; xx < DX; xx++) df->dist[df_idx(df, xx, yy, zz)] = (float)d[xx];
		}
	/* along y */
	for (int zz = 0; zz < DZ; zz++)
		for (int xx = 0; xx < DX; xx++) {
			for (int yy = 0; yy < DY; yy++) f[yy] = df->dist[df_idx(df, xx, yy, zz)];
			edt_1d(f, DY, d, v, z);
			for (int yy = 0; yy < DY; yy++) df->dist[df_idx(df, xx, yy, zz)] = (float)d[yy];
		}
	/* along z */
	for (int yy = 0; yy < DY; yy++)
		for (int xx = 0; xx < DX; xx++) {
			for (int zz = 0; zz < DZ; zz++) f[zz] = df->dist[df_idx(df, xx, yy, zz)];
			edt_1d(f, DZ, d, v, z);
			for (int zz = 0; zz < DZ; zz++) df->dist[df_idx(df, xx, yy, zz)] = (float)d[zz];
		}

	free(f); free(d); free(v); free(z);

	/* squared cell distance -> metric distance */
	for (long m = 0; m < N; m++) df->dist[m] = (float)(sqrt((double)df->dist[m]) * voxel);

	return 0;
}

void df_free(DistanceField *df)
{
	free(df->dist);
	df->dist = NULL;
	df->dim[0] = df->dim[1] = df->dim[2] = 0;
}
