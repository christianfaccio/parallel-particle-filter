#ifndef LIDAR_MAP_H_
#define LIDAR_MAP_H_

#include "world.h"
#include "pose.h"

/*
 * The localization map: a dense 3D voxel grid storing, at each cell, the
 * Euclidean distance to the nearest surface point of a reference cloud.
 *
 * This is the "likelihood field" measurement model (Thrun, Probabilistic
 * Robotics). The reference cloud is built once by scanning the ground-truth
 * world from several viewpoints and transforming the hits to the world frame.
 * The distance field is then produced with a separable exact Euclidean
 * distance transform (Felzenszwalb–Huttenlocher), which is far cheaper than a
 * brute-force cells x points scan and is itself a good future parallel kernel.
 *
 * Scoring a transformed scan point is then a single O(1) grid lookup.
 */

typedef struct {
	double origin[3];   /* world coords of voxel (0,0,0) center */
	double voxel;       /* edge length [m]                      */
	int    dim[3];      /* grid size in x, y, z                 */
	float *dist;        /* dim[0]*dim[1]*dim[2] distances [m]   */
} DistanceField;

/*
 * Build the reference cloud by scanning `w` from a grid x grid x 3 lattice of
 * interior viewpoints, accumulating world-frame hit points into `out`. Set a
 * dense LiDAR pattern with world_set_lidar() first so surfaces are sampled
 * below the voxel size. `out` must be initialized (scancloud_init) by the
 * caller, which also owns/frees it.
 */
void lidar_map_build_cloud(World *w, int grid, ScanCloud *out);

/*
 * Build a distance field over the cloud's bounding box (padded by `pad` m) at
 * the given voxel size. Returns 0 on success. Free with df_free.
 */
int  df_build(DistanceField *df, const ScanCloud *cloud, double voxel, double pad);
void df_free(DistanceField *df);

/* Distance [m] to the nearest surface at world point (x,y,z). Points outside
 * the grid are clamped to the boundary. static inline so the N_particles x
 * M_points scoring loop in pf_update_weights inlines it (no call overhead). */
static inline float df_lookup(const DistanceField *df, double x, double y, double z)
{
	int ix = (int)((x - df->origin[0]) / df->voxel + 0.5);
	int iy = (int)((y - df->origin[1]) / df->voxel + 0.5);
	int iz = (int)((z - df->origin[2]) / df->voxel + 0.5);
	if (ix < 0) ix = 0; else if (ix >= df->dim[0]) ix = df->dim[0] - 1;
	if (iy < 0) iy = 0; else if (iy >= df->dim[1]) iy = df->dim[1] - 1;
	if (iz < 0) iz = 0; else if (iz >= df->dim[2]) iz = df->dim[2] - 1;
	return df->dist[((long)iz * df->dim[1] + iy) * df->dim[0] + ix];
}

#endif /* LIDAR_MAP_H_ */
