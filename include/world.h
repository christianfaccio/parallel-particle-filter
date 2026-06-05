#ifndef WORLD_H_
#define WORLD_H_

#include "pose.h"

/*
 * Thin wrapper over the vendored lidarsim simulator (third_party/lidarsim).
 *
 * This header deliberately does NOT include lidarsim.h: lidarsim is a
 * header-with-definitions (not STB-style), so it is compiled in exactly one
 * translation unit (src/world.c). The rest of the project sees only the opaque
 * `World` handle and the plain SoA `ScanCloud` below.
 *
 * The `World` is the ground-truth environment. A scan taken at a sensor pose
 * returns the hit points in the SENSOR frame (vehicle-centric) — exactly what
 * a real LiDAR delivers — which the particle filter then re-projects through
 * each particle's hypothesized pose.
 */

/* Point cloud in struct-of-arrays layout (SIMD/GPU friendly for later). */
typedef struct {
	float *x;
	float *y;
	float *z;
	int    n;
	int    cap;
} ScanCloud;

void scancloud_init(ScanCloud *c);
void scancloud_reserve(ScanCloud *c, int n);
void scancloud_push(ScanCloud *c, float x, float y, float z);
void scancloud_free(ScanCloud *c);

/* Write the cloud to an ASCII .pcd file (PCL / CloudCompare / Open3D readable),
 * keeping every `stride`-th point (stride <= 1 writes all). */
void scancloud_save_pcd(const ScanCloud *c, const char *path, int stride);

typedef struct World World;   /* opaque */

/*
 * Build the ground-truth scene (primitives only — no meshes, to avoid
 * lidarsim's per-update mesh-BVH reallocation) and a programmatic spherical
 * LiDAR pattern with about n_az * n_el rays.
 */
World *world_create(int n_az, int n_el);
void   world_destroy(World *w);

/* Swap the LiDAR ray pattern (e.g. a dense pattern for one-time map building,
 * then back to the runtime pattern). */
void   world_set_lidar(World *w, int n_az, int n_el);

/*
 * Simulate a LiDAR scan from `sensor`. Points are returned in the sensor frame
 * (origin at the sensor, axes aligned with the sensor). `out` is cleared and
 * refilled; reuse it across calls to amortize allocation.
 */
void world_scan(World *w, const Pose *sensor, ScanCloud *out);

/* Axis-aligned bounding box of the scene geometry, in world coordinates. */
void world_bounds(const World *w, double min[3], double max[3]);

#endif /* WORLD_H_ */
