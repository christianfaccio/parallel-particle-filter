/*
 * world.c — the ONLY translation unit that includes lidarsim.
 *
 * lidarsim.h contains function *definitions* (not static inline, not STB-style),
 * so including it in more than one TU would cause multiple-definition link
 * errors. Everything lidarsim-specific is therefore confined here and exposed
 * through the opaque World / plain ScanCloud in world.h.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * Portability shim: lidarsim calls qsort_r with the glibc signature
 *   qsort_r(base, n, size, cmp(a,b,arg), arg)
 * which is correct in the project's Linux/glibc container but does not exist
 * on macOS (whose qsort_r has a different argument order). On Apple we redirect
 * the call to a plain qsort via a file-local trampoline. Single-threaded scene
 * build, so the static state is fine. On Linux this block compiles out.
 */
#ifdef __APPLE__
static int (*g_ls_cmp)(const void *, const void *, void *);
static void *g_ls_arg;
static int ls_cmp2(const void *a, const void *b) { return g_ls_cmp(a, b, g_ls_arg); }
static void ls_qsort_r(void *base, size_t n, size_t sz,
                       int (*cmp)(const void *, const void *, void *), void *arg)
{
	g_ls_cmp = cmp;
	g_ls_arg = arg;
	qsort(base, n, sz, ls_cmp2);
}
#define qsort_r ls_qsort_r
#endif

/* Vendored third-party header: silence its internal warnings, keep ours on. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "lidarsim.h"
#pragma GCC diagnostic pop

#include "world.h"

#include <math.h>

/* ------------------------------------------------------------------ */
/* ScanCloud (struct-of-arrays)                                       */
/* ------------------------------------------------------------------ */

void scancloud_init(ScanCloud *c)
{
	c->x = c->y = c->z = NULL;
	c->n = 0;
	c->cap = 0;
}

void scancloud_reserve(ScanCloud *c, int n)
{
	if (n <= c->cap) return;
	int cap = c->cap ? c->cap : 1024;
	while (cap < n) cap *= 2;
	c->x = realloc(c->x, (size_t)cap * sizeof(float));
	c->y = realloc(c->y, (size_t)cap * sizeof(float));
	c->z = realloc(c->z, (size_t)cap * sizeof(float));
	c->cap = cap;
}

void scancloud_push(ScanCloud *c, float x, float y, float z)
{
	scancloud_reserve(c, c->n + 1);
	c->x[c->n] = x;
	c->y[c->n] = y;
	c->z[c->n] = z;
	c->n++;
}

void scancloud_free(ScanCloud *c)
{
	free(c->x); free(c->y); free(c->z);
	scancloud_init(c);
}

void scancloud_save_pcd(const ScanCloud *c, const char *path, int stride)
{
	if (stride < 1) stride = 1;
	FILE *f = fopen(path, "w");
	if (!f) { perror(path); return; }

	int count = 0;
	for (int i = 0; i < c->n; i += stride) count++;

	fprintf(f, "# .PCD v.7 - Point Cloud Data file format\n");
	fprintf(f, "VERSION .7\nFIELDS x y z\nSIZE 4 4 4\nTYPE F F F\nCOUNT 1 1 1\n");
	fprintf(f, "WIDTH %d\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\n", count);
	fprintf(f, "POINTS %d\nDATA ascii\n", count);
	for (int i = 0; i < c->n; i += stride)
		fprintf(f, "%f %f %f\n", c->x[i], c->y[i], c->z[i]);

	fclose(f);
}

/* ------------------------------------------------------------------ */
/* World                                                              */
/* ------------------------------------------------------------------ */

struct World {
	scene   scn;
	lidar   lid;
	mat4x4 *world_tf;   /* per-object world transform, restored before each scan */
	int     n_obj;
	double  bmin[3], bmax[3];
};

/* Pose -> lidarsim rigid transform (rotation in top-left 3x3, translation in
 * column 3). Matches mat4x4_mul_vec3's row-major convention. */
static mat4x4 pose_to_mat4x4(const Pose *p)
{
	double R[9];
	quat_rotation_matrix(p->quat, R);
	mat4x4 m;
	for (int i = 0; i < 3; i++) {
		m.m[i][0] = (float)R[i*3 + 0];
		m.m[i][1] = (float)R[i*3 + 1];
		m.m[i][2] = (float)R[i*3 + 2];
		m.m[i][3] = (float)p->pos[i];
	}
	m.m[3][0] = 0.0f; m.m[3][1] = 0.0f; m.m[3][2] = 0.0f; m.m[3][3] = 1.0f;
	return m;
}

static void bounds_expand(World *w, double cx, double cy, double cz,
                          double hx, double hy, double hz)
{
	double lo[3] = {cx - hx, cy - hy, cz - hz};
	double hi[3] = {cx + hx, cy + hy, cz + hz};
	for (int k = 0; k < 3; k++) {
		if (lo[k] < w->bmin[k]) w->bmin[k] = lo[k];
		if (hi[k] > w->bmax[k]) w->bmax[k] = hi[k];
	}
}

static void add_sphere(World *w, double radius, double cx, double cy, double cz)
{
	scene_object obj = create_sphere((float)radius);
	obj.transform.m[0][3] = (float)cx;
	obj.transform.m[1][3] = (float)cy;
	obj.transform.m[2][3] = (float)cz;
	scene_add_object(&w->scn, obj);
	bounds_expand(w, cx, cy, cz, radius, radius, radius);
}

static void set_tri(triangle *t, vec3 a, vec3 b, vec3 c)
{
	t->vertices[0] = a; t->vertices[1] = b; t->vertices[2] = c;
}

/*
 * Room shell as a 12-triangle MESH rather than Box walls. lidarsim's Box (and
 * Plane) intersection orients a *transformed* primitive by the transpose of the
 * transform's rotation — fine when axis-aligned, wrong once the sensor (hence
 * the scene, via inv(L)) is rotated. Meshes transform their vertices and use
 * ray_triangle, which is correct under any rotation; Spheres are rotation-
 * invariant. So the world is built only from those two primitives.
 */
static void add_room_mesh(World *w, double hx, double hy, double zlo, double zhi)
{
	scene_object obj;
	obj.type = MESH;
	obj.mesh.triangle_count = 12;
	obj.mesh.original_triangles    = malloc(12 * sizeof(triangle));
	obj.mesh.transformed_triangles = malloc(12 * sizeof(triangle));
	obj.mesh.root = NULL;
	mat4x4_identity(&obj.transform);   /* room is already in world coords */

	vec3 v000 = cast_vec3(-hx,-hy,zlo), v100 = cast_vec3(hx,-hy,zlo);
	vec3 v110 = cast_vec3( hx, hy,zlo), v010 = cast_vec3(-hx,hy,zlo);
	vec3 v001 = cast_vec3(-hx,-hy,zhi), v101 = cast_vec3(hx,-hy,zhi);
	vec3 v111 = cast_vec3( hx, hy,zhi), v011 = cast_vec3(-hx,hy,zhi);

	triangle *T = obj.mesh.original_triangles;
	set_tri(&T[0],  v000,v100,v110); set_tri(&T[1],  v000,v110,v010); /* floor   */
	set_tri(&T[2],  v001,v101,v111); set_tri(&T[3],  v001,v111,v011); /* ceiling */
	set_tri(&T[4],  v000,v010,v011); set_tri(&T[5],  v000,v011,v001); /* x = -hx */
	set_tri(&T[6],  v100,v110,v111); set_tri(&T[7],  v100,v111,v101); /* x = +hx */
	set_tri(&T[8],  v000,v100,v101); set_tri(&T[9],  v000,v101,v001); /* y = -hy */
	set_tri(&T[10], v010,v110,v111); set_tri(&T[11], v010,v111,v011); /* y = +hy */

	scene_add_object(&w->scn, obj);
	bounds_expand(w, 0, 0, 0.5*(zlo+zhi), hx, hy, 0.5*(zhi-zlo));
}

/* A 20 x 20 x 3 m room (mesh shell) plus asymmetric spheres so the 6-DOF pose
 * is fully observable (a symmetric empty box would be ambiguous). */
static void build_scene(World *w)
{
	add_room_mesh(w, 10.0, 10.0, 0.0, 3.0);

	add_sphere(w, 0.8,  4.0,  3.0, 1.0);
	add_sphere(w, 0.6, -5.0, -2.0, 0.9);
	add_sphere(w, 1.0,  2.0, -6.0, 1.5);
	add_sphere(w, 0.5, -6.0,  5.0, 0.7);
	add_sphere(w, 0.7,  6.5, -4.0, 2.0);
}

/* Programmatic spherical LiDAR pattern: n_az azimuths x n_el elevation rings. */
static void build_lidar(World *w, int n_az, int n_el)
{
	const double el_min = -75.0 * M_PI / 180.0;
	const double el_max =  75.0 * M_PI / 180.0;
	uint32_t rc = (uint32_t)n_az * (uint32_t)n_el;

	lidar_free(&w->lid);   /* drop any previous pattern */
	w->lid.rays = malloc((size_t)rc * sizeof(ray));
	uint32_t idx = 0;
	for (int j = 0; j < n_el; j++) {
		double el = el_min + (el_max - el_min) * (j + 0.5) / n_el;
		for (int i = 0; i < n_az; i++) {
			double az = 2.0 * M_PI * (i + 0.5) / n_az;
			double dx = cos(el) * cos(az);
			double dy = cos(el) * sin(az);
			double dz = sin(el);
			/* keep components away from exactly zero so inv_dir stays finite
			 * and ray_aabb never forms 0*inf = NaN. */
			if (fabs(dx) < 1e-4) dx = (dx < 0 ? -1e-4 : 1e-4);
			if (fabs(dy) < 1e-4) dy = (dy < 0 ? -1e-4 : 1e-4);
			if (fabs(dz) < 1e-4) dz = (dz < 0 ? -1e-4 : 1e-4);
			w->lid.rays[idx].ori     = cast_vec3(0.0f, 0.0f, 0.0f);
			w->lid.rays[idx].dir     = cast_vec3((float)dx, (float)dy, (float)dz);
			w->lid.rays[idx].inv_dir = cast_vec3((float)(1.0/dx), (float)(1.0/dy), (float)(1.0/dz));
			idx++;
		}
	}
	w->lid.ray_count = rc;
	mat4x4_identity(&w->lid.transform);
}

World *world_create(int n_az, int n_el)
{
	World *w = calloc(1, sizeof(World));
	scene_init(&w->scn);
	w->scn.root = NULL;
	w->bmin[0] = w->bmin[1] = w->bmin[2] =  1e30;
	w->bmax[0] = w->bmax[1] = w->bmax[2] = -1e30;

	build_scene(w);

	/* snapshot each object's world transform (restored before every scan). */
	w->n_obj = w->scn.current_size;
	w->world_tf = malloc((size_t)w->n_obj * sizeof(mat4x4));
	for (int i = 0; i < w->n_obj; i++)
		w->world_tf[i] = w->scn.objects[i].transform;

	build_lidar(w, n_az, n_el);
	return w;
}

void world_set_lidar(World *w, int n_az, int n_el)
{
	build_lidar(w, n_az, n_el);
}

void world_destroy(World *w)
{
	if (!w) return;
	scene_free(&w->scn);        /* frees objects + BVH root */
	lidar_free(&w->lid);        /* frees rays                */
	free(w->world_tf);
	free(w);
}

void world_scan(World *w, const Pose *sensor, ScanCloud *out)
{
	/* restore world transforms (scene_update mutates them in place) */
	for (int i = 0; i < w->n_obj; i++)
		w->scn.objects[i].transform = w->world_tf[i];

	w->lid.transform = pose_to_mat4x4(sensor);
	scene_update(&w->scn, w->lid);     /* move scene into the sensor frame */

	if (w->scn.root) {                 /* scene_build leaks the old root otherwise */
		bvh_free_tree(w->scn.root);
		w->scn.root = NULL;
	}
	scene_build(&w->scn);

	pointcloud pc = cast_rays(w->lid, w->scn, 0, w->lid.ray_count);

	out->n = 0;
	scancloud_reserve(out, (int)pc.point_count);
	for (uint32_t k = 0; k < pc.point_count; k++)
		scancloud_push(out, pc.points[k].x, pc.points[k].y, pc.points[k].z);

	pointcloud_free(&pc);

	/* scene_update rebuilds each mesh's BVH without freeing the previous one;
	 * release them here so repeated scans don't leak. */
	for (int i = 0; i < w->n_obj; i++)
		if (w->scn.objects[i].type == MESH && w->scn.objects[i].mesh.root) {
			bvhmesh_free(w->scn.objects[i].mesh.root);
			w->scn.objects[i].mesh.root = NULL;
		}
}

void world_bounds(const World *w, double min[3], double max[3])
{
	for (int k = 0; k < 3; k++) { min[k] = w->bmin[k]; max[k] = w->bmax[k]; }
}
