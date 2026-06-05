/*
 * main.c — 6-DOF LiDAR Monte-Carlo Localization driver.
 *
 * Everything is simulated, so ground-truth poses are known exactly:
 *   1. build a ground-truth world (lidarsim primitives) and a LiDAR pattern;
 *   2. build the localization map (reference cloud -> 3D distance field) once;
 *   3. fly a scripted 6-DOF trajectory; at each step simulate a noisy scan and
 *      run predict -> update -> resample -> estimate;
 *   4. report the trajectory error (ATE) against ground truth.
 *
 * Usage: pf_serial [num_particles] [num_steps] [n_az] [n_el]
 */

#include "pose.h"
#include "world.h"
#include "lidar_map.h"
#include "particle_filter.h"
#include "helper_functions.h"   /* rand_normal */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Scripted ground-truth UAV pose at step t of T: a tilted circular orbit with
 * z bob, heading-aligned yaw, and gentle roll/pitch — exercises all 6 DOF. */
static Pose gt_pose(int t, int T)
{
	double theta = 2.0 * M_PI * 1.5 * t / T;   /* 1.5 loops over the run */
	Pose p;
	p.pos[0] = 4.0 * cos(theta);
	p.pos[1] = 4.0 * sin(theta);
	p.pos[2] = 1.5 + 0.4 * sin(2.0 * theta);
	double yaw   = theta + M_PI / 2.0;          /* tangent to the circle */
	double pitch = 0.10 * sin(theta);
	double roll  = 0.10 * cos(theta);
	euler_to_quat(roll, pitch, yaw, p.quat);
	return p;
}

static double pos_dist(const double a[3], const double b[3])
{
	double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

int main(int argc, char **argv)
{
	int    N      = (argc > 1) ? atoi(argv[1]) : 3000;   /* particles */
	int    T      = (argc > 2) ? atoi(argv[2]) : 200;    /* time steps */
	int    n_az   = (argc > 3) ? atoi(argv[3]) : 180;    /* LiDAR azimuths */
	int    n_el   = (argc > 4) ? atoi(argv[4]) : 16;     /* LiDAR elevation rings */

	/* noise / model parameters */
	const double sigma_scan = 0.02;   /* per-point sensor noise [m]         */
	const double sigma_meas = 0.20;   /* likelihood-field Gaussian std [m]  */
	const double init_pos_std[3] = {0.40, 0.40, 0.40};
	const double init_rot_std[3] = {0.05, 0.05, 0.05};   /* ~3 deg */
	const double proc_pos_std[3] = {0.03, 0.03, 0.03};
	const double proc_rot_std[3] = {0.010, 0.010, 0.010};
	const double odom_pos_std[3] = {0.02, 0.02, 0.02};
	const double odom_rot_std[3] = {0.005, 0.005, 0.005};

	/* pass/fail thresholds */
	const double max_trans_err = 0.5;             /* m   */
	const double max_rot_err   = 5.0 * M_PI/180;  /* rad */
	const int    warmup        = 10;              /* steps excluded from the mean */

	printf("6-DOF LiDAR MCL: N=%d particles, T=%d steps, LiDAR=%dx%d rays\n",
	       N, T, n_az, n_el);

	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC, &t0);

	/* --- world + map (built once) --- */
	World *w = world_create(n_az, n_el);

	/* Build the map with a DENSE LiDAR + many viewpoints so the distance field
	 * approximates true distance-to-surface, then restore the runtime pattern. */
	ScanCloud ref;
	scancloud_init(&ref);
	world_set_lidar(w, 720, 180);
	lidar_map_build_cloud(w, 5, &ref);
	world_set_lidar(w, n_az, n_el);

	DistanceField df;
	if (df_build(&df, &ref, /*voxel=*/0.15, /*pad=*/1.0) != 0) {
		fprintf(stderr, "Error: distance-field build failed\n");
		return -1;
	}
	printf("Map: reference cloud %d pts, distance field %dx%dx%d (%.2f m voxels)\n",
	       ref.n, df.dim[0], df.dim[1], df.dim[2], df.voxel);

	struct timespec tmap;
	clock_gettime(CLOCK_MONOTONIC, &tmap);

	/* --- run the filter along the trajectory --- */
	ParticleFilter pf = {0};
	ScanCloud scan;
	scancloud_init(&scan);

	double sum_trans = 0.0, sum_rot = 0.0;
	int    counted = 0;

	/* Optional point-cloud dump for a viewer (set PF_DUMP=1). We collect the
	 * ground-truth and estimated trajectories over the whole run, plus a single
	 * world-frame scan and the particle cloud at one "snapshot" step. */
	const char *dump = getenv("PF_DUMP");
	const int   viz_step = T / 4;
	ScanCloud gt_traj, est_traj, scan_world, particles_pc;
	scancloud_init(&gt_traj); scancloud_init(&est_traj);
	scancloud_init(&scan_world); scancloud_init(&particles_pc);

	for (int t = 0; t < T; t++) {
		Pose gt = gt_pose(t, T);

		/* simulate the live scan (sensor frame) + per-point noise */
		world_scan(w, &gt, &scan);
		for (int m = 0; m < scan.n; m++) {
			scan.x[m] += (float)rand_normal(0.0, sigma_scan);
			scan.y[m] += (float)rand_normal(0.0, sigma_scan);
			scan.z[m] += (float)rand_normal(0.0, sigma_scan);
		}

		if (!pf_initialized(&pf)) {
			pf_init(&pf, N, gt, init_pos_std, init_rot_std);
		} else {
			Pose prev = gt_pose(t - 1, T);
			Pose rel  = pose_relative(&prev, &gt);     /* control = GT motion */
			pose_perturb(&rel, odom_pos_std, odom_rot_std); /* noisy odometry */
			pf_prediction(&pf, rel, proc_pos_std, proc_rot_std);
		}

		pf_update_weights(&pf, &scan, &df, sigma_meas);
		pf_resample(&pf);

		Pose est;
		pf_estimate(&pf, &est);

		double terr = pos_dist(est.pos, gt.pos);
		double rerr = quat_geodesic(est.quat, gt.quat);
		if (t >= warmup) { sum_trans += terr; sum_rot += rerr; counted++; }

		if (dump) {
			scancloud_push(&gt_traj,  (float)gt.pos[0],  (float)gt.pos[1],  (float)gt.pos[2]);
			scancloud_push(&est_traj, (float)est.pos[0], (float)est.pos[1], (float)est.pos[2]);
			if (t == viz_step) {
				/* the live scan re-projected through GT -> world frame (overlays the map) */
				for (int m = 0; m < scan.n; m++) {
					double v[3] = { scan.x[m], scan.y[m], scan.z[m] }, wpt[3];
					pose_transform_point(&gt, v, wpt);
					scancloud_push(&scan_world, (float)wpt[0], (float)wpt[1], (float)wpt[2]);
				}
				for (int i = 0; i < pf.num_particles; i++)
					scancloud_push(&particles_pc,
					               (float)pf.particles[i].pose.pos[0],
					               (float)pf.particles[i].pose.pos[1],
					               (float)pf.particles[i].pose.pos[2]);
			}
		}

		if (t % 20 == 0 || t == T - 1)
			printf("  step %3d: scan %4d pts | trans err %.3f m | rot err %.2f deg\n",
			       t, scan.n, terr, rerr * 180.0 / M_PI);
	}

	clock_gettime(CLOCK_MONOTONIC, &t1);

	double mean_trans = counted ? sum_trans / counted : 0.0;
	double mean_rot   = counted ? sum_rot   / counted : 0.0;
	double map_s = (tmap.tv_sec - t0.tv_sec) + (tmap.tv_nsec - t0.tv_nsec)/1e9;
	double run_s = (t1.tv_sec - tmap.tv_sec) + (t1.tv_nsec - tmap.tv_nsec)/1e9;

	printf("\nMean ATE (after %d-step warmup): trans %.3f m | rot %.2f deg\n",
	       warmup, mean_trans, mean_rot * 180.0 / M_PI);
	printf("Timing: map build %.2f s | filter %.2f s (%.1f ms/step)\n",
	       map_s, run_s, 1000.0 * run_s / T);

	int ok = (mean_trans < max_trans_err) && (mean_rot < max_rot_err);
	printf("%s\n", ok ? "Success! Localization converged within thresholds."
	                   : "FAIL: error above thresholds.");

	if (dump) {
		int map_stride = ref.n / 500000 + 1;   /* keep map.pcd ~<=500k points */
		scancloud_save_pcd(&ref,          "map.pcd",        map_stride);
		scancloud_save_pcd(&scan_world,   "scan.pcd",       1);
		scancloud_save_pcd(&particles_pc, "particles.pcd",  1);
		scancloud_save_pcd(&gt_traj,      "gt_traj.pcd",    1);
		scancloud_save_pcd(&est_traj,     "est_traj.pcd",   1);
		printf("\nWrote PCDs (step %d snapshot): map.pcd (every %dth pt), "
		       "scan.pcd, particles.pcd, gt_traj.pcd, est_traj.pcd\n",
		       viz_step, map_stride);
	}

	/* cleanup */
	pf_free(&pf);
	scancloud_free(&scan);
	scancloud_free(&ref);
	scancloud_free(&gt_traj);
	scancloud_free(&est_traj);
	scancloud_free(&scan_world);
	scancloud_free(&particles_pc);
	df_free(&df);
	world_destroy(w);

	return ok ? 0 : 1;
}
