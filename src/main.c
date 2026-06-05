/*
 * Main function.
 * Reads data and runs the 2D particle filter.
 *
 * C port of the original main.cpp (Udacity "kidnapped vehicle" project).
 * Differences from the C++ reference:
 *   - std::vector<T>            -> T* + int count (allocated by the readers)
 *   - default_random_engine /  -> rand_normal() (Box-Muller on rand())
 *     normal_distribution<>
 *   - ostringstream filename    -> snprintf()
 *   - getError() returns        -> getError() writes into a caller-owned error[3]
 *     a static buffer
 */

#include "helper_functions.h"
#include "particle_filter.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* The original sets this inside ParticleFilter::init; our pf_init takes it as
 * an argument, so we keep the same value here. */
#define NUM_PARTICLES 729

int main(void)
{
	// Parameters related to grading.
	int    time_steps_before_lock_required = 100;  // steps before accuracy is checked.
	double max_runtime           = 45;             // Max allowable runtime to pass [sec]
	double max_translation_error = 1;              // Max allowable translation error [m]
	double max_yaw_error         = 0.05;           // Max allowable yaw error [rad]

	// Start timer.
	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);

	// Set up parameters here.
	double delta_t      = 0.1;  // Time elapsed between measurements [sec]
	double sensor_range = 50;   // Sensor range [m]

	/*
	 * Sigmas - just an estimate, usually comes from sensor uncertainty.
	 */
	double sigma_pos[3]      = {0.3, 0.3, 0.01}; // GPS uncertainty [x [m], y [m], theta [rad]]
	double sigma_landmark[2] = {0.3, 0.3};       // Landmark uncertainty [x [m], y [m]]

	// Read map data.
	Map map;
	if (!read_map_data("data/map_data.txt", &map)) {
		printf("Error: Could not open map file\n");
		return -1;
	}

	// Read control/position data (the readers allocate the arrays).
	control_s *position_meas = NULL;
	int count_control_data   = 0;
	if (!read_control_data("data/control_data.txt", &position_meas, &count_control_data)) {
		printf("Error: Could not open position/control measurement file\n");
		return -1;
	}

	// Read ground truth data.
	ground_truth *gt   = NULL;
	int count_gt_data  = 0;
	if (!read_gt_data("data/gt_data.txt", &gt, &count_gt_data)) {
		printf("Error: Could not open ground truth data file\n");
		return -1;
	}

	// Run the particle filter.
	int num_time_steps = count_control_data;
	ParticleFilter pf = {0};
	double total_error[3]    = {0, 0, 0};
	double cum_mean_error[3] = {0, 0, 0};

	for (int i = 0; i < num_time_steps; i++)
	{
		printf("Time step: %i\n", i);

		// Read landmark observations for the current time step.
		char filename[128];
		snprintf(filename, sizeof filename,
		         "data/observation/observations_%06d.txt", i + 1);

		LandmarkObs *observations = NULL;
		int n_observations        = 0;
		if (!read_landmark_data(filename, &observations, &n_observations)) {
			printf("Error: Could not open observation file %d\n", i + 1);
			return -1;
		}

		// Initialize the filter on the first time step, otherwise predict forward.
		if (!pf_initialized(&pf)) {
			double n_x     = rand_normal(0.0, sigma_pos[0]);
			double n_y     = rand_normal(0.0, sigma_pos[1]);
			double n_theta = rand_normal(0.0, sigma_pos[2]);
			pf_init(&pf, NUM_PARTICLES,
			        gt[i].x + n_x, gt[i].y + n_y, gt[i].theta + n_theta,
			        sigma_pos);
		} else {
			pf_prediction(&pf, delta_t, sigma_pos,
			              position_meas[i - 1].velocity,
			              position_meas[i - 1].yawrate);
		}

		// Simulate sensor noise on the (noiseless) observation data, in place.
		for (int j = 0; j < n_observations; j++) {
			observations[j].x += rand_normal(0.0, sigma_landmark[0]);
			observations[j].y += rand_normal(0.0, sigma_landmark[1]);
		}

		// Update the weights and resample.
		pf_update_weights(&pf, sensor_range, sigma_landmark,
		                  observations, n_observations, &map);
		pf_resample(&pf);

		// Find the highest-weight particle as the filter's best estimate.
		double highest_weight = 0.0;
		Particle best_particle = pf.particles[0];
		for (int b = 0; b < pf.num_particles; b++) {
			if (pf.particles[b].weight > highest_weight) {
				highest_weight = pf.particles[b].weight;
				best_particle  = pf.particles[b];
			}
		}

		double avg_error[3];
		getError(gt[i].x, gt[i].y, gt[i].theta,
		         best_particle.x, best_particle.y, best_particle.theta,
		         avg_error);

		for (int j = 0; j < 3; j++) {
			total_error[j]    += avg_error[j];
			cum_mean_error[j]  = total_error[j] / (double)(i + 1);
		}

		// Print the cumulative weighted error.
		printf("Cumulative mean weighted error: x %f y %f yaw %f\n",
		       cum_mean_error[0], cum_mean_error[1], cum_mean_error[2]);

		// If the error is too high after the lock-in point, report and exit.
		if (i >= time_steps_before_lock_required) {
			if (cum_mean_error[0] > max_translation_error ||
			    cum_mean_error[1] > max_translation_error ||
			    cum_mean_error[2] > max_yaw_error) {
				if (cum_mean_error[0] > max_translation_error) {
					printf("Your x error, %f is larger than the maximum allowable error, %f\n",
					       cum_mean_error[0], max_translation_error);
				} else if (cum_mean_error[1] > max_translation_error) {
					printf("Your y error, %f is larger than the maximum allowable error, %f\n",
					       cum_mean_error[1], max_translation_error);
				} else {
					printf("Your yaw error, %f is larger than the maximum allowable error, %f\n",
					       cum_mean_error[2], max_yaw_error);
				}
				free(observations);
				free(position_meas);
				free(gt);
				free(map.landmark_list);
				pf_free(&pf);
				return -1;
			}
		}

		free(observations);
	}

	// End timer and report runtime.
	clock_gettime(CLOCK_MONOTONIC, &end);
	double runtime = (end.tv_sec - start.tv_sec) +
	                 (end.tv_nsec - start.tv_nsec) / 1e9;  // seconds
	printf("Runtime (sec): %f\n", runtime);

	if (runtime < max_runtime && pf_initialized(&pf)) {
		printf("Success! Your particle filter passed!\n");
	} else if (!pf_initialized(&pf)) {
		printf("This is the starter code. You haven't initialized your filter.\n");
	} else {
		printf("Your runtime %f is larger than the maximum allowable runtime, %f\n",
		       runtime, max_runtime);
		free(position_meas);
		free(gt);
		free(map.landmark_list);
		pf_free(&pf);
		return -1;
	}

	// Clean up.
	free(position_meas);
	free(gt);
	free(map.landmark_list);
	pf_free(&pf);

	return 0;
}
