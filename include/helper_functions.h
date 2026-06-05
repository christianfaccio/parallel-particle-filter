#ifndef HELPER_FUNCTIONS_H_
#define HELPER_FUNCTIONS_H_

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "map.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double velocity; double yawrate; } control_s;       /* one control measurement */
typedef struct { double x, y, theta; } ground_truth;                 /* one ground-truth pose   */
typedef struct { int id; double x, y; } LandmarkObs;                 /* one landmark observation */

/* 
 * Euclidean distance between two 2D points. 
 */
static inline double dist(double x1, double y1, double x2, double y2) 
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

/* 
 * Absolute (x, y, theta) error, theta wrapped into [0, pi].
 * Reentrant: writes into caller-provided error[3] instead of a static buffer. 
 */
static inline void getError(double gt_x, double gt_y, double gt_theta,
				double pf_x, double pf_y, double pf_theta,
                		double error[3]) 
{
	error[0] = fabs(pf_x - gt_x);
      	error[1] = fabs(pf_y - gt_y);
      	error[2] = fabs(pf_theta - gt_theta);
      	error[2] = fmod(error[2], 2.0 * M_PI);
      	if (error[2] > M_PI) error[2] = 2.0 * M_PI - error[2];
}

/*
 * Distribution functions.
 */
static inline double rand_normal(double mean, double stddev)
{
	double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
	double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
	double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
	return mean + stddev * z;
}

/* 
 * File readers. Each returns true on success.
 * The *_data readers allocate the output array with malloc; the caller frees it.
 * read_map_data fills map->landmark_list (free with free(map.landmark_list)). 
 */
bool read_map_data(const char *filename, Map *map);
bool read_control_data(const char *filename, control_s **position_meas, int *count);
bool read_gt_data(const char *filename, ground_truth **gt, int *count);
bool read_landmark_data(const char *filename, LandmarkObs **observations, int *count);

#endif /* HELPER_FUNCTIONS_H_ */
