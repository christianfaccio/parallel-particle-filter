#ifndef HELPER_FUNCTIONS_H_
#define HELPER_FUNCTIONS_H_

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Gaussian sample ~ N(mean, stddev) via Box-Muller.
 *
 * static inline so it can live in the header without multiple-definition
 * errors across translation units. NOTE: uses rand(), which is NOT thread
 * safe — the parallel backends must replace this with a counter-based,
 * per-particle RNG before going multi-threaded.
 */
static inline double rand_normal(double mean, double stddev)
{
	double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
	double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
	double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
	return mean + stddev * z;
}

#endif /* HELPER_FUNCTIONS_H_ */
