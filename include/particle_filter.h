#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <stdbool.h>

#include "pose.h"
#include "world.h"       /* ScanCloud */
#include "lidar_map.h"   /* DistanceField */

/*
 * 6-DOF Monte-Carlo Localization with a LiDAR likelihood-field measurement
 * model. Each particle is a full rigid-body pose hypothesis; it is scored by
 * how well the live scan, re-projected through that pose, lands on mapped
 * surfaces (a distance-field lookup per scan point).
 */

typedef struct {
	Pose   pose;
	double weight;
} Particle;

typedef struct {
	int       num_particles;
	Particle *particles;
	double   *weights;       /* normalized copy, used by resampling */
	bool      is_initialized;
} ParticleFilter;

/*
 * Allocate and sample `n` particles ~ N(mean, std) in 6-DOF: position spread
 * pos_std[3] (world axes, m) and small-angle orientation spread rot_std[3]
 * (body axes, rad). Weights start uniform.
 */
void pf_init(ParticleFilter *self, int n, Pose mean,
             const double pos_std[3], const double rot_std[3]);

/* Free the arrays owned by the filter. */
void pf_free(ParticleFilter *self);

/*
 * Predict: compose the (noisy) body-frame relative motion `rel` onto every
 * particle. pos_std[3]/rot_std[3] are the process-noise sigmas.
 */
void pf_prediction(ParticleFilter *self, Pose rel,
                   const double pos_std[3], const double rot_std[3]);

/*
 * Update weights from the LiDAR likelihood field. For each particle and each
 * scan point: transform the point into the map frame by the particle pose,
 * look up the distance to the nearest surface, and accumulate a Gaussian
 * log-likelihood with std `sigma` [m]. Weights are normalized in place.
 *
 * This N_particles x M_points loop is the HPC hot kernel.
 */
void pf_update_weights(ParticleFilter *self, const ScanCloud *scan,
                       const DistanceField *df, double sigma);

/* Resample particles proportional to weight (multinomial / CDF scan). */
void pf_resample(ParticleFilter *self);

/* Filter estimate: weighted-mean position and weighted-average orientation. */
void pf_estimate(const ParticleFilter *self, Pose *out);

/* Whether pf_init has run. */
bool pf_initialized(const ParticleFilter *self);

#endif /* PARTICLE_FILTER_H_ */
