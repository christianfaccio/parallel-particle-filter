#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

typedef struct {
	int id;
	double x;
	double y;
	double theta;
	double weight;
} Particle;

typedef struct {
	int num_particles;
	Particle *particles;	// dynamically grown array
	bool is_initialized;
	double *weights;
} ParticleFilter;

/*
 * Initialize: allocate arrays, draw particles ~ N(mean, std) around (x,y,theta),
 * set all weights to 1. std[] = {std_x, std_y, std_theta}.
 */
void pf_init(ParticleFilter *self, int num_particles, double x, double y, double theta, double std[]);

/*
 * Free the arrays owned by the filter.
 */
void pf_free(ParticleFilter *self);

/*
 * Predict each particle forward with the CTRV (velocity + yaw-rate) motion model.
 */
void pf_prediction(ParticleFilter *self, double delta_t, double std_pos[], double velocity, double yaw_rate);

/*
 * Nearest-neighbour data association: assigns each observation the id of the 
 * closest predicted landmark. Modifies observations in place (No filter state -> no self).
 */
void pf_data_association(LandmarkObs *predicted, int n_predicted, LandmarkObs *observations, int n_observations);

/*
 * Update weights from measurement likelihood. std_landmark[] = {std_range, std_bearing}.
 */
void pf_update_weights(ParticleFilter *self, double sensor_range, double std_landmark[],
		LandmarkObs *observations, int n_observations, const Map *map_landmarks);

/*
 * Resample particles proportional to weight.
 */
void pf_resample(ParticleFilter *self);

/*
 * Append particle poses to a file.
 */
void pf_write(const ParticleFilter *self, const char *filename);

/*
 * Whether pf_init has run.
 */
bool pf_initialized(const ParticleFilter *self);

#endif /* PARTICLE_FILTER_H_ */
