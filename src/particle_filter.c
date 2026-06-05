/*
 * Particle Filter source code
 */

#include "particle_filter.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

void pf_init(ParticleFilter *self, int num_particles, double x, double y, double theta, double std[])
{
	double std_x     = std[0];
	double std_y     = std[1];
	double std_theta = std[2];

	// Allocate memory
	Particle *particles = malloc(num_particles * sizeof(Particle));
	double   *weights   = malloc(num_particles * sizeof(double));

	// Seed once for the whole run (reproducibility)
	srand(42);

	// Create particles and set values
	for (int i = 0; i < num_particles; i++)
	{
		Particle p;
		p.id     = i;
		p.x      = rand_normal(x, std_x);
		p.y      = rand_normal(y, std_y);
		p.theta  = rand_normal(theta, std_theta);
		p.weight = 1.0;

		particles[i] = p;
		weights[i]   = p.weight;
	}

	self->num_particles  = num_particles;
	self->particles      = particles;
	self->weights        = weights;
	self->is_initialized = true;
}

void pf_free(ParticleFilter *self)
{
	free(self->particles);
	free(self->weights);
	self->particles = NULL;
	self->weights   = NULL;
}

void pf_prediction(ParticleFilter *self, double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	int num_particles = self->num_particles;

	for (int i = 0; i < num_particles; i++)
	{
		Particle *p = &self->particles[i];

		// CTRV motion model (assumes yaw_rate != 0)
		double new_x     = p->x + (velocity / yaw_rate) * (sin(p->theta + yaw_rate * delta_t) - sin(p->theta));
		double new_y     = p->y + (velocity / yaw_rate) * (cos(p->theta) - cos(p->theta + yaw_rate * delta_t));
		double new_theta = p->theta + (yaw_rate * delta_t);

		// Add Gaussian process noise around the predicted mean
		p->x     = rand_normal(new_x, std_pos[0]);
		p->y     = rand_normal(new_y, std_pos[1]);
		p->theta = rand_normal(new_theta, std_pos[2]);
	}
}

void pf_data_association(LandmarkObs *predicted, int n_predicted,
                         LandmarkObs *observations, int n_observations)
{
	// For each observation, find the nearest predicted landmark and adopt its id.
	for (int j = 0; j < n_observations; j++)
	{
		double dist_min = DBL_MAX;
		for (int i = 0; i < n_predicted; i++)
		{
			double distance = dist(observations[j].x, observations[j].y,
			                       predicted[i].x, predicted[i].y);
			if (distance < dist_min)
			{
				dist_min = distance;
				observations[j].id = predicted[i].id;
			}
		}
	}
}

void pf_update_weights(ParticleFilter *self, double sensor_range, double std_landmark[],
                       LandmarkObs *observations, int n_observations, const Map *map_landmarks)
{
	(void)sensor_range;  // unused, kept for interface parity with the original

	double std_x       = std_landmark[0];
	double std_y       = std_landmark[1];
	double weights_sum = 0.0;
	int    num_particles = self->num_particles;

	for (int i = 0; i < num_particles; i++)
	{
		Particle *p = &self->particles[i];
		double wt = 1.0;

		// Convert each observation from the vehicle frame to the map frame
		for (int j = 0; j < n_observations; j++)
		{
			LandmarkObs current_obs = observations[j];
			LandmarkObs transformed_obs;

			transformed_obs.x  = (current_obs.x * cos(p->theta)) - (current_obs.y * sin(p->theta)) + p->x;
			transformed_obs.y  = (current_obs.x * sin(p->theta)) + (current_obs.y * cos(p->theta)) + p->y;
			transformed_obs.id = current_obs.id;

			// Nearest-landmark association
			single_landmark_s landmark = {0};
			double distance_min = DBL_MAX;

			for (int k = 0; k < map_landmarks->count; k++)
			{
				single_landmark_s cur_l = map_landmarks->landmark_list[k];
				double distance = dist(transformed_obs.x, transformed_obs.y, cur_l.x_f, cur_l.y_f);
				if (distance < distance_min)
				{
					distance_min = distance;
					landmark = cur_l;
				}
			}

			// Multivariate Gaussian likelihood
			double num   = exp(-0.5 * (pow(transformed_obs.x - landmark.x_f, 2) / pow(std_x, 2)
			                         + pow(transformed_obs.y - landmark.y_f, 2) / pow(std_y, 2)));
			double denom = 2.0 * M_PI * std_x * std_y;
			wt *= num / denom;
		}

		weights_sum += wt;
		p->weight = wt;
	}

	// Normalize weights into (0, 1]
	for (int i = 0; i < num_particles; i++)
	{
		Particle *p = &self->particles[i];
		p->weight /= weights_sum;
		self->weights[i] = p->weight;
	}
}

void pf_resample(ParticleFilter *self)
{
	int n = self->num_particles;

	// Build the cumulative distribution of the (normalized) weights.
	// cumsum is ~1.0 if pf_update_weights normalized, but we use it explicitly
	// so resampling is also correct on un-normalized weights.
	double *cdf = malloc(n * sizeof(double));
	double cumsum = 0.0;
	for (int i = 0; i < n; i++)
	{
		cumsum += self->weights[i];
		cdf[i] = cumsum;
	}

	Particle *resampled = malloc(n * sizeof(Particle));

	// Multinomial resampling with replacement: draw n indices, each with
	// probability proportional to its weight (equivalent to std::discrete_distribution).
	for (int i = 0; i < n; i++)
	{
		double u = ((double)rand() / RAND_MAX) * cumsum;

		// First particle whose cdf >= u. Linear scan here; this O(N) search is
		// exactly what becomes a parallel scan / binary search in the HPC backends.
		int idx = 0;
		while (idx < n - 1 && cdf[idx] < u)
			idx++;

		resampled[i] = self->particles[idx];
	}

	free(cdf);
	free(self->particles);
	self->particles = resampled;
}

void pf_write(const ParticleFilter *self, const char *filename)
{
	FILE *f = fopen(filename, "a");   // append (std::ios::app)
	if (!f)
		return;

	for (int i = 0; i < self->num_particles; i++)
	{
		fprintf(f, "%f %f %f\n",
		        self->particles[i].x,
		        self->particles[i].y,
		        self->particles[i].theta);
	}

	fclose(f);
}

bool pf_initialized(const ParticleFilter *self)
{
	return self->is_initialized;
}
