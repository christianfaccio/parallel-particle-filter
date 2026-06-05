/*
 * particle_filter.c — 6-DOF LiDAR Monte-Carlo Localization.
 */

#include "particle_filter.h"
#include "helper_functions.h"   /* rand_normal */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void pf_init(ParticleFilter *self, int n, Pose mean,
             const double pos_std[3], const double rot_std[3])
{
	self->num_particles = n;
	self->particles = malloc((size_t)n * sizeof(Particle));
	self->weights   = malloc((size_t)n * sizeof(double));

	srand(42);   /* reproducible baseline (single-threaded RNG) */

	for (int i = 0; i < n; i++) {
		Particle p;
		p.pose = mean;
		pose_perturb(&p.pose, pos_std, rot_std);
		p.weight = 1.0 / n;
		self->particles[i] = p;
		self->weights[i]   = p.weight;
	}
	self->is_initialized = true;
}

void pf_free(ParticleFilter *self)
{
	free(self->particles);
	free(self->weights);
	self->particles = NULL;
	self->weights   = NULL;
}

void pf_prediction(ParticleFilter *self, Pose rel,
                   const double pos_std[3], const double rot_std[3])
{
	for (int i = 0; i < self->num_particles; i++) {
		Particle *p = &self->particles[i];
		p->pose = pose_compose(&p->pose, &rel);   /* deterministic motion */
		pose_perturb(&p->pose, pos_std, rot_std); /* process noise        */
	}
}

void pf_update_weights(ParticleFilter *self, const ScanCloud *scan,
                       const DistanceField *df, double sigma)
{
	const int    N = self->num_particles;
	const int    M = scan->n;
	const double inv_two_sigma2 = 1.0 / (2.0 * sigma * sigma);

	double max_log = -INFINITY;

	/* --- hot kernel: N particles x M scan points --- */
	for (int i = 0; i < N; i++) {
		Particle *p = &self->particles[i];

		double R[9];
		quat_rotation_matrix(p->pose.quat, R);
		const double tx = p->pose.pos[0];
		const double ty = p->pose.pos[1];
		const double tz = p->pose.pos[2];

		double logw = 0.0;
		for (int m = 0; m < M; m++) {
			double vx = scan->x[m], vy = scan->y[m], vz = scan->z[m];
			double wx = R[0]*vx + R[1]*vy + R[2]*vz + tx;
			double wy = R[3]*vx + R[4]*vy + R[5]*vz + ty;
			double wz = R[6]*vx + R[7]*vy + R[8]*vz + tz;
			double d  = df_lookup(df, wx, wy, wz);
			logw -= d * d * inv_two_sigma2;
		}
		p->weight = logw;                 /* stash log-weight */
		if (logw > max_log) max_log = logw;
	}

	/* --- log-sum-exp normalization (max-shift for stability) --- */
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		double w = exp(self->particles[i].weight - max_log);
		self->particles[i].weight = w;
		sum += w;
	}
	double inv_sum = (sum > 0.0) ? 1.0 / sum : 0.0;
	for (int i = 0; i < N; i++) {
		self->particles[i].weight *= inv_sum;
		self->weights[i] = self->particles[i].weight;
	}
}

void pf_resample(ParticleFilter *self)
{
	int n = self->num_particles;

	/* cumulative distribution of the (normalized) weights */
	double *cdf = malloc((size_t)n * sizeof(double));
	double cumsum = 0.0;
	for (int i = 0; i < n; i++) {
		cumsum += self->weights[i];
		cdf[i] = cumsum;
	}

	Particle *resampled = malloc((size_t)n * sizeof(Particle));

	/* multinomial resampling with replacement; the O(N) linear scan here is
	 * exactly what becomes a parallel scan / binary search in the HPC backends. */
	for (int i = 0; i < n; i++) {
		double u = ((double)rand() / RAND_MAX) * cumsum;
		int idx = 0;
		while (idx < n - 1 && cdf[idx] < u) idx++;
		resampled[i] = self->particles[idx];
	}

	free(cdf);
	free(self->particles);
	self->particles = resampled;
}

void pf_estimate(const ParticleFilter *self, Pose *out)
{
	int n = self->num_particles;

	/* weighted mean position; track the highest-weight particle for the
	 * quaternion hemisphere reference. */
	out->pos[0] = out->pos[1] = out->pos[2] = 0.0;
	double wsum = 0.0;
	int    best = 0;
	for (int i = 0; i < n; i++) {
		double w = self->particles[i].weight;
		out->pos[0] += w * self->particles[i].pose.pos[0];
		out->pos[1] += w * self->particles[i].pose.pos[1];
		out->pos[2] += w * self->particles[i].pose.pos[2];
		wsum += w;
		if (w > self->particles[best].weight) best = i;
	}
	if (wsum > 0.0) {
		out->pos[0] /= wsum; out->pos[1] /= wsum; out->pos[2] /= wsum;
	}

	/* weighted quaternion average: accumulate with hemisphere alignment to the
	 * highest-weight particle, then renormalize (valid for tight clusters). */
	const double *ref = self->particles[best].pose.quat;
	double acc[4] = {0, 0, 0, 0};
	for (int i = 0; i < n; i++) {
		const double *q = self->particles[i].pose.quat;
		double dot = q[0]*ref[0] + q[1]*ref[1] + q[2]*ref[2] + q[3]*ref[3];
		double s = (dot < 0.0) ? -self->particles[i].weight : self->particles[i].weight;
		acc[0] += s*q[0]; acc[1] += s*q[1]; acc[2] += s*q[2]; acc[3] += s*q[3];
	}
	out->quat[0] = acc[0]; out->quat[1] = acc[1];
	out->quat[2] = acc[2]; out->quat[3] = acc[3];
	quat_normalize(out->quat);
}

bool pf_initialized(const ParticleFilter *self)
{
	return self->is_initialized;
}
