#ifndef POSE_H_
#define POSE_H_

/*
 * 6-DOF rigid-body pose: position + orientation.
 *
 * Orientation is stored as a unit quaternion (w, x, y, z) — no gimbal lock,
 * clean composition, and process noise is applied as a small-angle axis-angle
 * perturbation (see pose_perturb). Euler angles (roll, pitch, yaw) are used
 * only for input convenience and error reporting.
 *
 * This module is pure math and has NO dependency on lidarsim (so it can be
 * included anywhere). The conversion Pose -> lidarsim mat4x4 lives in world.c,
 * the only translation unit that knows about lidarsim types.
 */

typedef struct {
	double pos[3];   /* x, y, z in the world frame            */
	double quat[4];  /* unit quaternion w, x, y, z (body->world) */
} Pose;

/* ---- quaternion primitives (q = [w, x, y, z]) ---- */

/* Identity orientation at the origin. */
Pose pose_identity(void);

/* Normalize q in place to unit length. */
void quat_normalize(double q[4]);

/* out = a * b  (Hamilton product; rotation composition). */
void quat_mul(const double a[4], const double b[4], double out[4]);

/* Quaternion from an axis-angle rotation vector rv (= axis * angle [rad]). */
void quat_from_rotvec(const double rv[3], double q[4]);

/* roll/pitch/yaw (XYZ intrinsic) -> quaternion and back. */
void euler_to_quat(double roll, double pitch, double yaw, double q[4]);
void quat_to_euler(const double q[4], double *roll, double *pitch, double *yaw);

/* Rotate vector v by quaternion q: out = q * v * q^-1. */
void quat_rotate(const double q[4], const double v[3], double out[3]);

/* Row-major 3x3 rotation matrix R[9] equivalent to q (body->world). */
void quat_rotation_matrix(const double q[4], double R[9]);

/* Geodesic angle [rad] between two orientations (>= 0). */
double quat_geodesic(const double a[4], const double b[4]);

/* ---- pose operations ---- */

/* Apply a body-frame relative motion `rel` to `base` (base then rel):
 *   out.pos  = base.pos + R(base) * rel.pos
 *   out.quat = base.quat * rel.quat
 */
Pose pose_compose(const Pose *base, const Pose *rel);

/* Relative motion taking `from` to `to`, expressed in `from`'s body frame:
 *   rel = from^-1 * to   (inverse of pose_compose).
 */
Pose pose_relative(const Pose *from, const Pose *to);

/* Add zero-mean Gaussian noise: position std pos_std[3] (world axes) and a
 * small-angle rotation std rot_std[3] (body axes, rad). Uses rand_normal. */
void pose_perturb(Pose *p, const double pos_std[3], const double rot_std[3]);

/* Transform a point from body frame to world frame: out = R(p)*v + p.pos. */
void pose_transform_point(const Pose *p, const double v[3], double out[3]);

#endif /* POSE_H_ */
