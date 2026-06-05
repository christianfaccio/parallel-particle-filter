/*
 * 6-DOF pose / quaternion math. Pure C, no external dependencies beyond libm
 * and rand_normal (for pose_perturb).
 */

#include "pose.h"
#include "helper_functions.h"   /* rand_normal */

#include <math.h>

Pose pose_identity(void)
{
	Pose p = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0}};
	return p;
}

void quat_normalize(double q[4])
{
	double n = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
	if (n < 1e-12) {            /* degenerate -> identity */
		q[0] = 1.0; q[1] = q[2] = q[3] = 0.0;
		return;
	}
	double inv = 1.0 / n;
	q[0] *= inv; q[1] *= inv; q[2] *= inv; q[3] *= inv;
}

void quat_mul(const double a[4], const double b[4], double out[4])
{
	double w = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
	double x = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
	double y = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
	double z = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
	out[0] = w; out[1] = x; out[2] = y; out[3] = z;
}

void quat_from_rotvec(const double rv[3], double q[4])
{
	double angle = sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]);
	if (angle < 1e-12) {       /* small-angle -> identity */
		q[0] = 1.0; q[1] = q[2] = q[3] = 0.0;
		return;
	}
	double half = 0.5 * angle;
	double s = sin(half) / angle;   /* sin(half)/angle = (sin(half)/half)/2 */
	q[0] = cos(half);
	q[1] = rv[0] * s;
	q[2] = rv[1] * s;
	q[3] = rv[2] * s;
}

void euler_to_quat(double roll, double pitch, double yaw, double q[4])
{
	double cr = cos(roll  * 0.5), sr = sin(roll  * 0.5);
	double cp = cos(pitch * 0.5), sp = sin(pitch * 0.5);
	double cy = cos(yaw   * 0.5), sy = sin(yaw   * 0.5);
	q[0] = cr*cp*cy + sr*sp*sy;
	q[1] = sr*cp*cy - cr*sp*sy;
	q[2] = cr*sp*cy + sr*cp*sy;
	q[3] = cr*cp*sy - sr*sp*cy;
	quat_normalize(q);
}

void quat_to_euler(const double q[4], double *roll, double *pitch, double *yaw)
{
	double w = q[0], x = q[1], y = q[2], z = q[3];
	/* roll (x-axis) */
	double sinr_cosp = 2.0 * (w*x + y*z);
	double cosr_cosp = 1.0 - 2.0 * (x*x + y*y);
	*roll = atan2(sinr_cosp, cosr_cosp);
	/* pitch (y-axis), clamped at the poles */
	double sinp = 2.0 * (w*y - z*x);
	if (sinp >  1.0) sinp =  1.0;
	if (sinp < -1.0) sinp = -1.0;
	*pitch = asin(sinp);
	/* yaw (z-axis) */
	double siny_cosp = 2.0 * (w*z + x*y);
	double cosy_cosp = 1.0 - 2.0 * (y*y + z*z);
	*yaw = atan2(siny_cosp, cosy_cosp);
}

void quat_rotation_matrix(const double q[4], double R[9])
{
	double w = q[0], x = q[1], y = q[2], z = q[3];
	double xx = x*x, yy = y*y, zz = z*z;
	double xy = x*y, xz = x*z, yz = y*z;
	double wx = w*x, wy = w*y, wz = w*z;
	R[0] = 1.0 - 2.0*(yy + zz); R[1] = 2.0*(xy - wz);       R[2] = 2.0*(xz + wy);
	R[3] = 2.0*(xy + wz);       R[4] = 1.0 - 2.0*(xx + zz); R[5] = 2.0*(yz - wx);
	R[6] = 2.0*(xz - wy);       R[7] = 2.0*(yz + wx);       R[8] = 1.0 - 2.0*(xx + yy);
}

void quat_rotate(const double q[4], const double v[3], double out[3])
{
	double R[9];
	quat_rotation_matrix(q, R);
	out[0] = R[0]*v[0] + R[1]*v[1] + R[2]*v[2];
	out[1] = R[3]*v[0] + R[4]*v[1] + R[5]*v[2];
	out[2] = R[6]*v[0] + R[7]*v[1] + R[8]*v[2];
}

double quat_geodesic(const double a[4], const double b[4])
{
	/* angle = 2 * acos(|<a,b>|) for unit quaternions */
	double dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
	dot = fabs(dot);
	if (dot > 1.0) dot = 1.0;
	return 2.0 * acos(dot);
}

Pose pose_compose(const Pose *base, const Pose *rel)
{
	Pose out;
	double rp[3];
	quat_rotate(base->quat, rel->pos, rp);
	out.pos[0] = base->pos[0] + rp[0];
	out.pos[1] = base->pos[1] + rp[1];
	out.pos[2] = base->pos[2] + rp[2];
	quat_mul(base->quat, rel->quat, out.quat);
	quat_normalize(out.quat);
	return out;
}

Pose pose_relative(const Pose *from, const Pose *to)
{
	/* from^-1: conjugate quaternion (unit), negated rotated translation */
	double qinv[4] = { from->quat[0], -from->quat[1], -from->quat[2], -from->quat[3] };
	double dp[3] = { to->pos[0] - from->pos[0],
	                 to->pos[1] - from->pos[1],
	                 to->pos[2] - from->pos[2] };
	Pose rel;
	quat_rotate(qinv, dp, rel.pos);
	quat_mul(qinv, to->quat, rel.quat);
	quat_normalize(rel.quat);
	return rel;
}

void pose_perturb(Pose *p, const double pos_std[3], const double rot_std[3])
{
	p->pos[0] += rand_normal(0.0, pos_std[0]);
	p->pos[1] += rand_normal(0.0, pos_std[1]);
	p->pos[2] += rand_normal(0.0, pos_std[2]);

	double rv[3] = { rand_normal(0.0, rot_std[0]),
	                 rand_normal(0.0, rot_std[1]),
	                 rand_normal(0.0, rot_std[2]) };
	double dq[4];
	quat_from_rotvec(rv, dq);
	double q[4];
	quat_mul(p->quat, dq, q);   /* body-frame perturbation */
	quat_normalize(q);
	p->quat[0] = q[0]; p->quat[1] = q[1]; p->quat[2] = q[2]; p->quat[3] = q[3];
}

void pose_transform_point(const Pose *p, const double v[3], double out[3])
{
	double rv[3];
	quat_rotate(p->quat, v, rv);
	out[0] = rv[0] + p->pos[0];
	out[1] = rv[1] + p->pos[1];
	out[2] = rv[2] + p->pos[2];
}
