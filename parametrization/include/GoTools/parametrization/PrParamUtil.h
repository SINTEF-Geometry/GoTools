#ifndef PRPARAMUTIL_H
#define PRPARAMUTIL_H

#include "GoTools/utils/Array.h"
using Go::Vector2D;
using Go::Vector3D;

/** Find the three barycentric coordinates of the point (x,y) with respect
 * to the triangle with vertices (x0,y0), (x1,y1), (x2,y2).
 * The three vertices should not be collinear!
 */
void baryCoords(double x, double y, double x0, double y0, double x1, double y1,
                double x2, double y2, double& tau0, double& tau1, double& tau2);

/** Find the barycentric coordinates of the origin (0,0) with respect
 * to the triangle formed by three vectors (u0,v0),(u0,v1),(u0,v2).
 * This is bit more efficient than calling  baryCoords.  
 */
//M.F. Feb 97.
void baryCoords0(double& u0, double& v0, double& u1, double& v1,
                 double& u2, double& v2, double& tau0, double& tau1, double& tau2);


double det(const double& u1, const double& v1,
	   const double& u2, const double& v2);

/// Find the signed area of triangle [(x0,y0),(x1,y1),(x2,y2)],
/// where (x0,y0),(x1,y1),(x2,y2) are assumed to be ordered anticlockwise.
double area(const double& x0, const double& y0,
	    const double& x1, const double& y1,
	    const double& x2, const double& y2);

/// Find the positive area of a 3D triangle.
double area(const Vector3D& a,
	    const Vector3D& b,
	    const Vector3D& c);

/// Represent a vector v in polar coordinates r(cos(theta),sin(theta))
/// where 0 <= r < 2 pi.
void polarCoords(Vector2D v, double& r, double& theta);

/// Represent a vector (u,v) in polar coordinates r(cos(theta),sin(theta))
/// where 0 <= r < 2 pi.
void polarCoords(double u, double v, double& r, double& theta);

/** Return tangent of half the angle between vectors b-a and c-a
 * without using trig functions.
 * Use fact that tan(alpha/2) = (1-cos(alpha)) / sin(alpha).
 * and use scalar and dot products to get cos(alpha) and sin(alpha).
 */
//*  M.F. Apr. 2002.
double tanThetaOverTwo(const Vector3D& a,
		       const Vector3D& b,
		       const Vector3D& c);

double cotangent(const Vector3D& a, const Vector3D& b, const Vector3D& c);


#endif // PRPARAMUTIL_H
