#ifndef _CURVATUREUTILS_
#define _CURVATUREUTILS_

#include "GoTools/utils/Point.h"
#include <vector>

namespace Go
{
    /** Help functions in related to curvature.
     */

    /// Given position, first and second derivative
    /// of a curve passing through a point, compute
    /// the unit tangent, curvature vector and curvature 
    /// radius of this curve.
    double curvatureRadius(const std::vector<Point>& der,
			   std::vector<Point>& unitder);

    /// Computes the step length along a curve based on radius of curvature
    /// at a point on the curve, and an absolute tolerance.
    double stepLenFromRadius(double radius, double aepsge);

    /// To create the tangent length for interpolating a
    /// circular arc with an almost equi-oscillating Hermit qubic.
    double tanLenFromRadius(double radius, double angle);

    /// Given position, first and second derivative in both ends of
    /// an Hermite segment, compute parameter interval and tangent lengths
    /// in order to stay close to a circular segment.
    void getHermiteData(const std::vector<Point>& der1,
			const std::vector<Point>& der2, 
			double& parint, double& len1, double& len2);

} // End of namespace Go


#endif

