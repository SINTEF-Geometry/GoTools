//===========================================================================
//                                                                           
// File: HermiteInterpolator.h                                              
//                                                                           
// Created: 01.08.30
//                                                                           
// Author: Vibeke Skytt, SINTEF
//                                                                           
// Revision: 
//                                                                           
// Description: Does hermite spline interpolation.
//                                                                           
//===========================================================================

#ifndef _HERMITEINTERPOLATOR_H
#define _HERMITEINTERPOLATOR_H


#include "GoTools/geometry/Interpolator.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/utils/Point.h"

namespace Go
{

/** An Interpolator that generates a hermite spline
 *  curve through given points and tangents.
 */

class GO_API HermiteInterpolator : public Interpolator
{
 public:
    /// Constructor takes no arguments
    HermiteInterpolator() {}
 
    /// Virtual destructor enables safe inheritance
    virtual ~HermiteInterpolator();

    // inherited from Interpolator
    virtual const BsplineBasis& basis();

    /// Hermite interpolation of a sequence of points with
    /// associated tangents and parameter values
    /// \param num_points number of points to interpolate
    /// \param dimension dimension of points to interpolate (2D, 3D, etc..)
    /// \param param_start pointer to the start of the array where
    ///                    the parameter values of the points are stored.
    ///                    This should be a strictly increasing sequence
    ///                    of 'num_points' values.
    /// \param data_start pointer to the start of the array where the
    ///                   points and tangents to be interpolated are stored.
    ///                   Each point and tangent consist of 'dimension'
    ///                   coordinates, and each tangent is stored immediately
    ///                   after its corresponding point.
    /// \retval coefs The control points of the computed hermite interpolation
    ///               curve will be returned in this vector. (Use the basis() 
    ///               function to get the associated b-spline basis).
    virtual void interpolate(int num_points,
			     int dimension,
			     const double* param_start,
			     const double* data_start,
			     std::vector<double>& coefs);

    /// Hermite interpolation of a sequence of points with 
    /// associated tangents and parameter values.
    /// \param data This vector contains the set of points and 
    ///             tangents to be interpolated.  The size of the
    ///             vector is twice the total number of points.
    ///             Each entry on the form [2i] represents a point,
    ///             and the entry [2i+1] represents the associated
    ///             tangent.
    /// \param param This vector represent the parameterization of the
    ///              points, and should have one entry per point.
    ///              (Making it half the size of the 'data' vector).
    ///              The parameter sequence must be strictly increaing.
    /// \retval coefs The control points of the computed hermite 
    ///               interpolation curve will be returned in this 
    ///               vector.  (Use the basis() function to get the 
    ///               associated b-spline basis.
    void interpolate(const std::vector<Point>& data,
		     const std::vector<double>& param,
		     std::vector<double>& coefs);

 private:
    BsplineBasis basis_;
}; 

} // namespace Go

#endif // _HERMITEINTERPOLATOR_H


 
