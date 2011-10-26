//===========================================================================
//                                                                           
// File: SplineApproximator.h                                              
//                                                                           
// Created: Wed Oct 18 14:00:55 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SplineApproximator.h,v 1.6 2007-12-04 16:12:01 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SPLINEAPPROXIMATOR_H
#define _SPLINEAPPROXIMATOR_H


#include "GoTools/geometry/Interpolator.h"
#include "GoTools/geometry/BsplineBasis.h"

namespace Go
{


    /** An Interpolator that generates a spline curve approximating
     *  the given dataset in the least squares sense.
     */
class SplineApproximator : public Interpolator
{
public:
    /// Constructor takes no arguments
    SplineApproximator()
	: num_coefs_(0)
    {
    }

    /// Virtual destructor ensures safe inheritance
    virtual ~SplineApproximator();

    // inherited from Interpolator
    virtual const BsplineBasis& basis();

    /// The interpolating function, as inherited by \ref Interpolator.  
    /// Prior to calling this function, the user must have specified:
    /// - \em either the number of control points to use in the 
    ///   approximating curve by \ref setNumCoefs().
    /// - \em or/and directly specified the BsplineBasis to use
    ///   by \ref setSplineSpace().
    /// A default BsplineBasis will be generated if the user has only set
    /// the number of control points prior to calling interpolate().
    /// For parameter list, see \ref Interpolator.
    virtual void interpolate(int num_points,
			     int dimension,
			     const double* param_start,
			     const double* data_start,
			     std::vector<double>& coefs);
    
    /// Specify the number of basis functions / control points to use
    /// in the approximating curve.
    void setNumCoefs(int num) {
	num_coefs_ = num;
    }

    /// Directly specify the spline space in which to search for the 
    /// approximating function.
    void setSplineSpace(const BsplineBasis& basis) {
	basis_ = basis;
    }


private:
    int num_coefs_;
    BsplineBasis basis_;
};


} // namespace Go


#endif // _SPLINEAPPROXIMATOR_H

