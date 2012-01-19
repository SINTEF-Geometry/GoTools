//===========================================================================
//                                                                           
// File: PointOnCurve
//                                                                           
// Created: Dec. 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _POINT_ON_CURVE_H
#define _POINT_ON_CURVE_H

#include "GoTools/utils/Point.h"
#include <memory>
#include "GoTools/geometry/ParamCurve.h"

namespace Go
{

/** Represents a point on a curve. Lean functionality for the time being, but should
    be extended to cover for instance curvature computations  */

class PointOnCurve
    {
    public:
      /// Default constructor
      PointOnCurve();  

	/// Constructor taking a curve and a parameter value
	PointOnCurve(shared_ptr<ParamCurve> curve, double par);

	/// Constructor taking a curve and a point
	PointOnCurve(shared_ptr<ParamCurve> curve, Point pnt);

	/// Destructor
	~PointOnCurve();

	/// Evaluate
	/// If an input point is given, this is the point that will be returned even if it
	/// does not lie exactly on the curve
	Point getPos() const;

	/// Fetch the curve associated to the point
	shared_ptr<ParamCurve> getCurve() const
	    {
		return crv_;
	    }

	/// Fetch the curve parameter associated to the point
	double getPar() const
	    {
		return par_;
	    }

	/// The curve is evaluated including der derivatives
	void evaluate(int der, std::vector<Point>& deriv) const;

	/// Add limiting information to the current curve
	void setParInterval(double start, double end);

    private:
	Point point_;  // The point
	double par_; // Parameter value on the curve corresponding to this point
	double t1_, t2_;  // End parameters of current parameter interval
	shared_ptr<ParamCurve> crv_;  // The curve

};


} // namespace Go

#endif // _POINT_ON_CURVE_H

