/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#ifndef _PARAMCURVE_H
#define _PARAMCURVE_H

#include "GoTools/utils/config.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/utils/CompositeBox.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/utils/Values.h"
#include <vector>
#include <memory>


namespace Go
{

class SplineCurve;

    /** Base class for parametric curves in Go
     *
     */

class GO_API ParamCurve : public GeomObject
{
public:
    /// virtual destructor - ensures safe inheritance
    virtual ~ParamCurve();

    /// Evaluate the curve's position at a given parameter
    /// \param pt the evaluated position will be written to this Point
    /// \param tpar the parameter for which we wish to evaluate the curve
    virtual void point(Point& pt, double tpar) const = 0;

    /// Evaluate the curve's position and a certain number of derivatives
    /// at a given parameter.
    /// \param pts the evaluated position and derivatives (tangent, curvature vector, etc.)
    ///            will be written to this vector.  The first entry will be the position,
    ///            the second entry will be the first derivative, etc.  The size of this
    ///            vector must be set to 'derivs'+ 1 prior to calling this function.
    /// \param tpar the parameter for which we want to evaluate the curve
    /// \param derivs the number of derivatives we want to have calculated
    /// \param from_right specify whether we should calculate derivatives 'from the right'
    ///                   or 'from the left' (default is from the right).  This matters 
    ///                   only when the curve presents discontinuities in its derivatives.
    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const = 0;


    /// Evaluate the curve's position at a certain parameter
    /// \param tpar the parameter for which we want to evaluate the curve's position
    /// \return the curve's position for the parameter 'tpar'.
    /// NB: This function is implemented in terms of the ParamCurve's virtual
    /// 'point(...)' functions, but is itself not virtual.  If you make a concrete 
    /// subclass and wish to make this function visible to the user, you must
    /// put a 'using ParamCurve::point' statement in the class definition.
    Point point(double tpar) const;

    /// Evaluate the curve's position and a certain number of derivatives 
    /// at a given parameter.
    /// \param tpar the parameter for which we want to evaluate the curve
    /// \param derivs the number of derivatives we want to have calculated
    /// \param from_right specify whether we should calculate derivatives 'from the right'
    ///                   or 'from the left' (default is from the right).  This matters 
    ///                   only when the curve presents discontinuities in its derivatives.
    /// \return a STL vector containing the evaluated position and the specified number
    ///         of derivatives.  The first entry will be the position, the second entry
    ///         will be the first derivative, etc.  The total size of the returned vector
    ///         will be 'derivs' + 1.
    /// NB: This function is implemented in terms of the ParamCurve's virtual
    /// 'point(...)' functions, but is itself not virtual.  If you make a concrete 
    /// subclass and wish to make this function visible to the user, you must
    /// put a 'using ParamCurve::point' in the class definition.
    std::vector<Point> point(double tpar, int derivs, bool from_right = true) const;

    /// Evaluate points on a regular set of parameter values
    /// \param num number of values to evaluate
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param param upon function return, this vector holds all the numerical values
    ///              where evaluation has taken place
    virtual void uniformEvaluator(int num, std::vector<Point>& points, std::vector<double>& param) const;

    /// Query the start parameter of the curve
    /// \return the curve's start parameter
    virtual double startparam() const = 0;

    /// Query the end parameter of the curve
    /// \return the curve's end parameter
    virtual double endparam() const = 0;

    /// Set the parameter direction of the curve.  The curve's parameter interval will 
    /// always remain constant, but by flipping the parameter direction, the curve will
    /// be traced the opposite way when moving a parameter over the parameter interval.
    /// \param switchparam if true, and the curve is 2D, the x and y
    /// coordinates should be swapped. This is used when turning the
    /// orientation of bounded (trimmed) surfaces.
    virtual void reverseParameterDirection(bool switchparam = false) = 0;
    
    /// Linear reparametrization. The meaning is changed for elementary curves
    virtual void setParameterInterval(double t1, double t2) = 0;

    /// If the definition of this ParamCurve contains a SplineCurve describing its 
    /// spatial shape, then this function will return a pointer to this SplineCurve.
    /// Otherwise it will return a null pointer.
    /// The returned curve is NEWed, so the user is responsible
    /// for deleting it.
    /// This function may have side-effects.
    /// \return a pointer to a SplineCurve representation of the ParamCurve,
    /// if it exists.  Null pointer otherwise.
    virtual SplineCurve* geometryCurve() = 0;

    /// Return the spline curve associated to this curve, if any
    /// No copying
    virtual SplineCurve* getSplineCurve() 
    {
      return 0;  // Default behaviour
    }

    /// Query whether the curve is degenerate (collapsed into a single point).
    /// \param degenerate_epsilon the tolerance used in determine whether the curve
    ///        is degenerate.  A curve is considered degenerate if its total length 
    ///        is shorter than this value.
    /// \return \c true if the curve is degenerate, \c false otherwise.
    virtual bool isDegenerate(double degenerate_epsilon) = 0;

    /// Query whether the curve is closed. Periodic curves like
    /// circles and ellipses are closed.
    /// \return \c true if the curve is closed, \c false otherwise.
    virtual bool isClosed();

    /// Returns a curve which is a part of this curve.  The result is
    /// NEWed, so the user is responsible for deleting it.  NB: It is
    /// not guaranteed that the ParamCurve that is returned is of the
    /// same type as the curve itself. Thus, the returned curve might
    /// be a SplineCurve.
    /// \param from_par start value of parameter interval that will
    /// define the subcurve
    /// \param to_par end value of parameter interval that will define
    /// the subcurve
    /// \param fuzzy since subCurve works on those curves who are
    /// spline-based, this tolerance defines how close the start and
    /// end parameter must be to an existing knot in order to be
    /// considered \em on the knot.
    /// \return a pointer to a new subcurve which represents the part
    /// of the curve between 'from_par' and 'to_par'.  It will be
    /// spline-based and have a k-regular knotvector.  The user is
    /// responsible for deleting this subcurve when it is no longer
    /// needed.
    virtual ParamCurve* subCurve(double from_par, double to_par,
				 double fuzzy = DEFAULT_PARAMETER_EPSILON) const = 0;

    /// Split curve in a specified parameter value
    virtual
    std::vector<shared_ptr<ParamCurve> > 
      split(double param,
	    double fuzzy = DEFAULT_PARAMETER_EPSILON) const; 

    /// The clone-function is herited from GeomObject, but overridden here to get a 
    /// covariant return type (for those compilers that allow this).
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
//     virtual GeomObject* clone() const = 0;
// #else
//     virtual ParamCurve* clone() const = 0;
// #endif //_MSC_VER < 1300
// #else
    virtual ParamCurve* clone() const = 0;
// #endif

    /// Creates a DirectionCone which covers all tangent directions of this
    /// curve.
    /// \return the smallest DirectionCone containing all tangent directions of 
    ///         this curve.
    virtual DirectionCone directionCone() const = 0;
 
    /// Creates a composite box enclosing the curve. The composite box
    /// consists of an inner and an edge box. The inner box is
    /// supposed to be made from the interior of the curve, while the
    /// edge box is made from the start and end points. The default
    /// implementation simply makes both boxes identical to the
    /// regular bounding box.
    /// \return the CompositeBox enclosing the curve.
    virtual CompositeBox compositeBox() const;

    /// append a curve to this curve, with eventual reparametrization
    /// NB: This virtual member function currently only works for SplineCurves 
    /// and CurveOnSurfaces.  Moreover, 'this' curve and the 'cv' curve must 
    /// be of the same type.
    /// \param cv the curve to append to 'this' curve.
    /// \param reparam specify whether or not there should be reparametrization
    virtual void appendCurve(ParamCurve* cv, bool reparam=true) = 0;

    /// append a curve to this curve, with eventual reparmetrization
    /// \param cv the curve to append to 'this' curve.
    /// \param continuity the required continuity at the transition.  Can be G^(-1) 
    ///        and upwards.
    /// \param dist a measure of the local distorsion around the transition in order
    ///             to achieve the specified continuity.
    /// \param reparam specify whether or not there should be reparametrization
    /// \param tol relevant only for CurveOnSurface
    virtual void appendCurve(ParamCurve* cv,
			     int continuity, double& dist, bool reparam=true,
			     double tol = 1.0e-4) = 0;

    /// Estimate the length of the curve, by sampling it at a certain number of points
    /// and calculating the linear approximation to the curve through these points.
    /// \param numpts number of sample points used
    /// \return the estimated curve length
    double estimatedCurveLength(int numpts = 4) const;

    /// Estimate the length of an interval of the curve, by sampling it at a certain 
    /// number of points in the interval and calculating the linear approximation 
    /// through these points.
    /// \param tmin parameter at start of interval
    /// \param tmax parameter at end of interval
    /// \param numpts number of sample points used
    /// \return the estimated curve length
    double estimatedCurveLength(double tmin, double tmax, int numpts = 4) const;

    /// Compute the closest point from an interval of this curve to a specified point.
    /// \param pt point we want to find the closest point to
    /// \param tmin start parameter of search interval
    /// \param tmax end parameter of search interval
    /// \param clo_t upon function return, 'clo_t' will contain the parameter value of 
    ///              the closest point found.
    /// \param clo_pt upon function return, 'clo_pt' will contain the position of the
    ///               closest point found.
    /// \param clo_dist upon function return, 'clo_dist' will containn the distance 
    ///                 between 'pt' and the closest point found.
    /// \param seed pointer to initial guess value, provided by the user (can be 0,
    ///             for which the algorithm will determine a (hopefully) reasonable
    ///             choice).
    virtual void closestPoint(const Point& pt,
			      double tmin,
			      double tmax,
			      double& clo_t,
			      Point& clo_pt,
			      double& clo_dist,
			      double const *seed = 0) const = 0;

    /// Compute the closest point from this curve to a specified point, taking the whole
    /// curve into account (not just an interval of it).
    /// \param pt point we want to find the closest point to
    /// \param clo_t upon function return, 'clo_t' will contain the parameter value of 
    ///              the closest point found.
    /// \param clo_pt upon function return, 'clo_pt' will contain the position of the
    ///               closest point found.
    /// \param clo_dist upon function return, 'clo_dist' will containn the distance 
    ///                 between 'pt' and the closest point found.
    void closestPoint(const Point& pt, double& clo_t, Point& clo_pt, double& clo_dist) const;

    /// If the ParamCurve is divided up into logical segments, this function will return 
    /// the parameter value of the "next segment", starting from a parameter given by the user.
    /// If no division into logical segments exist, then it is the start- or end parameter that
    /// is returned.
    /// \param par the parameter from which we start the search for the next segment.
    /// \param forward whether to search forwards or backwards along the parameter domain.
    /// \param tol the tolerance to determine or not 'par' is already located on the start of 
    ///            the next segment.
    /// \return the parameter value of the next segment.
    virtual double nextSegmentVal(double par, bool forward, double tol) const;


    /// Compute the total length of this curve up to some tolerance
    /// \param tol the relative tolerance when approximating the length, i.e. stop iteration
    ///            when error becomes smaller than tol/(curve length)
    /// \return the length calculated
    virtual double length(double tol) = 0;

    /// Compute the length of a segment of this curve up to some tolerance
    /// \param tol the relative tolerance when approximating the length, i.e. stop iteration
    ///            when error becomes smaller than tol/(curve length)
    /// \param tstart the parameter value for the start point of the segment
    /// \param tend the parameter value for the end point of the segment
    /// \return the length calculated
    virtual double length(double tol, double tstart, double tend);

    /// Check if the curve is axis rotational. Only true if a connection
    /// to an axis rotational elementary curve exist
    /// The axis and rotational angle is only specified if the curve
    /// is actually rotational
    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle)
    {
      return false;  // Default behaviour, overriden for spline curves
      // bounded curves and some elementary curves
    }

    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle, double& radius)
    {
      return false;  // Default behaviour, overriden for spline curves
      // bounded curves and some elementary curves
    }

    virtual bool isBounded() const
    {
      return true;  // Is overridden when an unbounded curve is possible
    }
    
    /// Check if the curve is linear
    virtual bool isLinear(Point& dir, double tol);

    /// Check if the curves lies in a plane passing through a given axis
    virtual bool isInPlane(const Point& loc, const Point& axis,
			   double eps, Point& normal) const;

    /// Check if the curve lies in a plane with a given normal
    virtual bool isInPlane(const Point& norm,
			   double eps, Point& pos) const;

protected:
    void closestPointGeneric(const Point&   pt,
			     double    tmin,
			     double    tmax,
			     double guess_param,
			     double&   clo_t,
			     Point&  clo_pt,
			     double&   clo_dist) const;

    void s1771(Point pt,double aepsge,
	       double astart,double aend,double anext,double &cpos,int *jstat) const;

    void s1771_s9point(Point pt, std::vector<Point> val, Point diff,
			  double astart,double aend,int max_it,double *cnext,double *ad,
		       double adel,double *cdist,double aprev,int *jstat) const;

    double s1771_s9del(double *eco,double *eco1,double *eco2,int idim) const;
};

} // namespace Go



#endif // _PARAMCURVE_H

