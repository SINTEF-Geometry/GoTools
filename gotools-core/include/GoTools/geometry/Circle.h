//==========================================================================
//                                                                          
// File: Circle.h                                                            
//                                                                          
// Created: Thu Oct 16 12:56:57 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Circle.h,v 1.5 2009-01-28 08:00:59 vsk Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _CIRCLE_H
#define _CIRCLE_H


#include "GoTools/geometry/ElementaryCurve.h"


namespace Go {


/// \brief Class that represents a circle. It is a subclass of
/// ElementaryCurve and thus has a parametrization.
///
/// A Circle has a natural parametrization in terms of the angle \a t:
/// \b p(\a t) = \b C + \b r(cos(\a t) \b x + sin(\a t) \b y), where
/// \b C is the centre, \b r is the radius, and \b x and \b y are the
/// (local) axes. The parametrization is (in principle) bounded:
/// \f$0 \leq t \leq 2\pi\f$. The dimension is either 2 or 3.

class Circle : public ElementaryCurve
{
public:
    /// Default constructor. Constructs an uninitialized Circle which
    /// can only be assigned to or read into.
    Circle()
    {};

    /// Constructor. Input is the radius, the centre, the normal to the
    /// plane of the circle, and the (approximate) direction of the
    /// x-axis. The dimension must be either 2 or 3.
    Circle(double radius,
	   Point centre, Point normal, Point x_axis);

    /// virtual destructor - ensures safe inheritance
    virtual ~Circle();

    /// Read object from stream
    /// \param is stream from which object is read
    virtual void read (std::istream& is);
    /// Write object to stream
    /// \param os stream to which object is written
    virtual void write (std::ostream& os) const;


    // --- Functions inherited from GeomObject ---

    virtual BoundingBox boundingBox() const;
    
    virtual int dimension() const;

    virtual ClassType instanceType() const;

    static ClassType classType();

    virtual Circle* clone() const;

    // --- Functions inherited from ParamCurve ---

    virtual void point(Point& pt, double tpar) const;

    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    virtual double startparam() const;
    virtual double endparam() const;

    virtual void reverseParameterDirection(bool switchparam = false);
    
    /// Circle parametrized on [0, 2*M_PI). Allowing class to be
    /// defined on a section of the circle.
    virtual void setParameterInterval(double t1, double t2);

    virtual SplineCurve* geometryCurve();

    /// Creates a SplineCurve representation of the circle. The spline
    /// representation is geometrically exact, but the parametrization
    /// is only an approximation. If setParameterBounds() have been
    /// used to set bounds on the parametrization, then the endpoints
    /// of the resulting spline curve \em does have correct parameter
    /// values.
    virtual SplineCurve* createSplineCurve() const;


    virtual bool isDegenerate(double degenerate_epsilon);

    virtual 
    Circle* subCurve(double from_par, double to_par,
		   double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    virtual DirectionCone directionCone() const;
 
    virtual void appendCurve(ParamCurve* cv, bool reparam=true);

    virtual void appendCurve(ParamCurve* cv,
			     int continuity, double& dist, bool reparam=true);

    virtual void closestPoint(const Point& pt,
			      double tmin,
			      double tmax,
			      double& clo_t,
			      Point& clo_pt,
			      double& clo_dist,
			      double const *seed = 0) const;

    virtual double length(double tol);

    // --- Functions specific to Circle ---

    /// Set bounds for the parametrization of the Circle other than
    /// the default \f$[0, 2\pi]\f$. Requirements for a valid
    /// parametrization are: 1) The first parameter must be strictly
    /// less than the second, 2) Parameter values must be in
    /// \f$[-2\pi, 2\pi]\f$, and 3) (\e endpar - \e startpar) must not
    /// exceed \f$2\pi\f$.
    /// \param startpar start parameter
    /// \param endpar end parameter
    void setParamBounds(double startpar, double endpar);


protected:

    double radius_;
    Point centre_;
    Point normal_;
    Point vec1_;
    Point vec2_;

    double startparam_;
    double endparam_;

    void setSpanningVectors();
};


} // namespace Go


#endif // _CIRCLE_H

