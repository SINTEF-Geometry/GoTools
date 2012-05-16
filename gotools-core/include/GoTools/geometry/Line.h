//==========================================================================
//                                                                          
// File: Line.h                                                              
//                                                                          
// Created: Mon Sep  8 13:32:23 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Line.h,v 1.5 2009-01-28 08:01:00 vsk Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _LINE_H
#define _LINE_H


#include "GoTools/geometry/ElementaryCurve.h"


namespace Go {


/// \brief Class that represents a line. It is a subclass of
/// ElementaryCurve and thus has a parametrization.
///
/// A Line has a natural parametrization in terms of its location \b
/// C and direction vector \b V: p(t) = C + tV.  This
/// parametrization is unbounded: -\f$\infty < t < \infty\f$.

class Line : public ElementaryCurve
{
public:
    /// Default constructor. Constructs an uninitialized Line which
    /// can only be assigned to or read into.
    Line()
    {};

    /// Constructor. Input is point and direction. The default Line is
    /// unbounded in its parametrization. To bound it, use
    /// setParamBounds().
    Line(Point point, Point direction);

    /// Constructor. Bounded line
    Line(Point point, Point direction, double length);

    /// Constructor. Bounded line with parameterization
    Line(Point point1, Point point2, double par1, double par2);

    /// virtual destructor - ensures safe inheritance
    virtual ~Line();

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

    virtual Line* clone() const;

    // --- Functions inherited from ParamCurve ---

    virtual void point(Point& pt, double tpar) const;

    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    virtual double startparam() const;
    virtual double endparam() const;

    virtual void reverseParameterDirection(bool switchparam = false);
    
     /// Limit the curve by limiting the parameter interval
   virtual void setParameterInterval(double t1, double t2);

    virtual SplineCurve* geometryCurve();
    /// Fetch spline representation of curve
    virtual SplineCurve* createSplineCurve() const;

    virtual bool isDegenerate(double degenerate_epsilon);

    virtual 
    Line* subCurve(double from_par, double to_par,
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

    // --- Functions specific to Line ---

    /// Set bounds for the parametrization of the Line.
    /// \param startpar start parameter
    /// \param endpar end parameter
    virtual void setParamBounds(double startpar, double endpar);

    // Translate the curve along a given vector
    virtual void translateCurve(const Point& dir);

    /// Query if parametrization is bounded. Both upper and lower
    /// parameter bounds must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded() const;

protected:

    Point location_;
    Point dir_;

    double startparam_;
    double endparam_;

};


} // namespace Go


#endif // _LINE_H

