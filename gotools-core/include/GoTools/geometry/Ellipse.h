//===========================================================================
//                                                                           
// File: Ellipse.h                                                           
//                                                                           
// Created: Thu Jul  2 16:00:44 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _ELLIPSE_H
#define _ELLIPSE_H


#include "GoTools/geometry/ElementaryCurve.h"


namespace Go {


/// \brief Class that represents an ellipse. It is a subclass of
/// ElementaryCurve and thus has a parametrization.
///
/// An ellipse has a natural parametrization in terms of its values:
/// f(t) = C + (R_1 * cos(t))*x + (R_2 * sin(t)) * y. This
/// parametrization is bounded: 0 <= t <= 360 (degrees).

class Ellipse : public ElementaryCurve
{
public:
    /// Default constructor. Constructs an uninitialized Ellipse which
    /// can only be assigned to or read into.
    Ellipse()
    {};

    /// Constructor. Input is point, axis direction and lengths of the
    /// two semi-axis.
    Ellipse(Point centre, Point direction, Point normal,
	    double r1, double r2);

    /// virtual destructor - ensures safe inheritance
    virtual ~Ellipse();

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

    virtual Ellipse* clone() const;

    // --- Functions inherited from ParamCurve ---

    virtual void point(Point& pt, double tpar) const;

    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    virtual double startparam() const;
    virtual double endparam() const;

    virtual void reverseParameterDirection(bool switchparam = false);
    
    /// Limit the curve by limiting the parameter domain
    virtual void setParameterInterval(double t1, double t2);

    virtual SplineCurve* geometryCurve();
    /// Create the spline representation of this curve
    virtual SplineCurve* createSplineCurve() const;

    virtual bool isDegenerate(double degenerate_epsilon);

    virtual 
    Ellipse* subCurve(double from_par, double to_par,
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

    // --- Functions specific to Ellipse ---

    /// Set bounds for the parametrization of the Ellipse.
    /// \param startpar start parameter
    /// \param endpar end parameter
    void setParamBounds(double startpar, double endpar);


protected:

    Point centre_; // Center of the ellipse.
    Point vec1_; // The axis wrt r1_.
    Point vec2_; // The axis wrt r2_.
    Point normal_;

    double r1_; // semi_axis_1;
    double r2_; // semi_axis_2;

    double startparam_; // At least 0.0.
    double endparam_; // At most 2.0*M_PI.

    void setSpanningVectors();

};


} // namespace Go



#endif // _ELLIPSE_H

