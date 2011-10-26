//===========================================================================
//                                                                           
// File: Hyperbola.h                                                         
//                                                                           
// Created: Thu Jul  2 16:05:05 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _HYPERBOLA_H
#define _HYPERBOLA_H




#include "GoTools/geometry/ElementaryCurve.h"


namespace Go {


/// \brief Class that represents a hyperbola. It is a subclass of
/// ElementaryCurve and thus has a parametrization.
///
/// A hyperbola has a natural parametrization in terms of its values:
/// f(t) = C + (R_1 * cosh(t))*x + (R_2 * sinh(t)) * y. This
/// parametrization is ubounded: -\f$\infty < t < \infty\f$.

class Hyperbola : public ElementaryCurve
{
public:
    /// Default constructor. Constructs an uninitialized Hyperbola which
    /// can only be assigned to or read into.
    Hyperbola()
    {};

    /// Constructor. Input is location, direction and the two
    /// semi-axes. The default hyperbola is unbounded in its
    /// parametrization. To bound it, use setParamBounds().
    Hyperbola(Point location, Point direction, Point normal,
	      double r1, double r2);

    /// virtual destructor - ensures safe inheritance
    virtual ~Hyperbola();

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

    virtual Hyperbola* clone() const;

    // --- Functions inherited from ParamCurve ---

    virtual void point(Point& pt, double tpar) const;

    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    virtual double startparam() const;
    virtual double endparam() const;

    virtual void reverseParameterDirection(bool switchparam = false);
    
    virtual void setParameterInterval(double t1, double t2);

    virtual SplineCurve* geometryCurve();

    virtual bool isDegenerate(double degenerate_epsilon);

    virtual
    Hyperbola* subCurve(double from_par, double to_par,
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

    // --- Functions specific to Hyperbola ---

    /// Set bounds for the parametrization of the Hyperbola.
    /// \param startpar start parameter
    /// \param endpar end parameter
    void setParamBounds(double startpar, double endpar);

    /// Query if parametrization is bounded. Both upper and lower
    /// parameter bounds must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded();

    void setSpanningVectors();

protected:

    Point location_;
    Point vec1_; // The axis wrt r1_ (input direction). 
    Point vec2_; // The axis wrt r2_.
    Point normal_;

    double r1_; // semi_imag_axis;
    double r2_; // semi_axis;

    double startparam_;
    double endparam_;

};


} // namespace Go



#endif // _HYPERBOLA_H

