//===========================================================================
//                                                                           
// File: Parabola.h                                                          
//                                                                           
// Created: Thu Jul  2 16:04:53 2009                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARABOLA_H
#define _PARABOLA_H




#include "GoTools/geometry/ElementaryCurve.h"


namespace Go {


/// \brief Class that represents a parabola. It is a subclass of
/// ElementaryCurve and thus has a parametrization.
///
/// An parabola has a natural parametrization in terms of its values:
/// f(t) = C + F*(t^2 * x + 2 * t * y). The
/// parametrization is unbounded: -\f$\infty < t < \infty\f$.

class Parabola : public ElementaryCurve
{
public:
    /// Default constructor. Constructs an uninitialized Parabola which
    /// can only be assigned to or read into.
    Parabola()
    {};

    /// Constructor. Input is location, direction and focal dist. The
    /// default Parabola is unbounded in its parametrization. To bound
    /// it, use setParamBounds().
    Parabola(Point location, Point direction, Point normal,
	     double focal_dist);

    /// virtual destructor - ensures safe inheritance
    virtual ~Parabola();

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

    virtual Parabola* clone() const;

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
    Parabola* subCurve(double from_par, double to_par,
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

    // --- Functions specific to Parabola ---

    /// Set bounds for the parametrization of the Parabola.
    /// \param startpar start parameter
    /// \param endpar end parameter
    virtual void setParamBounds(double startpar, double endpar);

    /// Query if parametrization is bounded. Both upper and lower
    /// parameter bounds must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded();

    // Translate the curve along a given vector
    virtual void translateCurve(const Point& dir);

   /// In 3D, the spanning vectors vec1_, vec2_, and the vector
    /// normal_ defines a right-handed coordinate system. 
     void setSpanningVectors();

    /// Check if the parabola lies in a plane with a given normal
    virtual bool isInPlane(const Point& norm,
			   double eps, Point& pos) const;

protected:

    Point location_;
    Point vec1_; // The axes.
    Point vec2_;
    Point normal_;

    double f_; // Focal dist.

    double startparam_;
    double endparam_;

};


} // namespace Go



#endif // _PARABOLA_H

