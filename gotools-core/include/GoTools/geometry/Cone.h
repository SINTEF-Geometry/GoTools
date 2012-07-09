//==========================================================================
//                                                                          
// File: Cone.h                                                              
//                                                                          
// Created: Tue Nov 18 15:18:34 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Cone.h,v 1.4 2009-02-19 15:52:54 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _CONE_H
#define _CONE_H


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Line.h"


namespace Go
{


class SplineSurface;


/// \brief Class that represents a cone. It is a subclass of
/// ElementarySurface, and thus has a parametrization and is
/// non-selfintersecting.
///
/// A Cone has a natural parametrization in terms of the angle \a
/// u and distance \a v:
/// \b p(\a u, \a v)
/// = \b C + (\a R + \a v tan \f$\alpha\f$)((cos \a u) \b x + (sin \a u) \b y) + \a v \b z,
/// where \b C is a position vector, \a R is the radius, \b x, \b
/// y and \b z are the (local) axes, and \f$\alpha\f$ is the cone angle.
/// The parametrization is bounded by: \f$0 \leq u \leq 2\pi\f$,
/// \f$-\infty < v < \infty\f$.  The dimension is 3.

class Cone : public ElementarySurface
{
public:
    /// Default constructor. Constructs an uninitialized Cone which
    /// can only be assigned to or read into.
    Cone()
    {};

    /// Constructor. Input is the radius, the location, the direction of
    /// the z-axis and the (possibly approximate) x-axis. The local
    /// coordinate axes are normalized even if \c z_axis and/or \c
    /// x_axis are not unit vectors.
    Cone(double radius, Point location, Point z_axis, Point x_axis,
         double cone_angle, bool isSwapped = false);

    /// Virtual destructor - ensures safe inheritance
    virtual ~Cone();

    /// read object from stream
    /// \param is stream from which object is read
    virtual void read (std::istream& is);
    /// write object to stream
    /// \param os stream to which object is written
    virtual void write (std::ostream& os) const;

    // Inherited from GeomObject
    virtual int dimension() const;

    // Inherited from GeomObject
    virtual ClassType instanceType() const;

    // Inherited from GeomObject
    static ClassType classType()
    { return Class_Cone; }

    /// Return empty box if infinite cone
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Cone* clone() const;

    // --- Functions inherited from ParamSurface ---

    const RectDomain& parameterDomain() const;

    std::vector<CurveLoop> allBoundaryLoops(double degenerate_epsilon
                                            = DEFAULT_SPACE_EPSILON) const;

    DirectionCone normalCone() const;
    DirectionCone tangentCone(bool pardir_is_u) const;

    void point(Point& pt, double upar, double vpar) const;
    void point(std::vector<Point>& pts, 
               double upar, double vpar,
               int derivs,
               bool u_from_right = true,
               bool v_from_right = true,
               double resolution = 1.0e-12) const;

    void normal(Point& n, double upar, double vpar) const;

    std::vector<shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const;

    std::vector<shared_ptr<ParamSurface> >
    subSurfaces(double from_upar, double from_vpar,
                double to_upar, double to_vpar,
                double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    double nextSegmentVal(int dir, double par, bool forward, double tol) const;

    void closestPoint(const Point& pt,
                      double&        clo_u,
                      double&        clo_v, 
                      Point&       clo_pt,
                      double&        clo_dist,
                      double         epsilon,
                      const RectDomain* domain_of_interest = NULL,
                      double   *seed = 0) const;

    void closestBoundaryPoint(const Point& pt,
                              double&        clo_u,
                              double&        clo_v, 
                              Point&       clo_pt,
                              double&        clo_dist,
                              double epsilon,
                              const RectDomain* rd = NULL,
                              double *seed = 0) const;

    void getBoundaryInfo(Point& pt1, Point& pt2,
                         double epsilon, SplineCurve*& cv,
                         SplineCurve*& crosscv, double knot_tol = 1e-05) const;

    void turnOrientation();

    void reverseParameterDirection(bool direction_is_u);

    //void swapParameterDirection();

    bool isDegenerate(bool& b, bool& r,
                      bool& t, bool& l, double tolerance) const;


    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    // --- Functions specific to Cone ---

    /// Radius at the location
    double getRadius() const
    { return radius_; }
    
    /// Point on cone axis
    Point getLocation() const
    { return location_;	}

    /// Local coordinate axes. The z_axis corresponds to the cone axis
    void getCoordinateAxes(Point& x_axis, Point& y_axis, Point& z_axis) const
    {
        x_axis = x_axis_;
        y_axis = y_axis_;
        z_axis = z_axis_;
    }

    /// The opening angle of the cone
    double getConeAngle() const
    { return cone_angle_; }

    /// Limit the cone surface by limiting the parameter domain
    void setParameterBounds(double from_upar, double from_vpar,
                            double to_upar, double to_vpar);

    /// Set parameter bounds in the \a u direction. \a u is the "angular"
    /// direction.
    void setParamBoundsU(double from_upar, double to_upar);

    /// Set parameter bounds in the \a v direction. \a v is the "linear"
    /// direction.
    virtual void setParamBoundsV(double from_vpar, double to_vpar);

    /// Query if parametrization is bounded. Only the \a v parameter
    /// direction is queried. The \a u parameter is always bounded.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded() const;

    /// Check if the surface is closed in the first parameter direction
    bool isClosed() const;

    /// Return the part of the cone surface limited by the parameter bounds
    Cone* subSurface(double from_upar, double from_vpar,
                         double to_upar, double to_vpar,
                         double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Create a SplineSurface representation of the cone.
    virtual SplineSurface* geometrySurface() const;

    /// Create a SplineSurface representation of the cone.
    virtual SplineSurface*  createSplineSurface() const;

    shared_ptr<Line> getLine(double upar) const; 

    shared_ptr<ElementaryCurve> 
      getElementaryParamCurve(ElementaryCurve* space_crv, double tol) const;

    // Confirm that this surface is axis rotational
    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle);

protected:

    double radius_;
    Point location_;
    Point z_axis_;
    Point x_axis_;
    Point y_axis_;
    double cone_angle_;

    RectDomain domain_;
    mutable RectDomain orientedDomain_; // Takes isSwapped_ flag into account

    void setCoordinateAxes();

};

} // namespace Go


#endif // _CONE_H

