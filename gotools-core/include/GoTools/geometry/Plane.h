//===========================================================================
//                                                                           
// File: Plane.h
//                                                                           
// Created: April 2008
//                                                                           
// Author: Vibeke Skytt, SINTEF
//                                                                           
// Revision: $Id: Plane.h,v 1.11 2009-02-10 13:17:08 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PLANE_H
#define _PLANE_H


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/utils/RotatedBox.h"


namespace Go
{


class SplineSurface;


/// \brief Class that represents a plane. It is a subclass of
/// ElementarySurface, and thus has a parametrization and is
/// non-selfintersecting.
///
/// A Plane has a natural parametrization in terms of its location \b
/// C and spanning vectors \b x and \b y: p(u, v) = C + ux + vy.  This
/// parametrization might be unbounded: -\f$\infty < u,v < \infty\f$.

class Plane : public ElementarySurface
{
public:
    /// Default constructor. Constructs an uninitialized Plane which
    /// can only be assigned to or read into.
    Plane()
    {};

    /// Constructor. Input is location and normal
    Plane(Point location, Point normal);

    /// Constructor. Input is location, normal and (local,
    /// approximate) x-axis
    Plane(Point location, Point normal, Point x_axis);

    /// Constructor. Input is coefficients of implicit equation +
    /// point suggestion
    Plane(double a, double b, double c, double d);

    /// virtual destructor - ensures safe inheritance
    virtual ~Plane();

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
    { return Class_Plane; }

    /// Return empty box if infinite plane
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Plane* clone() const
    { return new Plane(location_, normal_); }


    // --- Functions inherited from ParamSurface ---

    const Domain& parameterDomain() const;

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

    void swapParameterDirection();

    bool isDegenerate(bool& b, bool& r,
                      bool& t, bool& l, double tolerance) const;


    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    // --- Functions specific to Plane ---

    /// Point in plane
    Point getPoint()
    { return location_; }

    /// Plane normal
    Point getNormal()
    { return normal_; }
    
    /// Vectors in plane
    void getSpanningVectors(Point& axis1, Point& axis2)
    {
        axis1 = vec1_;
        axis2 = vec2_;
    }

    /// Projection of the point pnt in the plane
    Point projectPoint(const Point& pnt) const;

    /// Distance between the point pnt and the plane
    double distance(const Point& pnt) const;

    /// Restrict the plane by restricting the parameter domain. It is
    /// initially infinite
    void setParameterBounds(double from_upar, double from_vpar,
                            double to_upar, double to_vpar);

    /// Fetch a part of the plane
    Plane* subSurface(double from_upar, double from_vpar,
                      double to_upar, double to_vpar,
                      double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Create a SplineSurface representation of the Plane.
    virtual SplineSurface* geometrySurface() const;

    /// Create a SplineSurface representation of the Plane.
    virtual SplineSurface*  createSplineSurface() const;

    /// Query if parametrization is bounded. All four parameter bounds
    /// must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded() const;

    /// Return the result from intersecting the unbounded plane with a
    /// rotated bounding box (having axis[0]=vec1_, axis[1]=vec2_,
    /// axis[2]=normal_). Useful for visualizing the (unbounded) plane.
    /// If intersection is empty, the returned plane is the NULL
    /// pointer.  The rotated box may be created from a boundingbox by
    /// defining the coordinate system and the 8 corner points of the
    /// bd_box.
    Plane* intersect(const RotatedBox& bd_box) const;

protected:

    Point location_;

    // The vectors vec1_, vec2_, and normal_ define a right-handed
    // coordinate system.
    Point normal_;
    Point vec1_;
    Point vec2_;

    RectDomain domain_;

    void setSpanningVectors();

};

} // namespace Go



#endif // _PLANE_H

