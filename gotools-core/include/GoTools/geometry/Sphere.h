//==========================================================================
//                                                                          
// File: Sphere.h                                                            
//                                                                          
// Created: Tue Nov 18 15:17:42 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Sphere.h,v 1.4 2009-02-05 13:56:33 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _SPHERE_H
#define _SPHERE_H


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Circle.h"


namespace Go
{


class SplineSurface;


/// \brief Class that represents a sphere. It is a subclass of
/// ElementarySurface, and thus has a parametrization and is
/// non-selfintersecting.
///
/// A Sphere has a natural parametrization in terms of the angles \a
/// u and \a v:
/// \b p(\a u, \a v)
/// = \b C + \a R cos \a v((cos \a u) \b x + sin(\a u) \b y) + \a R (sin\a v) \b z,
/// where \b C is a position vector, \a R is the radius, and \b x, \b
/// y and \b z are the (local) axes. The parametrization is
/// bounded by: \f$0 \leq u \leq 2\pi\f$, \f$-\frac{\pi}{2} < v < \frac{\pi}{2}\f$.
/// The dimension is 3.

class Sphere : public ElementarySurface
{
public:
    /// Default constructor. Constructs an uninitialized Sphere which
    /// can only be assigned to or read into.
    Sphere()
    {};

    /// Constructor. Input is the radius, the location, the direction
    /// of the z-axis and the (possibly approximate) x-axis. The local
    /// coordinate axes are normalized even if \c z_axis and/or \c
    /// x_axis are not unit vectors.
    Sphere(double radius, Point location, Point z_axis, Point x_axis);

    /// Virtual destructor - ensures safe inheritance
    virtual ~Sphere();

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
    { return Class_Sphere; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Sphere* clone() const
    { return new Sphere(radius_, location_, z_axis_, x_axis_); }


    // --- Functions inherited from ParamSurface ---

    const RectDomain& parameterDomain() const;

    CurveLoop outerBoundaryLoop(double degenerate_epsilon
				= DEFAULT_SPACE_EPSILON) const;
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

    // --- Functions specific to Sphere ---

    double getRadius() const
    { return radius_; }
    
    Point getLocation() const
    { return location_;	}

    void getCoordinateAxes(Point& x_axis, Point& y_axis, Point& z_axis) const
    {
	x_axis = x_axis_;
	y_axis = y_axis_;
	z_axis = z_axis_;
    }

    void setParameterBounds(double from_upar, double from_vpar,
			    double to_upar, double to_vpar);

    Sphere* subSurface(double from_upar, double from_vpar,
		       double to_upar, double to_vpar,
		       double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    virtual SplineSurface* geometrySurface() const;

    /// Create a SplineSurface representation of the Sphere.
    virtual SplineSurface*  createSplineSurface() const;

    /// Get the circle along the latitude for a given v parameter.
    /// \param vpar v parameter
    /// \return A circle for the corresponding v parameter. If the v
    /// parameter is bounded, only a segment of a full circle is
    /// returned.
    shared_ptr<Circle> getLatitudinalCircle(double vpar) const;

    /// Get the circle along a longitude for a given u parameter.
    /// \param upar u parameter
    /// \return A circle for the corresponding u parameter. If the u
    /// parameter is bounded, only a segment of a full circle is
    /// returned.
    shared_ptr<Circle> getLongitudinalCircle(double upar) const;

protected:

    double radius_;
    Point location_;
    Point z_axis_;
    Point x_axis_;
    Point y_axis_;

    RectDomain domain_;

    void setCoordinateAxes();

};

} // namespace Go


#endif // _SPHERE_H

