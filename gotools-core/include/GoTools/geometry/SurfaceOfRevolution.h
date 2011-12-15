//==========================================================================
//                                                                          
// File: SurfaceOfRevolution.h                                               
//                                                                          
// Created: Wed Oct 21 12:51:21 2009                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id$
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _SURFACEOFREVOLUTION_H
#define _SURFACEOFREVOLUTION_H


#include "GoTools/geometry/ParamSurface.h"



namespace Go {


class SplineCurve;
class SplineSurface;


/// \brief Class that represents a surface of revolution. A
/// SurfaceOfRevolution is swept out by a SplineCurve that is rotated
/// around an axis with a complete revolution, and is thereby a
/// parametric surface. 
///
/// The parametrization of the surface is given in terms of a location
/// \b C, an axis line \b V, and the spline curve \b \f$\lambda(v)\f$
/// with parameter \a v:
///
/// \f$ \sigma(u, v)
/// = C + (\lambda(v) - C) \cos u
/// + ((\lambda(v) - C) \cdot V) V (1 - \cos u)
/// + V \times (\lambda(v) - C) \sin u
/// \f$
///
/// The parameter \a u is bounded by: \f$0 \leq u \leq 2\pi\f$. The
/// axis \b V is normalized.
/// 
/// The curve \f$\lambda\f$ must be such that it doesn't lead to a
/// self-intersecting surface.


class SurfaceOfRevolution : public ParamSurface
{
public:
    /// Default constructor. Constructs an uninitialized
    /// SurfaceOfRevolution which can only be assigned to or read
    /// into.
    SurfaceOfRevolution()
    {};

    /// Constructor. Input is the location and normalized direction of
    /// the axis, and the SplineCurve that is swept out by the
    /// revolution.
    SurfaceOfRevolution(Point location, Point axis_dir,
			std::shared_ptr<SplineCurve> curve);

    /// Virtual destructor - ensures safe inheritance
    virtual ~SurfaceOfRevolution();

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
    static ClassType classType();

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual SurfaceOfRevolution* clone() const;


    // --- Functions inherited from ParamSurface ---

    const RectDomain& parameterDomain() const;
    virtual RectDomain containingDomain() const;

    virtual bool inDomain(double u, double v) const;

    virtual Point closestInDomain(double u, double v) const;

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

    std::vector<std::shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const;

    std::vector<std::shared_ptr<ParamSurface> >
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

    virtual double area(double tol) const;

    bool isDegenerate(bool& b, bool& r,
		      bool& t, bool& l, double tolerance) const;


    /// Check for parallel and anti parallel partial derivatives in
    /// surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, 
				      double tol) const;


    /// Return surface corners, geometric and parametric points in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    // --- Functions specific to SurfaceOfRevolution ---

    Point getLocation() const
    { return location_;	}

    Point getAxisDir() const
    { return axis_dir_;	}

    std::shared_ptr<SplineCurve> getCurve() const
    { return curve_; }

    void setParameterBounds(double from_upar, double from_vpar,
			    double to_upar, double to_vpar);

    SurfaceOfRevolution* subSurface(double from_upar, double from_vpar,
				    double to_upar, double to_vpar,
				    double fuzzy
				    = DEFAULT_PARAMETER_EPSILON) const;

    // Is "geometrySurface()" a good name for this function? @jbt
    virtual SplineSurface* geometrySurface() const;

private:
    Point location_;
    Point axis_dir_;
    std::shared_ptr<SplineCurve> curve_;

    RectDomain domain_;
    void setDefaultDomain();

};


} // namespace Go


#endif // _SURFACEOFREVOLUTION_H
