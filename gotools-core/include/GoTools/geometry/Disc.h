//===========================================================================
//
// File : Disc.h
//
// Created: Thu Nov  5 10:42:01 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#ifndef __DISC_H
#define __DISC_H


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Circle.h"


namespace Go
{


  class SplineSurface;


  /// \brief Class that represents a circular disc. It is a subclass of
  /// ElementarySurface, and has a natural parametrization by polar
  /// coordinates in terms of a radius \a r and angle \a v:
  /// \b p(\a r, \a v) = \a C + (cos \v) \b x + (\sin \v) \b y,
  /// where \b C is the centre position vector and \b x and \b y are
  /// the (local) axes. The parametrization is bounded by:
  /// \f$0 \leq r \leq R\f$ and \f$0 \leq 2\pi\f$, where \a R is the
  /// disc radius. The dimension is 2 or 3.
  ///
  /// A disc also holds degeneracy information for representing it as
  /// a SplineSurface object. There are two ways to parametrize.
  /// One is by polar coordinates, then parametrization of the
  /// SplineSurface object coincides with the one for the Disc
  /// object itself, i.e. self.point() gives (almoast) the same as
  /// geometrySurface()->point(). The only degeneracy point then is
  /// the centre. The other way is by splitting the boundary into
  /// four curves that become the boundary curves of the spline
  /// surface. Then the four meeting points of the curves
  /// become degeneracy points. This will be the default
  /// SplineSurface object representation.

  class Disc : public ElementarySurface
  {

  public:

    /// Default constructor. Constructs an uninitialized Disc which
    /// can only be assigned to or read into.
    Disc() { }

    /// Constructor. Input is the disc centre, disc radius, the x_axis
    /// and the normal (possibly approximate, only dummy in
    /// two-dimensional space). The local coordinate axes are
    /// normalized even if \c x_axis and/or \c normal are not unit
    /// vectors.
    Disc(Point centre, double radius, Point x_axis, Point normal);

    /// Virtual destructor - ensures safe inheritance
    virtual ~Disc() { }


    // Inherited from GeomObject
    /// read object from stream
    /// \param is stream from which object is read
    virtual void read (std::istream& is);

    // Inherited from GeomObject
    /// write object to stream
    /// \param os stream to which object is written
    virtual void write (std::ostream& os) const;

    // Inherited from GeomObject
    virtual int dimension() const;

    // Inherited from GeomObject
    virtual ClassType instanceType() const;

    // Inherited from GeomObject
    static ClassType classType()
    { return Class_Disc; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Disc* clone() const;

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

    void swapParameterDirection();

    bool isDegenerate(bool& b, bool& r,
		      bool& t, bool& l, double tolerance) const;


    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    /// Query if parametrization is bounded. All four parameter bounds
    /// must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded() const;

    /// Return the part of the surface limited by the given parameter bounds
    Disc* subSurface(double from_upar, double from_vpar,
		     double to_upar, double to_vpar,
		     double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Create a SplineSurface representation of the Disc.
    virtual SplineSurface* geometrySurface() const;

    /// Create a SplineSurface representation of the Disc.
    virtual SplineSurface*  createSplineSurface() const;

    virtual void setParameterBounds(double from_upar, double from_vpar,
				    double to_upar, double to_vpar);

    /// The NURBS representation of the disc will have degenereracy in the centre
    void useCentreDegen()
    { centre_degen_ = true; }

    /// The NURBS representation of the disc will have degenererate corners
    void useCornerDegen()
    { centre_degen_ = false; }


  private:

    Point centre_;
    double radius_;
    Point
      x_axis_,
      y_axis_,
      z_axis_;
    bool centre_degen_;  // If true, geometrySurface() gives a SplineSurface with
                         // degenerecy in the centre
                         // If false, let boundary curves of the SplineSurface
                         // lie on the boundary
    double degen_angles_[4];  // The angle parameter value giving the four degeneracy
                              // points on the boundary (only when center_degen = false)

    RectDomain domain_;

    void setCoordinateAxes();
    void setParameterDomain(double from_upar, double from_vpar,
			    double to_upar, double to_vpar);
    void setDefaultDomain();
    Circle boundaryCircle() const;

  };    // Class Disc


} // namespace Go



#endif    // #ifndef __DISC_H
