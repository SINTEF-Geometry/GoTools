//===========================================================================
//
// File : SphereVolume.h
//
// Created: Mon Nov  9 11:05:00 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#ifndef __SPHEREVOLUME_H
#define __SPHEREVOLUME_H



#include "GoTools/trivariate/ElementaryVolume.h"


namespace Go
{


  class SplineVolume;


  /// \brief Class that represents a solid sphere. It is a subclass of
  /// ElementaryVolume, and thus has a parametrization.
  ///
  /// A SphereVolume has a natural parametrization by spherical coordinates in terms
  /// the distance \a u from the centre, the elevation angle \a v and the azimuth
  /// angle \a w:
  /// \b p(\a u, \a v, \a w)
  /// = \b C + \a u cos \a w((cos \a v) \b x + sin(\a v) \b y) + \a u (sin\a w) \b z,
  /// where \b C is a position vector, and \b x, \b y and \b z are the (local) axes.
  /// The parametrization is bounded by: \f$0 \leq u \leq R\f$, \f$0 \leq v \leq 2\pi\f$,
  /// \f$-\frac{\pi}{2} < w < \frac{\pi}{2}\f$, where \a R is the radius. The dimension is 3.

  class SphereVolume : public ElementaryVolume
  {
  public:
    /// Default constructor. Constructs an uninitialized SphereVolume which
    /// can only be assigned to or read into.
    SphereVolume() { }

    /// Constructor. Input is the radius, the location, the direction of
    /// the z-axis and the (possibly approximate) x-axis.
    SphereVolume(double radius, Point location, Point z_axis, Point x_axis);

    /// Virtual destructor - ensures safe inheritance
    virtual ~SphereVolume();

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
    { return Class_SphereVolume; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual SphereVolume* clone() const
    { return new SphereVolume(radius_, location_, z_axis_, x_axis_); }


    // --- Functions inherited from ParamVolume ---

    DirectionCone tangentCone(int pardir) const;

    const Array<double,6> parameterSpan() const;

    void point(Point& pt, double upar, double vpar, double wpar) const;

    void point(std::vector<Point>& pts, 
	       double upar, double vpar, double wpar,
	       int derivs,
	       bool u_from_right = true,
	       bool v_from_right = true,
	       bool w_from_right = true,
	       double resolution = 1.0e-12) const;

    double nextSegmentVal(int dir, double par, bool forward, double tol) const;

    void closestPoint(const Point& pt,
		      double&        clo_u,
		      double&        clo_v,
		      double&        clo_w,
		      Point&         clo_pt,
		      double&        clo_dist,
		      double         epsilon,
		      double   *seed = 0) const;

    void reverseParameterDirection(int pardir);

    void swapParameterDirection(int pardir1, int pardir2);

    virtual std::vector<shared_ptr<ParamSurface> > 
	getAllBoundarySurfaces() const;

    virtual void translate(const Point& vec);

    // --- Functions inherited from ElementaryVolume ---

    SplineVolume* geometryVolume() const;


  private:

    double radius_;
    Point location_;
    Point z_axis_;
    Point x_axis_;
    Point y_axis_;

    void setCoordinateAxes();

  };    // Class SphereVolume



} // namespace Go



#endif    // #ifndef __SPHEREVOLUME_H
