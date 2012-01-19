//===========================================================================
//
// File : TorusVolume.h
//
// Created: Wed Nov 11 08:26:10 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#ifndef __TORUSVOLUME_H
#define __TORUSVOLUME_H



#include "GoTools/trivariate/ElementaryVolume.h"


namespace Go
{


  class SplineVolume;


  /// \brief Class that represents a solid torus, maybe with empty interior
  /// like a circular pipe, and/or a section (not full revolution along the
  /// main axis). It is a subclass of ElementaryVolume, and has a natural
  /// parametrization in terms of a radius \a u and two angles \a v and \a w:
  /// \b p(\a u, \a v, \a w)
  /// = \b C + (\a R + \a u cos \a w)((cos \a v) \b x + (sin \a v) \b y) + (\a u sin \a w) \b z,
  /// where \b C is a position vector, \a R is the major radius, and \b x,
  /// \b y and \b z are the (local) axes. The parametrization for a full
  /// solid torus is bounded by: \f$0 \leq u \leq r\f$,
  /// \f$0 \leq v, w \leq 2\pi\f$, where \a r is the minor axis.
  /// A bended pipe will have a positive minimal limit for \a u, while a
  /// torus segment will have other limits for \a v. The dimension is 3.

  class TorusVolume : public ElementaryVolume
  {
  public:
    /// Default constructor. Constructs an uninitialized TorusVolume which
    /// can only be assigned to or read into.
    TorusVolume() { }

    /// Constructor. Input is the major and minor radii, the location,
    /// the direction of the z-axis and the (possibly approximate)
    /// x-axis.
    TorusVolume(double major_radius, double minor_radius, 
		Point location, Point z_axis, Point x_axis);

    /// Virtual destructor - ensures safe inheritance
    virtual ~TorusVolume();

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
    { return Class_TorusVolume; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual TorusVolume* clone() const;


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


    // --- Own functions ---
    /// Restrict the size of the torus in one parameter direction
    void setParameters(double from_par, double to_par, int pardir);

    /// A NURBS representation will be a disc with degenereracy in the centre
    /// rotated around the centre
    void useCentreDegen()
    { centre_degen_ = true; }

    /// A NURBS representation will be a disc with degenerate corners rotated
    /// around the centre
    void useCornerDegen()
    {
      centre_degen_ = false;
      minor_radius_min_ = 0.0;
    }

  private:

    Point location_;
    double major_radius_;
    double minor_radius_min_, minor_radius_max_;
    double rev_angle_min_, rev_angle_max_;
    Point
      z_axis_,
      x_axis_,
      y_axis_;
    bool centre_degen_;  // If true, geometryVolume() gives a SplineVolume with
                         // degenerecy at the normal axis (if radius_min_ = 0.0)
                         // If false, let boundary curves of the SplineSurface
                         // lie on the boundary. In this case, radius_min_ = 0.0,
                         // and the torus is considered to be solid in the centre
    double degen_angles_[4];  // The angle parameter value giving the four degeneracy
                              // points on the boundary (only when center_degen = false)

    void setCoordinateAxes();

  };    // Class TorusVolume



} // namespace Go



#endif    // #ifndef __TORUSVOLUME_H
