//===========================================================================
//
// File : CylinderVolume.h
//
// Created: Tue Nov 10 08:10:21 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#ifndef __CYLINDERVOLUME_H
#define __CYLINDERVOLUME_H



#include "GoTools/trivariate/ElementaryVolume.h"


namespace Go
{


  class SplineVolume;


  /// \brief Class that represents a solid cylinder, maybe with empty interior
  /// like a straight tube. It is a subclass of
  /// ElementaryVolume, and has a natural parametrization in terms of a radius \a u,
  /// an angle \a v, and distance \a w:
  /// \b p(\a u, \a v, \a w)
  /// = \b C + \a u(cos(\a v) \b x + sin(\a v) \b y) + \a w \b z,
  /// where \b C is a position vector, and \b x, \b y and \b z are the (local) axes.
  /// The parametrization is bounded by: \f$R \leq u \leq S\f$,
  /// \f$0 \leq v \leq 2\pi\f$, \f$-\infty < w < \infty\f$, where \f$ R and \f$ S
  /// are the inner and outer radius. The dimension is 3.

  class CylinderVolume : public ElementaryVolume
  {
  public:
    /// Default constructor. Constructs an uninitialized CylinderVolume which
    /// can only be assigned to or read into.
    CylinderVolume() { }

    /// Constructor. Input is the centre, outer radius, the normal and the
    /// x_axis (possibly approximate). The inner radius is set to 0, and
    /// the cylinder is set to have infinite length in both directions
    CylinderVolume(Point centre, double radius, Point normal, Point x_axis);

    /// Virtual destructor - ensures safe inheritance
    virtual ~CylinderVolume();

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
    { return Class_CylinderVolume; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual CylinderVolume* clone() const;


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

    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v,
			      double&        clo_w,
			      Point&         clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      double   *seed = 0) const;

    void reverseParameterDirection(int pardir);

    void swapParameterDirection(int pardir1, int pardir2);

    virtual std::vector<std::shared_ptr<ParamSurface> > 
	getAllBoundarySurfaces() const;

    virtual void translate(const Point& vec);

    // --- Functions inherited from ElementaryVolume ---

    SplineVolume* geometryVolume() const;


    // --- Own functions ---

    void setParameters(double from_par, double to_par, int pardir);

    void useCentreDegen()
    { centre_degen_ = true; }

    void useCornerDegen()
    {
      centre_degen_ = false;
      radius_min_ = 0.0;
    }

  private:

    Point centre_;
    double radius_min_, radius_max_;
    double height_min_, height_max_;
    Point
      x_axis_,
      y_axis_,
      z_axis_;
    bool centre_degen_;  // If true, geometryVolume() gives a SplineVolume with
                         // degenerecy at the normal axis (if radius_min_ = 0.0)
                         // If false, let boundary curves of the SplineSurface
                         // lie on the boundary. In this case, radius_min_ = 0.0,
                         // and the cylinder is considered to be solid in the centre
    double degen_angles_[4];  // The angle parameter value giving the four degeneracy
                              // points on the boundary (only when center_degen = false)

    void setCoordinateAxes();

  };    // Class CylinderVolume



} // namespace Go



#endif    // #ifndef __CYLINDERVOLUME_H
