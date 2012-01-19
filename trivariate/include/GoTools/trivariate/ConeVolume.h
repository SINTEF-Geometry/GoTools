//===========================================================================
//
// File : ConeVolume.h
//
// Created: Tue Nov 10 12:12:09 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#ifndef __CONEVOLUME_H
#define __CONEVOLUME_H



#include "GoTools/trivariate/ElementaryVolume.h"


namespace Go
{


  class SplineVolume;


  /// \brief Class that represents a solid cone. It is a subclass of
  /// ElementaryVolume, and has a natural parametrization in terms of a radius \a u,
  /// an angle \a v, and distance \a w:
  /// \b p(\a u, \a v, \a w)
  /// = \b C + u (\a R + \a w tan \f$\alpha\f$)((cos \a v) \b x + (sin \a v) \b y) + \a w \b z,
  /// where \b C is the cone apex, \a R is the radius when \a w = 0, \f$\alpha\f$
  /// is the cone angle, and \b x, \b y and \b z are the (local) axes.
  /// The parametrization is bounded by: \f$0 \leq u \leq 1\f$,
  /// \f$0 \leq v \leq 2\pi\f$, \f$-\infty < w < \infty\f$. The dimension is 3.

  class ConeVolume : public ElementaryVolume
  {
  public:
    /// Default constructor. Constructs an uninitialized ConeVolume which
    /// can only be assigned to or read into.
    ConeVolume() { }

    /// Constructor. Input is the radius, the location, the direction of
    /// the z-axis and the (possibly approximate) x-axis.
    ConeVolume(double radius, Point location,
	       Point z_axis, Point x_axis,
	       double cone_angle);

    /// Virtual destructor - ensures safe inheritance
    virtual ~ConeVolume();

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
    { return Class_ConeVolume; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual ConeVolume* clone() const;


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

    std::vector<shared_ptr<ParamSurface> > 
	getAllBoundarySurfaces() const;

    virtual void translate(const Point& vec);

    // --- Functions inherited from ElementaryVolume ---

    SplineVolume* geometryVolume() const;


    // --- Own functions ---

    /// Restrict the size of the torus in one parameter direction
    void setParameters(double from_par, double to_par, int pardir);

     /// A NURBS representation will be a disc with degenereracy in the centre
    /// swept linearily and shrinked or extended
   void useCentreDegen()
    { centre_degen_ = true; }

    /// A NURBS representation will be a disc with degenerate corners 
    /// swept linearily and shrinked or extended. Note that extra degeneracies
   /// occur if the cone tip is included in the volume.
    void useCornerDegen()
    {
      centre_degen_ = false;
    }


  private:

    double radius_;
    Point location_;
    Point z_axis_;
    Point x_axis_;
    Point y_axis_;
    double cone_angle_;

    double height_min_, height_max_;

    bool centre_degen_;  // If true, geometryVolume() gives a SplineVolume with
                         // degenerecy at the normal axis.
                         // If false, let boundary curves of the SplineSurface
                         // lie on the boundary.
    double degen_angles_[4];  // The angle parameter value giving the four degeneracy
                              // points on the boundary (only when center_degen = false)

    void setCoordinateAxes();

  };    // Class ConeVolume



} // namespace Go



#endif    // #ifndef __CONEVOLUME_H
