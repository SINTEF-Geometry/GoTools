//===========================================================================
//
// File : Parallelepiped.h
//
// Created: Fri Nov  6 14:21:19 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#ifndef __PARALLELEPIPED_H
#define __PARALLELEPIPED_H



#include "GoTools/trivariate/ElementaryVolume.h"


namespace Go
{


  class SplineVolume;


  /// \brief Class that represents a solid parallelepiped. It is a subclass of
  /// ElementaryVolume, and thus has a parametrization.
  ///
  /// A Parallelepiped has a natural parametrization in terms
  /// of the distances \a u, \a v and \a w along the edges:
  /// \b p(\a u, \a v, \a w) = \b C + \a u \b x + \a v \b y + \a w \b z.
  /// The parametrization is bounded by: \f$0 \leq u \leq a\f$, \f$0 \leq v \leq b\f$,
  /// \f$0 \leq w \leq c\f$, where \a a, \a b, \a c are the edge lengths. The dimension is 3.

  class Parallelepiped : public ElementaryVolume
  {
  public:
    /// Default constructor. Constructs an uninitialized Parallelepiped which
    /// can only be assigned to or read into.
    Parallelepiped() { }

    /// Constructor
    Parallelepiped(Point corner,
		   Point dir_u, Point dir_v, Point dir_w,
		   double len_u, double len_v, double len_w);

    /// virtual destructor - ensures safe inheritance
    virtual ~Parallelepiped();

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
    { return Class_Parallelepiped; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Parallelepiped* clone() const
    { return new Parallelepiped(corner_, dir_u_, dir_v_, dir_w_, length_u_, length_v_, length_w_); }


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

    virtual std::vector<std::shared_ptr<ParamSurface> > 
	getAllBoundarySurfaces() const;

    virtual void translate(const Point& vec);

    // --- Functions inherited from ElementaryVolume ---

    SplineVolume* geometryVolume() const;



  private:
    Point corner_;                         // Lower, left, closest corner
    Point dir_u_, dir_v_, dir_w_;         // Edge directions, normalized
    double length_u_, length_v_, length_w_; // Edge lengths



  };    // Class Parallelepiped



} // namespace Go



#endif    // #ifndef __PARALLELEPIPED_H
