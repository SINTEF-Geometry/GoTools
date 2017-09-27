/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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

    virtual std::vector<shared_ptr<ParamSurface> > 
	getAllBoundarySurfaces(bool do_clear = false) const;

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
