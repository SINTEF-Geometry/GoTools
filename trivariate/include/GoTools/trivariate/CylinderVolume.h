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

    virtual std::vector<shared_ptr<ParamSurface> > 
	getAllBoundarySurfaces() const;

    virtual void translate(const Point& vec);

    // --- Functions inherited from ElementaryVolume ---

    SplineVolume* geometryVolume() const;


    // --- Own functions ---

    /// Restrict the size of the torus in one parameter direction
    void setParameters(double from_par, double to_par, int pardir);

    /// A NURBS representation will be a disc with degenereracy in the centre
    /// swept linearily 
   void useCentreDegen()
    { centre_degen_ = true; }

    /// A NURBS representation will be a disc with degenerate corners 
    /// swept linearily 
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
