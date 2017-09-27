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
  /// \b p(\a r, \a v) = \b C + \a r((cos \a v) \b x + (sin \a v) \b y),
  /// where \b C is the centre position vector and \b x and \b y are
  /// the (local) axes. The parametrization is bounded by:
  /// \f$0 \leq r \leq R\f$ and \f$0 \leq v \leq 2\pi\f$, where \a R is the
  /// disc radius. The dimension is 2 or 3.
  ///
  /// A disc also holds degeneracy information for representing it as
  /// a SplineSurface object. There are two ways to parametrize:
  /// 
  /// 1. <b>Polar coordinates:</b>
  /// The parametrization of the
  /// SplineSurface object almost coincides with the one for the Disc
  /// object itself, i.e. this->point() gives almost the same as
  /// geometrySurface()->point(). The only degeneracy point then is
  /// the centre. Note that the corresponing SplineSurface
  /// representation is not an identity mapping of the Disc because
  /// of a reparametrization in the angular direction. This is the
  /// default SplineSurface representation.
  /// 
  /// 2. <b>Rounded square:</b>
  /// Done by splitting the boundary into
  /// four curves that become the boundary curves of the spline
  /// surface. Then the four meeting points of the curves
  /// become degeneracy points.

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
    Disc(Point centre, double radius, Point x_axis, Point normal,
        bool isSwapped = false);

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

    bool isDegenerate(bool& b, bool& r,
		      bool& t, bool& l, double tolerance) const;


    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    /// Query if parametrization is bounded. All four parameter bounds
    /// must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded() const;

    /// Check if the surface is closed.
    bool isClosed(bool& closed_dir_u, bool& closed_dir_v) const;

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

    /// The NURBS representation of the disc will have degeneracy in the centre
    void useCentreDegen()
    { centre_degen_ = true; }

    /// The NURBS representation of the disc will have degenerate corners
    void useCornerDegen()
    { centre_degen_ = false; }

    virtual void enlarge(double len1, double len2, double len3, double len4);

    virtual Point location() const
    {
      return centre_;
    }

    virtual Point direction() const
    {
      return z_axis_;
    }

    virtual Point direction2() const
    {
      return x_axis_;
    }


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
    mutable RectDomain orientedDomain_; // Takes isSwapped_ into account

    void setCoordinateAxes();
    void setDefaultDomain();
    Circle boundaryCircle() const;

  };    // Class Disc


} // namespace Go



#endif    // #ifndef __DISC_H
