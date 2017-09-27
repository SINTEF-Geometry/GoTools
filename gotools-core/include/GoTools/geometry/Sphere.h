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

#ifndef _SPHERE_H
#define _SPHERE_H


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Circle.h"


namespace Go
{


class SplineSurface;

/// \brief Class that represents a sphere. It is a subclass of
/// ElementarySurface, and thus has a parametrization and is
/// non-selfintersecting.
///
/// A Sphere has a natural parametrization in terms of the angles \a
/// u and \a v:
/// \b p(\a u, \a v)
/// = \b C + \a R cos \a v((cos \a u) \b x + sin(\a u) \b y) + \a R (sin\a v) \b z,
/// where \b C is a position vector, \a R is the radius, and \b x, \b
/// y and \b z are the (local) axes. The parametrization is
/// bounded by: \f$0 \leq u \leq 2\pi\f$, \f$-\frac{\pi}{2} < v < \frac{\pi}{2}\f$.
/// The dimension is 3.

class Sphere : public ElementarySurface
{
public:
    /// Default constructor. Constructs an uninitialized Sphere which
    /// can only be assigned to or read into.
    Sphere()
    {};

    /// Constructor. Input is the radius, the location, the direction
    /// of the z-axis and the (possibly approximate) x-axis. The local
    /// coordinate axes are normalized even if \c z_axis and/or \c
    /// x_axis are not unit vectors.
    Sphere(double radius, Point location, Point z_axis, Point x_axis,
        bool isSwapped = false);

    /// Virtual destructor - ensures safe inheritance
    virtual ~Sphere();

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
    { return Class_Sphere; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Sphere* clone() const;


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

    // Fetch the parameter curve in this sphere corresponding to a
    // given space curve provided the parameter curve is known to be 
    // represented by a line
    virtual shared_ptr<ElementaryCurve> 
      getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
			      const Point* start_par_pt = NULL, const Point* end_par_pt = NULL) const;

    bool isDegenerate(bool& b, bool& r,
		      bool& t, bool& l, double tolerance) const;


    /// Query if parametrization is bounded. All four parameter bounds
    /// must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded() const;

    /// Check if the sphere is closed
    bool isClosed(bool& closed_dir_u, bool& closed_dir_v) const;

    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    // --- Functions specific to Sphere ---

    /// Sphere radius
    double getRadius() const
    { return radius_; }
    
    /// Sphere centre
    Point getLocation() const
    { return location_;	}

    /// Local coordinate axes
    void getCoordinateAxes(Point& x_axis, Point& y_axis, Point& z_axis) const
    {
	x_axis = x_axis_;
	y_axis = y_axis_;
	z_axis = z_axis_;
    }

    /// Limit the surface by limiting the parameter domain
    virtual void setParameterBounds(double from_upar, double from_vpar,
				    double to_upar, double to_vpar);

    /// Pick part of sphere limited by parameter domain information
    Sphere* subSurface(double from_upar, double from_vpar,
		       double to_upar, double to_vpar,
		       double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// A SplineSurface representation of the sphere.
    virtual SplineSurface* geometrySurface() const;

    /// Create a SplineSurface representation of the sphere.
    virtual SplineSurface*  createSplineSurface() const;

    /// Get the circle along the latitude for a given v parameter.
    /// \param vpar v parameter
    /// \return A circle for the corresponding v parameter. If the v
    /// parameter is bounded, only a segment of a full circle is
    /// returned.
    shared_ptr<Circle> getLatitudinalCircle(double vpar) const;

    /// Get the circle along a longitude for a given u parameter.
    /// \param upar u parameter
    /// \return A circle for the corresponding u parameter. If the u
    /// parameter is bounded, only a segment of a full circle is
    /// returned.
    shared_ptr<Circle> getLongitudinalCircle(double upar) const;


    // Confirm that this surface is axis rotational
    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle);
    
    /// Radius in a specified location, 
    virtual double radius(double u, double v) const
    {
      return radius_;
    }

    virtual Point location() const
    {
      return location_;
    }

    virtual Point direction() const
    {
      return z_axis_;
    }

    virtual Point direction2() const
    {
      return x_axis_;
    }

    virtual void enlarge(double len1, double len2, double len3, double len4);

protected:

    double radius_;
    Point location_;
    Point z_axis_;
    Point x_axis_;
    Point y_axis_;

    RectDomain domain_;
    mutable RectDomain orientedDomain_; // Takes isSwapped_ into account

    void setCoordinateAxes();

};

} // namespace Go


#endif // _SPHERE_H

