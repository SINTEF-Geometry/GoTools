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

#ifndef _TORUS_H
#define _TORUS_H


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/Circle.h"


namespace Go
{


class SplineSurface;


/// \brief Class that represents a torus. It is a subclass of
/// ElementarySurface, and thus has a parametrization. A torus may be
/// degenerate. Then the minor radius is greater than the major
/// radius, and it is in principle selfintersecting.
///
/// A nondegenerate Torus has a natural parametrization in terms of
/// the two angles \a u and \a v:
/// \b p(\a u, \a v)
/// = \b C + (\a R + \a r cos \a v)((cos \a u) \b x + (sin \a u) \b y) + (\a r sin \a v) \b z,
/// where \b C is a position vector, \a R is the major radius, \a r is
/// the minor radius, and \b x, \b y and \b z are the (local) axes.
/// The parametrization is bounded by: \f$0 \leq u, v \leq 2\pi\f$.
/// The dimension is 3.
///
/// If the Torus is degenerate, there is an additional boolean
/// parameter, \c select_outer, which selects the inner or outer
/// (non-selfintersecting) part of the surface. This restricts the
/// parametrization by: \f$ - \phi \leq v \leq \phi \f$,
/// if \c select_outer is \c true, and
/// \f$\phi \leq v \leq 2\pi - \phi \f$, if \c select_outer is \c
/// false. \c select_outer affects the direction of the normal, which
/// always points outwards.

class Torus : public ElementarySurface
{
public:
    /// Default constructor. Constructs an uninitialized Torus which
    /// can only be assigned to or read into.
    Torus()
    {};

    /// Constructor. Input is the major and minor radii, the location,
    /// the direction of the z-axis and the (possibly approximate)
    /// x-axis. The local coordinate axes are normalized even if \c
    /// z_axis and/or \c x_axis are not unit vectors. If degenerate,
    /// one may also provide the select_outer flag.
    Torus(double major_radius, double minor_radius, 
	  Point location, Point z_axis, Point x_axis,
	  bool select_outer = true, bool isSwapped = false);

    /// Virtual destructor - ensures safe inheritance
    virtual ~Torus();

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
    { return Class_Torus; }

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Torus* clone() const;

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

    /// The normal to the torus. The normal always points \a outwards
    /// from the surface, as seen from the circle of radius \c
    /// major_radius_. Thus, if the torus is degenerate and \c
    /// select_outer_ is \c false, then the direction of the normal is
    /// opposite to
    /// \f$\partial p/ \partial u\times\partial p\partial v\f$.
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


    /// Query if parametrization is bounded. All four parameter bounds
    /// must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    bool isBounded() const;

    /// Check if the surface is closed
    bool isClosed(bool& closed_dir_u, bool& closed_dir_v) const;

    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    // --- Functions specific to Torus ---

    /// Fetch major radius
    double getMajorRadius() const
    { return major_radius_; }
    
    /// Fetch major radius
    double getMinorRadius() const
    { return minor_radius_; }
    
    /// Fetch centre 
    Point getLocation() const
    { return location_;	}

    /// Local coordinate axes
    void getCoordinateAxes(Point& x_axis, Point& y_axis, Point& z_axis) const
    {
	x_axis = x_axis_;
	y_axis = y_axis_;
	z_axis = z_axis_;
    }

    /// Fetch the \c select_outer flag. If the torus is not degenerate,
    /// this function has no effect. If the torus is degenerate, and
    /// the value of the flag is changed, this will change the
    /// parameter domain of the surface to the appropriate default.
   bool getSelectOuter() const
    { return select_outer_; }

    /// Set the \c select_outer flag. If the torus is not degenerate,
    /// this function has no effect. If the torus is degenerate, and
    /// the value of the flag is changed, this will change the
    /// parameter domain of the surface to the appropriate default.
    void setSelectOuter(bool select_outer);

    /// arccos(-major_radius/minor_radius). Used in degeneracy checks and to
    /// set the parameter domain of the surface
    double getPhi() const
    { return phi_; }

    /// Query if the Torus is degenerate. Note that a Torus may be
    /// degnerate, although the parameter domain might not represent a
    /// degnerate patch. Thus, isDegenerateTorus() might return \c
    /// true and isDegerate() might return \c false at the same time.
    bool isDegenerateTorus() const
    { return is_degenerate_torus_; }

    /// Limit the surface by limiting the parameter domain
    virtual void setParameterBounds(double from_upar, double from_vpar,
				    double to_upar, double to_vpar);

    /// Return the part of the torus surface limited by the parameter bounds
    Torus* subSurface(double from_upar, double from_vpar,
		      double to_upar, double to_vpar,
		      double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Create a SplineSurface representation of the Torus.
    virtual SplineSurface* geometrySurface() const;

    /// Create a SplineSurface representation of the Torus.
    virtual SplineSurface*  createSplineSurface() const;

    /// Get the major circle for a given v parameter.
    /// \param vpar v parameter
    /// \return A circle for the corresponding v parameter. If the v
    /// parameter is bounded, only a segment of a full circle is
    /// returned.
    shared_ptr<Circle> getMajorCircle(double vpar) const;

    /// Get the minor circle for a given u parameter.
    /// \param upar u parameter
    /// \return A circle for the corresponding u parameter. If the u
    /// parameter is bounded, only a segment of a full circle is
    /// returned.
    shared_ptr<Circle> getMinorCircle(double upar) const;

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

    double major_radius_;
    double minor_radius_;
    Point location_;
    Point z_axis_;
    Point x_axis_;
    Point y_axis_;
    bool is_degenerate_torus_; // Function of the radii
    bool select_outer_;
    double phi_; // Function of the radii

    RectDomain domain_;
    mutable RectDomain orientedDomain_; // Takes isSwapped_ into account

    void setCoordinateAxes();
    void setDegenerateInfo();
    void setDefaultDomain();

};

} // namespace Go


#endif // _TORUS_H

