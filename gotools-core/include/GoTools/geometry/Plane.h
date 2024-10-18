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

#ifndef _PLANE_H
#define _PLANE_H


#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/utils/RotatedBox.h"


namespace Go
{


class SplineSurface;


/// \brief Class that represents a plane. It is a subclass of
/// ElementarySurface, and thus has a parametrization and is
/// non-selfintersecting.
///
/// A Plane has a natural parametrization in terms of its location \b
/// C and spanning vectors \b x and \b y: p(u, v) = C + ux + vy.  This
/// parametrization might be unbounded: -\f$\infty < u,v < \infty\f$.

class Plane : public ElementarySurface
{
public:
    /// Default constructor. Constructs an uninitialized Plane which
    /// can only be assigned to or read into.
    Plane()
    {};

    /// Constructor. Input is location and normal
    Plane(Point location, Point normal,
        bool isSwapped = false);

    /// Constructor. Input is location, normal and (local,
    /// approximate) x-axis
    Plane(Point location, Point normal, Point x_axis,
        bool isSwapped = false);

    /// Constructor. Input is coefficients of implicit equation +
    /// point suggestion
    Plane(double a, double b, double c, double d,
        bool isSwapped = false);

    /// Copy constructor
    Plane& operator= (const Plane& other);

    /// virtual destructor - ensures safe inheritance
    virtual ~Plane();

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
    { return Class_Plane; }

    /// Return empty box if infinite plane
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual Plane* clone() const;


    // --- Functions inherited from ParamSurface ---

    const RectDomain& parameterDomain() const;

    virtual DirectionCone normalCone() const;
    virtual DirectionCone tangentCone(bool pardir_is_u) const;

    virtual void point(Point& pt, double upar, double vpar) const;
    virtual void point(std::vector<Point>& pts, 
               double upar, double vpar,
               int derivs,
               bool u_from_right = true,
               bool v_from_right = true,
               double resolution = 1.0e-12) const;

    virtual void normal(Point& n, double upar, double vpar) const;

    virtual std::vector<shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const;

    virtual std::vector<shared_ptr<ParamSurface> >
    subSurfaces(double from_upar, double from_vpar,
                double to_upar, double to_vpar,
                double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    virtual void closestPoint(const Point& pt,
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

    virtual bool isDegenerate(bool& b, bool& r,
			      bool& t, bool& l, double tolerance) const;


    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    virtual shared_ptr<ElementaryCurve> 
      getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
			      const Point* start_par_pt = NULL, const Point* end_par_pt = NULL) const;
    // --- Functions specific to Plane ---

    /// Point in plane
    Point getPoint()
    { return location_; }

    /// Plane normal. NB: This function returns the "defining normal" and
    /// does not take swapped parameter directions into account! To get the
    /// "oriented normal", use Point::normal()
    Point getNormal()
    { return normal_; }
    
    /// Vectors in plane
    void getSpanningVectors(Point& axis1, Point& axis2)
    {
        axis1 = vec1_;
        axis2 = vec2_;
    }

    /// Projection of the point pnt in the plane
    Point projectPoint(const Point& pnt) const;

    /// Distance between the point pnt and the plane
    double distance(const Point& pnt) const;

    /// Restrict the plane by restricting the parameter domain. It is
    /// initially infinite
    virtual void setParameterBounds(double from_upar, double from_vpar,
                            double to_upar, double to_vpar);

    /// Fetch parameter bounds. NB! Not oriented
    virtual RectDomain getParameterBounds() const
    {
      return parbound_;
    }

    /// set the parameter domain to a given rectangle
    /// \param u1 new min. value of first parameter span
    /// \param u2 new max. value of first parameter span
    /// \param v1 new min. value of second parameter span
    /// \param v2 new max. value of second parameter span
    virtual void setParameterDomain(double u1, double u2, double v1, double v2);

    /// Fetch a part of the plane
    Plane* subSurface(double from_upar, double from_vpar,
                      double to_upar, double to_vpar,
                      double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Create a SplineSurface representation of the Plane.
    virtual SplineSurface* geometrySurface() const;

    /// Create a SplineSurface representation of the Plane.
    virtual SplineSurface*  createSplineSurface() const;

    /// Create a non-rational spline surface from the surface
    virtual SplineSurface* createNonRationalSpline(double eps) const;

    /// Query if parametrization is bounded. All four parameter bounds
    /// must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    virtual bool isBounded() const;

     /// Check if the plane is closed. Virtual function - always false.
    virtual bool isClosed(bool& closed_dir_u, bool& closed_dir_v) const;

    /// Return the result from intersecting the unbounded plane with a
    /// rotated bounding box (having axis[0]=vec1_, axis[1]=vec2_,
    /// axis[2]=normal_). Useful for visualizing the (unbounded) plane.
    /// If intersection is empty, the returned plane is the NULL
    /// pointer.  The rotated box may be created from a boundingbox by
    /// defining the coordinate system and the 8 corner points of the
    /// bd_box.
    Plane* intersect(const RotatedBox& bd_box) const;

    /// Inherited from elementary surface
    /// Confirm that the surface is a plane and return the plane normal
    virtual bool isPlanar(Point& normal, double tol);

    /// The surface is linear in all directions. Fetch the 
    /// parameter directions
    virtual bool isLinear(Point& dir1, Point& dir2, double tol);

    virtual Point location() const
    {
      return location_;
    }

    virtual Point direction() const
    {
      return normal_;
    }

    virtual Point direction2() const
    {
      return vec1_;
    }

    virtual void enlarge(double len1, double len2, double len3, double len4);

    virtual void translate(const Point& vec);

protected:

    Point location_;

    // The vectors vec1_, vec2_, and normal_ define a right-handed
    // coordinate system.
    Point normal_;
    Point vec1_;
    Point vec2_;

    RectDomain parbound_;
    RectDomain domain_;
    mutable RectDomain orientedDomain_; // Takes isSwapped_ flag into account

    // vec1_ is projected onto plane (defined by normal_), vec2_ set accordingly.
    void setSpanningVectors();

    // This version recomputes the vec1_ vector. Avoids tolerance issues.
    void setSpanningVectorsSafe();

};

} // namespace Go



#endif // _PLANE_H

