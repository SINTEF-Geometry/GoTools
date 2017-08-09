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

#ifndef _SURFACEOFREVOLUTION_H
#define _SURFACEOFREVOLUTION_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/Circle.h"


namespace Go {


class SplineCurve;
class SplineSurface;


/// \brief Class that represents a surface of revolution. A
/// SurfaceOfRevolution is swept out by a SplineCurve that is rotated
/// around an axis with a complete revolution, and is thereby a
/// parametric surface. 
///
/// The parametrization of the surface is given in terms of a location
/// \b C, an axis line \b V, and the spline curve \b \f$\lambda(v)\f$
/// with parameter \a v:
///
/// \f$ \sigma(u, v)
/// = C + (\lambda(v) - C) \cos u
/// + ((\lambda(v) - C) \cdot V) V (1 - \cos u)
/// + V \times (\lambda(v) - C) \sin u
/// \f$
///
/// The parameter \a u is bounded by: \f$0 \leq u \leq 2\pi\f$. The
/// axis \b V is normalized.
/// 
/// The curve \f$\lambda\f$ must be such that it doesn't lead to a
/// self-intersecting surface.


class SurfaceOfRevolution : public ParamSurface
{
public:
    /// Default constructor. Constructs an uninitialized
    /// SurfaceOfRevolution which can only be assigned to or read
    /// into.
    SurfaceOfRevolution()
    {};

    /// Constructor. Input is the location and normalized direction of
    /// the axis, and the SplineCurve that is swept out by the
    /// revolution.
    SurfaceOfRevolution(Point location, Point axis_dir,
			shared_ptr<SplineCurve> curve,
                        bool isSwapped = false);

    /// Virtual destructor - ensures safe inheritance
    virtual ~SurfaceOfRevolution();

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
    static ClassType classType();

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual SurfaceOfRevolution* clone() const;


    // --- Functions inherited from ParamSurface ---

    const RectDomain& parameterDomain() const;
    virtual RectDomain containingDomain() const;

    /// Check if a parameter pair lies inside the domain of this surface
    virtual bool inDomain(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies inside the domain of this surface
    /// return value = 0: outside
    ///              = 1: internal
    ///              = 2: at the boundary
    virtual int inDomain2(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies at the boundary of this surface
    virtual bool onBoundary(double u, double v, double eps=1.0e-4) const;

    /// Return the parameter value in the domain of this surface closest
    /// to the parameter pair (u,v).
    virtual Point closestInDomain(double u, double v) const;

    CurveLoop outerBoundaryLoop(double degenerate_epsilon
				= DEFAULT_SPACE_EPSILON) const;
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

    void turnOrientation();

    void reverseParameterDirection(bool direction_is_u);

    void swapParameterDirection();

    virtual double area(double tol) const;

    //bool isDegenerate(bool& b, bool& r,
		  //    bool& t, bool& l, double tolerance) const;


    /// Check for parallel and anti parallel partial derivatives in
    /// surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, 
				      double tol) const;


    /// Return surface corners, geometric and parametric points in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    // --- Functions specific to SurfaceOfRevolution ---

    /// Point on axis of revolution
    Point getLocation() const
    { return location_;	}

    /// Direction of axis of revolution
    Point getAxisDir() const
    { return axis_dir_;	}

    /// Generating curve, i.e. the curve that are rotated around the axis
    shared_ptr<SplineCurve> getCurve() const
    { return curve_; }
    
    /// Limit the surface by limiting the parameter domain
    void setParameterBounds(double from_upar, double from_vpar,
			    double to_upar, double to_vpar);

    /// Pick part of surface
    SurfaceOfRevolution* subSurface(double from_upar, double from_vpar,
				    double to_upar, double to_vpar,
				    double fuzzy
				    = DEFAULT_PARAMETER_EPSILON) const;

    // Is "geometrySurface()" a good name for this function? @jbt
    /// Return spline representation of the surface of revolution
    SplineSurface* geometrySurface() const;
    /// Create a SplineSurface representation of the surface of revolution
    SplineSurface* createSplineSurface() const;

    /// Get the circle for a given v parameter.
    /// \param vpar v parameter
    /// \return A circle for the corresponding v parameter. If the v
    /// parameter is bounded, only a segment of a full circle is
    /// returned.
    shared_ptr<Circle> getCircle(double vpar) const;


    virtual bool isSwapped() const;

private:
    Point location_;
    Point axis_dir_;
    shared_ptr<SplineCurve> curve_;

    RectDomain domain_;
    mutable RectDomain orientedDomain_; // Takes isSwapped_ into account

    bool isSwapped_;

    // Helper function to be used in functions like point(), where
    // we need to take the isSwapped_ flag into account.
    void getOrientedParameters(double& u, double&v) const;

};


} // namespace Go


#endif // _SURFACEOFREVOLUTION_H
