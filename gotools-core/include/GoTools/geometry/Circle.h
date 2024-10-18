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

#ifndef _CIRCLE_H
#define _CIRCLE_H


#include "GoTools/geometry/ElementaryCurve.h"


namespace Go {


/// \brief Class that represents a circle. It is a subclass of
/// ElementaryCurve and thus has a parametrization.
///
/// A Circle has a natural parametrization in terms of the angle \a t:
/// \b p(\a t) = \b C + \b r(cos(\a t) \b x + sin(\a t) \b y), where
/// \b C is the centre, \b r is the radius, and \b x and \b y are the
/// (local) axes. The parametrization is (in principle) bounded:
/// \f$0 \leq t \leq 2\pi\f$. The dimension is either 2 or 3.

class Circle : public ElementaryCurve
{
public:
    /// Default constructor. Constructs an uninitialized Circle which
    /// can only be assigned to or read into.
  Circle();

    /// Constructor. Input is the radius, the centre, the normal to the
    /// plane of the circle, and the (approximate) direction of the
    /// x-axis. The dimension must be either 2 or 3.
    Circle(double radius,
	   Point centre, Point normal, Point x_axis,
           bool isReversed = false);

    /// Copy constructor
    Circle& operator= (const Circle& other);
    
    /// virtual destructor - ensures safe inheritance
    virtual ~Circle();

    /// Read object from stream
    /// \param is stream from which object is read
    virtual void read (std::istream& is);
    /// Write object to stream
    /// \param os stream to which object is written
    virtual void write (std::ostream& os) const;


    // --- Functions inherited from GeomObject ---

    virtual BoundingBox boundingBox() const;
    
    virtual int dimension() const;

    virtual ClassType instanceType() const;

    static ClassType classType();

    virtual Circle* clone() const;

    // --- Functions inherited from ParamCurve ---

    virtual void point(Point& pt, double tpar) const;

    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    virtual double startparam() const;
    virtual double endparam() const;

    //virtual void reverseParameterDirection(bool switchparam = false);
    
    /// The full circle is always parametrized on [0, 2*M_PI), it does
    /// not make sense to reparametrize. For picking a subset of the
    /// range use setParamBounds() instead.
    virtual void setParameterInterval(double t1, double t2);

    virtual SplineCurve* geometryCurve();

    /// Creates a SplineCurve representation of the circle. The spline
    /// representation is geometrically exact, but the parametrization
    /// is only an approximation. If setParameterBounds() have been
    /// used to set bounds on the parametrization, then the endpoints
    /// of the resulting spline curve \em does have correct parameter
    /// values.
    virtual SplineCurve* createSplineCurve() const;

    /// Approximating with a non-rational spline curve
    SplineCurve* createNonRationalSpline(double eps) const;

    virtual bool isDegenerate(double degenerate_epsilon);

    virtual 
    Circle* subCurve(double from_par, double to_par,
		   double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    virtual DirectionCone directionCone() const;
 
    virtual void appendCurve(ParamCurve* cv, bool reparam=true);

    virtual void appendCurve(ParamCurve* cv,
			     int continuity, double& dist, bool reparam=true,
			     double tol = 1.0e-4);

    virtual void closestPoint(const Point& pt,
			      double tmin,
			      double tmax,
			      double& clo_t,
			      Point& clo_pt,
			      double& clo_dist,
			      double const *seed = 0) const;

    virtual double length(double tol);

    // --- Functions specific to Circle ---

    Point getCentre() const
    {
      return centre_;
    }

    Point getNormal() const
    {
      return normal_;
    }

    Point getXAxis() const
    {
      return vec1_;
    }

    double getRadius() const
    {
      return radius_;
    }

    /// Set bounds for the parametrization of the Circle other than
    /// the default \f$[0, 2\pi]\f$. Requirements for a valid
    /// parametrization are: 1) The first parameter must be strictly
    /// less than the second, 2) Parameter values must be in
    /// \f$[-2\pi, 2\pi]\f$, and 3) (\e endpar - \e startpar) must not
    /// exceed \f$2\pi\f$.
    /// \param startpar start parameter
    /// \param endpar end parameter
    virtual void setParamBounds(double startpar, double endpar);

    // Translate the curve along a given vector
    virtual void translateCurve(const Point& dir);

    bool isClosed() const;

    // Confirm that the curve is axis rotational
    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle);

    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle, double& radius);

    /// If the curve is 2 dimensional, x and y coordinates will be swapped.
    /// Used when curve is a parameter curve.
    virtual void swapParameters2D();

    /// Check if the line lies in a plane passing through a given axis
    virtual bool isInPlane(const Point& loc, const Point& axis,
			   double eps, Point& normal) const;

    /// Check if the circle lies in a plane with a given normal
    virtual bool isInPlane(const Point& norm,
			   double eps, Point& pos) const;

protected:

    double radius_;
    Point centre_;
    Point normal_;
    Point vec1_;
    Point vec2_;

    double parbound1_;
    double parbound2_;
    double startparam_;
    double endparam_;

    void setSpanningVectors();
};


} // namespace Go


#endif // _CIRCLE_H

