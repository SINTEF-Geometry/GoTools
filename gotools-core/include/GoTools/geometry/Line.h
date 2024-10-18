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

#ifndef _LINE_H
#define _LINE_H


#include "GoTools/geometry/ElementaryCurve.h"


namespace Go {


/// \brief Class that represents a line. It is a subclass of
/// ElementaryCurve and thus has a parametrization.
///
/// A Line has a natural parametrization in terms of its location \b
/// C and direction vector \b V: p(t) = C + tV.

class Line : public ElementaryCurve
{
public:
    /// Default constructor. Constructs an uninitialized Line which
    /// can only be assigned to or read into.
    Line()
    {};

    /// Constructor. Input is point and direction. The default Line is
    /// unbounded in its parametrization. To bound it, use
    /// setParamBounds().
    Line(Point point, Point direction, bool isReversed = false);

    /// Constructor. Bounded line
    Line(Point point, Point direction, double length,
        bool isReversed = false);

    /// Constructor. Bounded line with parameterization
    Line(Point point1, Point point2, double par1, double par2,
        bool isReversed = false);

    /// Copy constructor
    Line& operator= (const Line& other);
    
    /// virtual destructor - ensures safe inheritance
    virtual ~Line();

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

    virtual Line* clone() const;

    // --- Functions inherited from ParamCurve ---

    virtual void point(Point& pt, double tpar) const;

    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    virtual double startparam() const;
    virtual double endparam() const;

    //virtual void reverseParameterDirection(bool switchparam = false);
    
     /// Limit the curve by limiting the parameter interval
   virtual void setParameterInterval(double t1, double t2);

    virtual SplineCurve* geometryCurve();
    /// Fetch spline representation of curve
    virtual SplineCurve* createSplineCurve() const;

    virtual bool isDegenerate(double degenerate_epsilon);

    virtual 
    Line* subCurve(double from_par, double to_par,
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

    // --- Functions specific to Line ---

    /// Point on line
    Point getPoint()
    { return location_; }

    /// Direction vector
    Point getDirection()
    { return dir_; }

    /// Set bounds for the parametrization of the Line.
    /// \param startpar start parameter
    /// \param endpar end parameter
    virtual void setParamBounds(double startpar, double endpar);

    // Translate the curve along a given vector
    virtual void translateCurve(const Point& dir);

    /// Query if parametrization is bounded. Both upper and lower
    /// parameter bounds must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    virtual bool isBounded() const;

    /// Confirm that the curve is linear
    virtual bool isLinear(Point& dir, double tol);

    /// If the curve is 2 dimensional, x and y coordinates will be swapped.
    /// Used when curve is a parameter curve.
    virtual void swapParameters2D();

    /// Check if the line lies in a plane passing through a given axis
    virtual bool isInPlane(const Point& loc, const Point& axis,
			   double eps, Point& normal) const;

    /// Check if the line lies in a plane with a given normal
    virtual bool isInPlane(const Point& norm,
			   double eps, Point& pos) const;

protected:

    Point location_;
    Point dir_;

    double parbound1_;
    double parbound2_;
    double startparam_;
    double endparam_;

};


} // namespace Go


#endif // _LINE_H

