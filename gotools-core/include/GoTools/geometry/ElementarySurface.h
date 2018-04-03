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

#ifndef _ELEMENTARYSURFACE_H
#define _ELEMENTARYSURFACE_H



#include "GoTools/geometry/ParamSurface.h"


namespace Go
{


class SplineSurface;
class ElementaryCurve;


/// \brief ElementarySurface is a base class for elementary surfaces
/// like planes and cylinders. Such surfaces have natural
/// parametrizations and ElementartSurface is therefore a subclass of
/// ParamSurface. These surfaces are non-self-intersecting.
class ElementarySurface : public ParamSurface
{
public:
    ElementarySurface();

    /// Virtual destructor, enables safe inheritance.
    virtual ~ElementarySurface();

    // --- Functions inherited from GeomObject ---
    virtual ElementarySurface* clone() const = 0;

    // --- Functions inherited from ParamSurface ---

    virtual RectDomain containingDomain() const;

    virtual CurveLoop outerBoundaryLoop(double degenerate_epsilon
					  = DEFAULT_SPACE_EPSILON) const;

    virtual std::vector<CurveLoop> allBoundaryLoops(double degenerate_epsilon
						      = DEFAULT_SPACE_EPSILON) const;

    virtual Point closestInDomain(double u, double v) const;

    /// Check if a parameter pair lies inside the domain of this surface
    virtual bool inDomain(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies inside the domain of this surface
    /// return value = 0: outside
    ///              = 1: internal
    ///              = 2: at the boundary
    virtual int inDomain2(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies at the boundary of this surface
    virtual bool onBoundary(double u, double v, double eps=1.0e-4) const;

    virtual double area(double tol) const;

    /// Return surface corners, geometric and parametric points
    /// in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    virtual SplineSurface* asSplineSurface();

     /// When the surface is divided up into logical segments, this function will return 
    /// the parameter value of the "next segment", starting from a parameter given by the user.
    /// In this case no logical segments exist, then it is the start- or end parameter that
    /// is returned.   
    virtual double nextSegmentVal(int dir, double par, bool forward, 
				  double tol) const;

    virtual bool isDegenerate(bool& b, bool& r,
			      bool& t, bool& l, double tolerance) const = 0;

    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const = 0;

    /// Query if parametrization is bounded. All four parameter bounds
    /// must be finite for this to be true.
    /// \return \a true if bounded, \a false otherwise
    virtual bool isBounded() const = 0;

    // // --- Functions specific to ElementarySurface ---

    virtual bool isClosed(bool& closed_dir_u, bool& closed_dir_v) const;

    virtual SplineSurface* geometrySurface() const = 0;
    /// Create a SplineSurface representation of the elementary surface
    virtual SplineSurface* createSplineSurface() const = 0;

    /// Limit the surface by limiting the parameter domain
    virtual void setParameterBounds(double from_upar, double from_vpar,
				    double to_upar, double to_vpar) = 0;

    virtual void turnOrientation();
    virtual void reverseParameterDirection(bool direction_is_u);
    virtual void swapParameterDirection();

    //virtual bool isReversedU() const;
    //virtual bool isReversedV() const;
    virtual bool isSwapped() const;

    /// Return associated elementary surface, if any
    virtual ElementarySurface* elementarySurface()
    {
      return this;
    }

    // ---  Functions in ElementarySurface  ---
    /// Fetch the parameter curve in the domain of the elementary surface
    /// corresponding to a given elementary curve in geometry space
    /// if this curve has a simpler elementary representation.
    /// Otherwise, nothing is returned
    virtual shared_ptr<ElementaryCurve> 
      getElementaryParamCurve(ElementaryCurve* space_crv, double tol,
			      const Point* start_par_pt = NULL, const Point* end_par_pt = NULL) const;

    /// Radius in a specified location, default -1
    virtual double radius(double u, double v) const
    {
      return -1.0;
    }

    virtual Point location() const
    {
      Point dummy;
      return dummy;
    }

    virtual Point direction() const
    {
      Point dummy;
      return dummy;
    }

    virtual Point direction2() const
    {
      Point dummy;
      return dummy;
    }

    virtual void enlarge(double len1, double len2, double len3, double len4) = 0;

protected:
    double ptol_;  // Tolerance used in decisions on parameter range
    // when no nother tolerance information is available

    //bool isReversedU_;
    //bool isReversedV_;
    bool isSwapped_;

    // Helper function to be used in functions like Cylinder::point(), where
    // we need to take the isSwapped_ flag into account.
    void getOrientedParameters(double& u, double&v) const;

};


} // namespace Go



#endif // _ELEMENTARYSURFACE_H

