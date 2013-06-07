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

#ifndef _SPLINE2FUNCTIONINT_H
#define _SPLINE2FUNCTIONINT_H


#include "GoTools/intersections/Param2FunctionInt.h"


namespace Go {


 /// This class represents the "intersection object" of a spline
 /// surface of dimension 1.

class Spline2FunctionInt : public Param2FunctionInt {
public:
    /// Constructor.
    /// \param surf the parametric 1-dimensional surface defining the
    /// object.
    explicit Spline2FunctionInt(shared_ptr<SplineSurface> surf);

    /// Constructor.
    /// \param surf the parametric 1-dimensional surface defining the
    /// object.
    /// \param parent the parent object to this object.
    explicit Spline2FunctionInt(shared_ptr<SplineSurface> surf, 
				Param2FunctionInt *parent);

    /// Destructor.
    virtual ~Spline2FunctionInt(){};

    /// Return an intersection object for the input surface.  This
    /// object is used as the parent for the intersection object.
    /// \param surf the parametric surface defining the intersection
    /// object.
    virtual shared_ptr<Param2FunctionInt>
    makeIntFunction(shared_ptr<ParamSurface> surf);
    
    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const;

    /// Return true if the object has any critical parameter values or
    /// inner knots in the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const;

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir,
						 bool sort=false) const;

    /// Return the critical parameter values and inner knots for
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const;

//     void getLengthAndWiggle(double *length, double *wiggle);

    // pardir = 0 || 1
    // @@sbr A trimmed sf need not have a startparam/endparam ...
    /// Return the start parameter value in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \return The start parameter in the specified direction.
    virtual double startParam(int pardir) const;

    /// Return the end parameter in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \return The end parameter in the specified direction.
    virtual double endParam(int pardir) const;

    // We want to know if there exists par direction in which the
    // function is monotone (guarantees at most one solution for a
    // continuous function).  We also would like to know in what
    // direction the function is monotone (will then move in
    // orthogonal direction).
    /// Return true if the surface is monotone in any direction.
    /// \param dir the direction in which the surface is
    /// monotone. Pertains only if surface is monotone.
    /// \return Whether or not the surface is monotone.
    virtual bool monotone(Point& dir, double tol=1.0e-15) const;

    // Functions used for coincidence testing

    // Find the interval in which t lies, if t is within tol of an
    // existing knot, t will be changed to that knot.
    /// Return the knot interval for which t lies inside, moving the
    /// value t if it lies close to a knot.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param t the parameter value
    /// \param tol the parametric tolerance deciding if the input
    /// parameter t should be moved.
    virtual int knotIntervalFuzzy(int pardir, double& t, double tol) const;

    /// Return the value of the knot next to the input parameter par.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(int pardir, double par, bool forward) const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<shared_ptr<BoundaryFunctionInt> >& bd_objs);

    /// Return a 3-dimensional visualization spline surface.
    /// \return The pointer to the 3-dimensional spline surface (u, v,
    /// surf_(u, v)).
    shared_ptr<SplineSurface> surface3D();

    /// Return the gradient of the function S as a two-dimensional
    /// spline surface: (dS/du, dS/dv).
    /// \return Pointer to the gradient spline surface
    shared_ptr<SplineSurface> createGradSurface() const;

protected:
    // Data members
    shared_ptr<SplineSurface> spsf_; // shared_ptr to this
					    // curve

};


} // namespace Go


#endif // _SPLINE2FUNCTIONINT_H

