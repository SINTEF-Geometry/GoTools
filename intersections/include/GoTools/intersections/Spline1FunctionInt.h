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

#ifndef _SPLINE1FUNCTIONINT_H
#define _SPLINE1FUNCTIONINT_H


#include "GoTools/intersections/Param1FunctionInt.h"


namespace Go {


class SplineCurve;


/// Class that represents the "intersection object" of a spline curve
/// of dimension 1.

class Spline1FunctionInt : public Param1FunctionInt {
public:
    /// Constructor.
    /// Input curve must be of type SplineCurve. This is not checked
    /// run-time, so we rely on the user to obey this rule.
    /// \param curve the parametric 1-dimensional curve defining the
    /// object.
    explicit Spline1FunctionInt(shared_ptr<ParamCurve> curve);

    /// Constructor.
    /// Input curve must be of type SplineCurve. This is not checked
    /// run-time, so we rely on the user to obey this rule.
    /// \param curve the parametric 1-dimensional curve defining the
    /// object. Can be either a curve or a surface.
    /// \param parent the parent object to this object.
    explicit Spline1FunctionInt(shared_ptr<ParamCurve> curve, 
				ParamFunctionInt *parent);

    /// Destructor.
    virtual ~Spline1FunctionInt() {};

    /// Return an intersection object for the input curve.
    /// Input curve must be of type SplineCurve. This is not checked
    /// run-time, so we rely on the user to obey this rule.  This
    /// object is used as the parent for the intersection object.
    /// \param curve the parametric curve defining the intersection
    /// object.
    virtual shared_ptr<Param1FunctionInt> 
    makeIntFunction(shared_ptr<ParamCurve> curve);
    
    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const;

    /// Return true if the object has any critical parameter values or
    /// inner knots.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const;

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir,
						 bool sort = false) const;

    /// Return the critical parameter values and inner knots for
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const; 

    /// Return the size of the geometric sample mesh in the specified
    /// direction.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    virtual int getMeshSize(int dir);

    /// Return the corresponding mesh parameter.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param idx the mesh idx in the specified direction. Indexing
    /// starts at 0.
    virtual double paramFromMesh(int dir, int idx);

    /// Return the geometric sample mesh for the spline function.
    /// \return The geometric sample mesh for the spline function.
    virtual std::vector<double>::iterator getMesh();

    // Accepts non-strict monotonicity, as long as sf is not totally flat.
    /// Return true if the curve is monotone.
    /// \param dir the direction in which the object is monotone. Is
    /// not of interest here as the curve has only 1 parameter
    /// direction.
    /// \return Whether or not the curve is monotone.
    virtual bool monotone(Point& dir, double tol=1.0e-15) const;

    // Functions used for coincidence testint

    // Find the interval in which t lies, if t is within tol of an
    // existing knot, t will be changed to that knot.
    /// Return the knot interval for which t lies inside, moving the
    /// value t if it lies close to a knot.
    /// \param t the parameter value
    /// \param tol the parametric tolerance deciding if the input
    /// parameter t should be moved.
    virtual int knotIntervalFuzzy(double& t, double tol) const;

    /// Return the value of the knot next to the input parameter par.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(double par, bool forward) const;

protected:

    // Data members
    shared_ptr<SplineCurve> spcv_; // Same cv as in
					  // Param1FunctionInt, casted
					  // to spline.

};


} // namespace Go


#endif // _SPLINE1FUNCTIONINT_H

    
