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

#ifndef _PARAM0FUNCTIONINT_H
#define _PARAM0FUNCTIONINT_H


#include "GoTools/intersections/ParamFunctionInt.h"


namespace Go {


/// This is a class that represents the "intersection object" of a
/// scalar (a "0-variate" function).

class Param0FunctionInt : public ParamFunctionInt {
public:
    /// Constructor.
    /// \param C the constant defining the object.
    explicit Param0FunctionInt(double C);

    /// Constructor.
    /// \param C the constant defining the object.
    /// \param parent the parent object to this object.
    explicit Param0FunctionInt(double C,
			       ParamFunctionInt *parent);

    /// Destructor.
    virtual ~Param0FunctionInt();

    /// Evaluate the object in the input parameter.
    /// \param res the Point to be returned.
    /// \param par the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    virtual void point(Point& res, const double *par) const;

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param pt the Point to be returned.
    /// \param tpar the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    /// \param derivs the number of derivatives to calculate.
    /// \param from_right if true the evaluation is to be performed
    /// from the right side of the parameter value.
    /// \param resolution tolerance used when determining whether
    /// parameters are located at special values of the parameter
    /// domain (in particualar: knot values in case of spline objects.
    virtual void point(std::vector<Point>& pt, 
		       const double* tpar, 
		       int derivs,  
		       const bool* from_right = 0,
		       double resolution = 1.0e-12) const
    {
	DEBUG_ERROR_IF(derivs != 0,
		       "Asked for derivatives of a standalone point.");
	pt.resize(1);
	point(pt[0], tpar);
    }

    /// Return pointer to this object.
    virtual Param0FunctionInt* getParam0FunctionInt();

    /// The number of parameters in the object.
    virtual int numParams() const;

    /// Return an estimate on the size and wiggle of the object.
    /// \param length the approximative length of the object in the
    /// corresponding parameter directions.  The size of the array
    /// should be equal to numParams().
    /// \param wiggle a scalar representing the wiggle of the object
    /// in the corresponding parameter directions. The size of the
    /// array should be equal to numParams().
    virtual void getLengthAndWiggle(double *length, double *wiggle);

    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const;

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir,
						 bool sort = false) const;

    /// Return true if the object has any critical parameter values in
    /// the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalVals(int pardir) const;

    /// Return the critical parameter values in the specified
    /// direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalVals(int pardir) const;

    /// Return true if the object has any critical parameter values or
    /// inner knots in the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const;

    /// Return the critical parameter values and inner knots for
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const;

    /// Return true if we are allowed to divide in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question.
    virtual bool canDivide(int pardir);

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

    /// Return the geometric sample mesh for the parametric function.
    virtual std::vector<double>::iterator getMesh();

    /// Return the start parameter value in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double startParam(int pardir) const;

    /// Return the end parameter in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double endParam(int pardir) const;

    /// Return true if the specified point lies within eps from the
    /// boundary.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    /// \param eps the tolerance defining boundary neighbourhood.
    virtual bool boundaryPoint(const double* par, double eps) const;

    /// Subdivide the object in the specified parameter direction and
    /// parameter value.
    /// \param pardir direction in which to subdive. Indexing starts
    /// at 0.
    /// \param par parameter in which to subdivide.
    /// \param subdiv_objs the subparts of this object. Of the same
    /// geometric dimension as this object.
    /// \param bd_objs the boundaries between the returned \a
    /// subdiv_objs. Of geometric dimension 1 less than this object.
    virtual void subdivide(int pardir, double par, 
			   std::vector<shared_ptr<ParamFunctionInt> >& subdiv_objs,
			   std::vector<shared_ptr<ParamFunctionInt> >& bd_objs);

    /// Return the CompositeBox for the parametric object.
    virtual CompositeBox compositeBox() const;

    /// Return true if the object is monotone in any direction.
    /// \param dir the direction in which the object is
    /// monotone. Relevant only if object is monotone.

    // The returned dir is of interest for the 2par-case only.
    // Accepts non-strict monotonicity, as long as sf is not totally flat.
    virtual bool monotone(Point& dir, double tol=1.0e-15) const; // = 0;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void
    getBoundaryObjects(std::vector<shared_ptr<BoundaryFunctionInt> >& bd_objs);

    /// Return the dimension of the geometric space.
    int dimension();

    // Functions used from coincidence checking
 
    // Find the interval in which t lies, if t is within tol of an
    // existing knot, t will be changed to that knot.

    /// Return the knot interval for which t lies inside, possibly
    /// moving the value t if it lies close to a knot.
    /// \param t the parameter value
    /// \param tol the parametric tolerance deciding if the input
    /// parameter should be moved.
    virtual int knotIntervalFuzzy(double& t, double tol) const;
    
    /// Return the value of the knot closest to the input parameter.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    virtual double nextSegmentVal(double par, bool forward) const;

    /// Return the scalar value of this object.
    double getValue() const
    { return C_; }

protected:

    double C_; // The 0-parameter "function".

    int dim_; // Space dimension. 1.

    ParamFunctionInt* parentcurve_; // @@sbr Or surface ...

};


} // namespace Go


#endif // _PARAM0FUNCTIONINT_H

