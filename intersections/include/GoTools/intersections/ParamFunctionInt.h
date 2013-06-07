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

#ifndef _PARAMFUNCTIONINT_H
#define _PARAMFUNCTIONINT_H


#include "GoTools/intersections/ParamObjectInt.h"
#include "GoTools/intersections/BoundaryFunctionInt.h"
#include <memory>


namespace Go {


class ParamFunctionInt;
class Param0FunctionInt;
class Param1FunctionInt;
class Param2FunctionInt;


/// This class is a base class providing an interface to the
/// parametric "intersection objects" with 1-dimensional range.

class ParamFunctionInt : public ParamObjectInt {
public:
    /// Destructor.
    virtual ~ParamFunctionInt(){};

    /// If the object is not of type Param0FunctionInt, return NULL,
    /// otherwise return the object.
    virtual Param0FunctionInt* getParam0FunctionInt()
    { return NULL; } // To be overloaded in Param0FuntionInt.

    /// If the object is not of type Param1FunctionInt, return NULL,
    /// otherwise return the object.
    virtual Param1FunctionInt* getParam1FunctionInt()
    { return NULL; } // To be overloaded in Param1FuntionInt.

    /// If the object is not of type Param2FunctionInt, return NULL,
    /// otherwise return the object.
    virtual Param2FunctionInt* getParam2FunctionInt()
    { return NULL; } // To be overloaded in Param2FuntionInt.

    /// Subdivide the object in the specified parameter direction and
    /// parameter value.
    /// \param pardir direction in which to subdive. Indexing starts
    /// at 0.
    /// \param par parameter in which to subdivide.
    /// \param subdiv_objs The subparts of this object. Of the same
    /// geometric dimension as this object.
    /// \param bd_objs the boundaries between the returned \a
    /// subdiv_objs. Of geometric dimension 1 less than this object.
    virtual void subdivide(int pardir, double par, 
			   std::vector<shared_ptr<ParamFunctionInt> >& subdiv_objs,
			   std::vector<shared_ptr<ParamFunctionInt> >& bd_objs) = 0;

    /// Returns true if the object is monotone in any direction.
    /// \param dir the direction in which the object is
    /// monotone. Relevant only if object is monotone.

    // The returned dir is of interest for the 2par-case only.
    // Accepts non-strict monotonicity, as long as sf is not totally
    // flat.
    virtual bool monotone(Point& dir, double tol=1.0e-15) const = 0;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<shared_ptr<BoundaryFunctionInt> >& bd_objs) = 0;

    /// Return the size of the geometric sample mesh in the specified
    /// direction.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    virtual int getMeshSize(int dir) = 0;

    /// Return the corresponding mesh parameter.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param idx the mesh idx in the specified direction. Indexing
    /// starts at 0.
    virtual double paramFromMesh(int dir, int idx) = 0;

    /// Return the geometric sample mesh for the parametric function.
    virtual std::vector<double>::iterator getMesh() = 0;

protected:
    std::vector<shared_ptr<BoundaryFunctionInt> > boundary_obj_;

};

    
} // namespace Go


#endif // _PARAMFUNCTIONINT_H

