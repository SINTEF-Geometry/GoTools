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

#ifndef _INTERSECTORFUNCOBJ_H
#define _INTERSECTORFUNCOBJ_H


#include "GoTools/intersections/Intersector.h"


namespace Go {


class ParamFunctionInt;


/// This class is an interface class for finding the intersection
/// between a parametric function (with range R) and a constant.

class IntersectorFuncConst : public Intersector {
public:

    /// Constructor.
    /// The last two variables are relevant only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param func of type ParamFunctionInt.
    /// \param C of type Param0FunctionInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index of the parameter that
    /// was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    IntersectorFuncConst(shared_ptr<ParamFunctionInt> func,
			 shared_ptr<ParamFunctionInt> C,
			 shared_ptr<GeoTol> epsge,
			 Intersector* prev = 0,
			 int eliminated_parameter = -1,
			 double eliminated_value = 0);

    /// Destructor
    virtual ~IntersectorFuncConst();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);  @bsp

    friend class IntersectorAlgPar;

protected:

    shared_ptr<ParamFunctionInt> func_int_;
    shared_ptr<ParamFunctionInt> C_;

    // @@sbr Currently we need intersection between a
    // Param1FunctionInt & ParamPointInt.  Possibly introduce
    // additional object (replacing ParamPointInt) in this branch of
    // the tree.  NB: The order of the objects ot input is not
    // arbitrary!  The knowledge of what is the 'first object' and the
    // 'second object' can be used internally, and must be consistent
    // with the parent Intersector.  The second object could be a
    // ParamPointInt or Param1FunctionInt.
    virtual void print_objs();

    virtual shared_ptr<Intersector>
    lowerOrderIntersector(shared_ptr<ParamFunctionInt> obj1,
			  shared_ptr<ParamFunctionInt> obj2,
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0) = 0;

    virtual int getBoundaryIntersections();

    virtual int performInterception();

    virtual int simpleCase();

    virtual bool isLinear();

    virtual bool complexityReduced()
    {
	// Default behaviour, continue recursion
	return true;
    }

    virtual void handleComplexity()
    {
	// Default, do nothing. Must be implemented together with
	// complexityReduced
    }

    virtual int checkCoincidence() = 0;

    virtual void microCase() = 0;
    
    virtual int updateIntersections() = 0;

    virtual int linearCase();

    virtual int doSubdivide() = 0;

    virtual void printDebugInfo();

private:

};


} // namespace Go


#endif // _INTERSECTORFUNCOBJ_H

