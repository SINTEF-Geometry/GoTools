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

#ifndef _PAR1FUNCINTERSECTOR_H
#define _PAR1FUNCINTERSECTOR_H


#include "GoTools/intersections/IntersectorFuncConst.h"


namespace Go {


/// This class is performing intersections between a 1-dimensional
/// parametric curve and a constant.

class Par1FuncIntersector : public IntersectorFuncConst {
public:

//     Par1FuncIntersector(shared_ptr<ParamObjectInt> func,
// 			shared_ptr<ParamObjectInt> C,
// 			double epsge,
// 			Intersector* prev = 0,
// 			int eliminated_parameter = -1,
// 			double eliminated_value = 0);			

    /// Constructor.
    /// One of the objects should refer to a 1D-curve, the other a
    /// constant (this is not checked compile-time, so we rely on the
    /// user to obey this rule).  The last two variables are relevant
    /// only if the parent has one more parameter than the Intersector
    /// to be constructed.
    /// \param func of type Param1FunctionInt.
    /// \param C of type Param0FunctionInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index (0 or 1) of the
    /// parameter that was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    Par1FuncIntersector(shared_ptr<ParamFunctionInt> func,
			shared_ptr<ParamFunctionInt> C,
			shared_ptr<GeoTol> epsge,
			Intersector *prev = 0,
			int eliminated_parameter = -1,
			double eliminated_value = 0);

    /// Destructor
    virtual ~Par1FuncIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 1; }
	
protected:
    // Data members

    virtual shared_ptr<Intersector> 
    lowerOrderIntersector(shared_ptr<ParamFunctionInt> obj1,
			  shared_ptr<ParamFunctionInt> obj2,
			  Intersector* prev = 0,
			  int eliminated_parameter = -1,
			  double eliminated_value = 0);

    virtual int checkCoincidence();

    virtual void microCase();
    
    virtual int updateIntersections();

    virtual int repairIntersections()
    { return 0; }

    virtual int doSubdivide();

private:

    // We return par in obj1_ of closest pt between objects.
    // The corr dist and valid range are also returned.
    void doIterate(double& par1, double& dist, double& tmin, double& tmax);

    int getSubdivisionParameter(int dir, double& par);

    void getIterationSeed(double seed[]);

};


} // namespace Go


#endif // _PAR1FUNCINTERSECTOR_H

