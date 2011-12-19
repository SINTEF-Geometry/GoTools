//===========================================================================
//                                                                           
// File: Par0FuncIntersector.h                                               
//                                                                           
// Created: Fri Oct  1 10:35:00 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Par0FuncIntersector.h,v 1.8 2006-11-03 14:43:14 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PAR0FUNCINTERSECTOR_H
#define _PAR0FUNCINTERSECTOR_H


#include "GoTools/intersections/IntersectorFuncConst.h"


namespace Go {


/// This class is performing intersections between two constants.

class Par0FuncIntersector : public IntersectorFuncConst {
public:

//     Par0FuncIntersector(shared_ptr<ParamObjectInt> func,
// 			shared_ptr<ParamObjectInt> C,
// 			double epsge,
// 			Intersector* prev = 0,
// 			int eliminated_parameter = -1,
// 			double eliminated_value = 0);			

    /// Constructor.
    /// Both objects should refer to constants (this is not checked
    /// compile-time, so we rely on the user to obey this rule).  The
    /// last two variables are relevant only if the parent has one
    /// more parameter than the Intersector to be constructed.
    /// \param func of type Param0FunctionInt.
    /// \param C of type Param0FunctionInt.
    /// \param epsge the associated tolerance.
    /// \param prev the "parent" Intersector (0 if there is no
    /// parent).
    /// \param eliminated_parameter the index (0) of the parameter
    /// that was removed from the parent \a prev.
    /// \param eliminated_value the value of the parameter that was
    /// removed from the parent \a prev.
    Par0FuncIntersector(shared_ptr<ParamFunctionInt> func,
			shared_ptr<ParamFunctionInt> C,
			shared_ptr<GeoTol> epsge,
			Intersector *prev = 0,
			int eliminated_parameter = -1,
			double eliminated_value = 0);

    /// Destructor.
    virtual ~Par0FuncIntersector();

//     // Validation of given intersection results
//     virtual void validate(int level, ValidationStat status);

    /// Return the number of parameter directions for the object.
    /// \return the number of parameter directions
    virtual int numParams() const
    { return 0; }
	
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

};


} // namespace Go


#endif // _PAR0FUNCINTERSECTOR_H

