//===========================================================================
//                                                                           
// File: IntersectorFuncConst.h                                                
//                                                                           
// Created: Tue Sep 21 16:19:09 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: IntersectorFuncConst.h,v 1.14 2006-03-10 10:11:51 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

