//===========================================================================
//                                                                           
// File: Par1FuncIntersector.h                                               
//                                                                           
// Created: Tue Sep 21 13:54:16 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Par1FuncIntersector.h,v 1.9 2006-11-03 14:43:14 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

