//==========================================================================
//                                                                          
// File: BoundaryFunctionInt.h                                               
//                                                                          
// Created: Mon Mar 13 13:42:40 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: BoundaryFunctionInt.h,v 1.1 2006-03-13 12:46:05 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _BOUNDARYFUNCTIONINT_H
#define _BOUNDARYFUNCTIONINT_H


#include "GoTools/utils/config.h"
#include <memory>


namespace Go {


class ParamFunctionInt;


/// This struct is a helper struct that bundles boundary information
/// for an object of type ParamFunctionInt.

struct BoundaryFunctionInt {

    shared_ptr<ParamFunctionInt> bd_obj_;
    int pardir_; // 2dim case: 0 || 1 (dir of par_, i.e. opp that of
		 // bd_obj_).
    double par_;

    /// Constructor
    /// \param bd_obj shared pointer to the object of interest
    /// \param dir parameter direction that specifies the boundary
    /// \param par the value of the relevant parameter at the boundary
    BoundaryFunctionInt(shared_ptr<ParamFunctionInt> bd_obj,
			int dir, double par)
    {
	bd_obj_ =  bd_obj;
	pardir_ = dir;
	par_ = par;
    }

    /// Get the parameter direction that specifies the boundary
    /// \return the parameter direction
    int getDir()
    { return pardir_; }

    /// Get the value of the relevant parameter at the boundary
    /// \return the parameter value
    double getPar()
    { return par_; }

    /// Get the object which the boundary belongs to
    /// \return shared pointer to the object
    shared_ptr<ParamFunctionInt> getObject()
    { return bd_obj_; }
};


} // namespace Go


#endif // _BOUNDARYFUNCTIONINT_H

