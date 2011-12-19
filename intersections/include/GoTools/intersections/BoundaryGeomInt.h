//==========================================================================
//                                                                          
// File: BoundaryGeomInt.h                                                   
//                                                                          
// Created: Fri Mar 10 17:07:42 2006                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: BoundaryGeomInt.h,v 1.1 2006-03-13 12:37:58 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _BOUNDARYGEOMINT_H
#define _BOUNDARYGEOMINT_H


#include <memory>


namespace Go {


class ParamGeomInt;


/// This struct is a helper struct that bundles boundary information
/// for an object of type ParamGeomInt.

struct BoundaryGeomInt {

    shared_ptr<ParamGeomInt> bd_obj_;
    int pardir_; // 2dim case: 0 || 1 (dir of par_, i.e. opp that of
		 // bd_obj_).
    double par_;

    /// Constructor
    /// \param bd_obj shared pointer to the object of interest
    /// \param dir parameter direction that specifies the boundary
    /// \param par the value of the relevant parameter at the boundary
    BoundaryGeomInt(shared_ptr<ParamGeomInt> bd_obj,
		    int dir, double par)
    {
	bd_obj_ =  bd_obj;
	pardir_ = dir;
	par_ = par;
    }

    /// Destructor
    ~BoundaryGeomInt() { }

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
    shared_ptr<ParamGeomInt> getObject()
    { return bd_obj_; }
};


} // namespace Go


#endif // _BOUNDARYGEOMINT_H

