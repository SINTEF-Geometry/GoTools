//===========================================================================
//                                                                           
// File: ParamFunctionInt.h                                                  
//                                                                           
// Created: Wed Sep 22 12:35:33 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: ParamFunctionInt.h,v 1.16 2006-09-01 12:25:26 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

