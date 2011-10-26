//===========================================================================
//                                                                           
// File: ttlPoint.h
//                                                                           
// Created: July 9 2001
//                                                                           
// Author: Sverre Briseid <sbr@math.sintef.no>
//                                                                           
// Revision: $Id: ttlPoint.h,v 1.8 2008-07-10 17:54:10 kfp Exp $
//                                                                           
// Description: Point suited for triangulation (setting up neighbourhood).
//                                                                           
//===========================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//===========================================================================
#ifndef _TTLPOINT_H_
#define _TTLPOINT_H_

#include "GoTools/compositemodel/ftPointSet.h"

namespace Go
{

  /** Utility class for representing a ftSamplePoint with additional parameter value.
   * Useful when dealing with sets of surfaces treated as one.
   */
  class ttlPoint
  {

  public:
    /// Constructor
  ttlPoint(PointIter pnt_iter, double x, double y, double z = 0.0)
    : par_val_(Vector3D(x, y, z)), pnt_iter_(pnt_iter) {}
    /// Destructor
    ~ttlPoint() {}
  
    inline const PointIter& pnt_iter() const { return pnt_iter_; };
    inline double x() const { return par_val_[0]; } ;
    inline double y() const { return par_val_[1]; } ;
    inline double z() const { return par_val_[2]; } ;


  private:

    // As point is to be included in a larger set, we must set new parameter value.
    Vector3D par_val_; // z coordinate is by default set to 0.
    PointIter pnt_iter_; // We need access to the actual object.

  };


} // namespace Go

#endif
