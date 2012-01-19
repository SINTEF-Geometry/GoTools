//===========================================================================
//                                                                           
// File: Tesselator.h                                                      
//                                                                           
// Created: Fri Apr 27 14:37:14 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: Tesselator.h,v 1.3 2009-05-13 07:32:51 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#ifndef TESSELATOR_H
#define TESSELATOR_H

#include "GoTools/utils/config.h"

namespace Go
{

/** Tesselator: super class for tesselators.
*/

class GO_API Tesselator
{
public:
  /// Destructor
    virtual ~Tesselator();
    /// Perform tesselation
    virtual void tesselate() = 0;
};

} // namespace Go




#endif // TESSELATOR_H

