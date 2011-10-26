//===========================================================================
//                                                                           
// File: ftPoint.C                                                           
//                                                                           
// Created: Wed May 10 11:30:51 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ftPoint.C,v 1.2 2008-04-02 07:18:29 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftSurface.h"

namespace Go
{

Point ftPoint::normal() const
{
    ALWAYS_ERROR_IF(surface_ == 0, "There is no underlying surface!");

    return surface_->normal(u_, v_);
}

} // namespace Go
