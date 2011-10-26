//===========================================================================
//                                                                           
// File: NoopTesselator.h                                                  
//                                                                           
// Created: Tue Feb 19 16:37:59 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvNoopTesselator.h,v 1.1 2007-04-17 12:25:38 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _NOOPTESSELATOR_H
#define _NOOPTESSELATOR_H

#include "GoTools/tesselator/Tesselator.h"

namespace Go
{

/** Documentation ...
    etc
 */
class NoopTesselator : public Tesselator
{
public:
    virtual ~NoopTesselator();
    virtual void tesselate();
};


} // namespace Go




#endif // _NOOPTESSELATOR_H

