//===========================================================================
//                                                                           
// File: DefaultDataHandler.h                                                
//                                                                           
// Created: Fri Jan  4 14:11:53 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: DefaultDataHandler.h,v 1.1 2007-04-17 12:25:30 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _DEFAULTDATAHANDLER_H
#define _DEFAULTDATAHANDLER_H

#include "GoTools/viewlib/DataHandler.h"

/** DefaultDataHandler: 
    etc
 */

class DefaultDataHandler : public DataHandler
{
public:
    /// Default constructor. Initializes factory (by registering classes).
    DefaultDataHandler();
    virtual ~DefaultDataHandler();

    virtual void create(shared_ptr<Go::GeomObject> obj,
			const gvColor& col, int id);
};

#endif // _DEFAULTDATAHANDLER_H

