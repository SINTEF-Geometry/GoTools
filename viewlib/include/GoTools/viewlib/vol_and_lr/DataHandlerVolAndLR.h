//===========================================================================
//                                                                           
// File: DataHandlerVolAndLR.h                                               
//                                                                           
// Created: Fri Feb  8 16:21:46 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _DATAHANDLERVOLANDLR_H
#define _DATAHANDLERVOLANDLR_H

#include "GoTools/viewlib/DefaultDataHandler.h"

class DataHandlerVolAndLR : public DefaultDataHandler
{
public:
    /// Default constructor. Initializes factory (by registering classes).
    DataHandlerVolAndLR();
    virtual ~DataHandlerVolAndLR();

    virtual void create(shared_ptr<Go::GeomObject> obj,
			const gvColor& col, int id);
};

#endif // _DATAHANDLERVOLANDLR_H

