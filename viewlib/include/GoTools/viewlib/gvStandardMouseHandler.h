//Added by qt3to4:
#include <QEvent>
//===========================================================================
//                                                                           
// File: gvStandardMouseHandler.h                                           
//                                                                           
// Created: Mon Apr 30 13:30:55 2001                                         
//                                                                           
// Author: Jens Olav Nygaard <jnygaard@sintef.math.no>
//                                                                           
// Revision: $Id: gvStandardMouseHandler.h,v 1.1 2007-04-17 12:25:44 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================






#ifndef GV_STANDARD_MOUSEHANDLER_H_INCLUDED





#if 0

// fra 'Attic':

#include "GoTools/viewlib/gvEventHandler.h"
class gvView;

class gvStandardMouseHandler : public gvEventHandler
{
public:
    gvStandardMouseHandler(gvView* view_to_control);
    virtual void handleEvent(QEvent* e);
};

#endif

/** Documentation ...
    etc
 */

class gvStandardMouseHandler
{
public:
    gvStandardMouseHandler();
    virtual ~gvStandardMouseHandler();

    // Insert other members...

};











#define GV_STANDARD_MOUSEHANDLER_H_INCLUDED

#endif // end of #ifndef GV_STANDARD_MOUSEHANDLER_H_INCLUDED
