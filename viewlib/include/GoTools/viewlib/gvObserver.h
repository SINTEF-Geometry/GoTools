//===========================================================================
//                                                                           
// File: gvObserver.h                                                        
//                                                                           
// Created: Tue Jul  3 12:56:05 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvObserver.h,v 1.1 2007-04-17 12:25:39 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVOBSERVER_H
#define _GVOBSERVER_H

/** Documentation ...
    etc
 */

class gvObserver
{
public:
    virtual ~gvObserver();
    virtual void observedChanged() = 0;
    
};


#endif // _GVOBSERVER_H

