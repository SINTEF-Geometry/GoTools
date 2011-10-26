//===========================================================================
//                                                                           
// File: gvPropertySheet.h                                                   
//                                                                           
// Created: Thu Nov 29 13:50:30 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvPropertySheet.h,v 1.1 2007-04-17 12:25:41 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVPROPERTYSHEET_H
#define _GVPROPERTYSHEET_H

class QWidget;
class gvObserver;

/** Documentation ...
    etc
 */

class gvPropertySheet
{
public:
    virtual ~gvPropertySheet();
    virtual void createSheet(QWidget* parent, gvObserver* observer) = 0;
};

#endif // _GVPROPERTYSHEET_H

