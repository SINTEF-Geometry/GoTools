//===========================================================================
//                                                                           
// File: CurveResolutionSheet.h                                          
//                                                                           
// Created: Thu Jan 31 13:05:30 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CurveResolutionSheet.h,v 1.2 2007-05-02 14:39:23 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CURVERESOLUTIONSHEET_H
#define _CURVERESOLUTIONSHEET_H

#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/viewlib/ui_CurveResolutionSheet_form.h"

#include <QtCore/qobject.h>


/** Ui_CurveResolutionSheet:
 */

class CurveResolutionSheet : public QObject //, public gvPropertySheet
{

Q_OBJECT

public:
    CurveResolutionSheet(int res = 500)
	: form_(0), obs_(0), res_(res), def_high_res_(5000)
    {}
    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void ok();

signals:
 void return_value(int); // Return new resolution.


private:
    Ui::CurveResolutionSheet_form* form_;
    gvObserver* obs_;
    int res_;
    int def_high_res_; // Used in setHighRes().

private slots:

    // Set default high resolution (currently 5000, optimized mode an advantage).
    void setHighRes();

};

#endif // _SPLINECURVEPROPERTYSHEET_H

