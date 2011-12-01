//===========================================================================
//                                                                           
// File: SplineCurvePropertySheet.h                                          
//                                                                           
// Created: Thu Jan 31 13:05:30 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SplineCurvePropertySheet.h,v 1.2 2007-05-02 14:39:24 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SPLINECURVEPROPERTYSHEET_H
#define _SPLINECURVEPROPERTYSHEET_H

#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/viewlib/ui_SplineCurvePropertySheet_form.h"
#include "GoTools/tesselator/CurveTesselator.h"

#include <QObject>

class gvData;
//class CurveTesselator;
class gvCurvePaintable;
// class SplineCurvePropertySheet_form;

/** Documentation ...
    etc
 */
class SplineCurvePropertySheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
  SplineCurvePropertySheet(Go::CurveTesselator* tess, gvCurvePaintable* pable)
	: tess_(tess), pable_(pable), form_(0), obs_(0)
    {}
    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void apply();
 

private:
    Go::CurveTesselator* tess_;
    gvCurvePaintable* pable_;
    Ui::SplineCurvePropertySheet_form* form_;
    gvObserver* obs_;
};

#endif // _SPLINECURVEPROPERTYSHEET_H

