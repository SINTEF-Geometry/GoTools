//===========================================================================
//                                                                           
// File: RectangularSurfacePropertySheet.h                                        
//                                                                           
// Created: Mon Jan  7 13:26:03 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectangularSurfacePropertySheet.h,v 1.2 2007-05-02 14:39:23 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _RECTANGULARSURFACEPROPERTYSHEET_H
#define _RECTANGULARSURFACEPROPERTYSHEET_H


#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/viewlib/ui_RectangularSurfacePropertySheet_form.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"

#include <QObject>

class gvData;
//class RectangularSurfaceTesselator;
class gvRectangularSurfacePaintable;
// class RectangularSurfacePropertySheet_form;
//class ParamSurface;

/** Documentation ...
    etc
 */

class RectangularSurfacePropertySheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
    RectangularSurfacePropertySheet(Go::RectangularSurfaceTesselator* tess,
				    gvRectangularSurfacePaintable* pable,
				    shared_ptr<Go::ParamSurface>& surf)
  : tess_(tess), pable_(pable), form_(), obs_(0), surf_(surf)
    {}

    virtual ~RectangularSurfacePropertySheet();

    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void apply();
 

private:
    Go::RectangularSurfaceTesselator* tess_;
    gvRectangularSurfacePaintable* pable_;
    Ui::RectangularSurfacePropertySheet_form* form_;
    gvObserver* obs_;
    shared_ptr<Go::ParamSurface> surf_;
};


#endif // _RECTANGULARSURFACEPROPERTYSHEET_H

