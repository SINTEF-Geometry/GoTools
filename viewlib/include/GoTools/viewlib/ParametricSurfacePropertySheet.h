//===========================================================================
//                                                                           
// File: ParametricSurfacePropertySheet.h                                        
//                                                                           
// Created:
//                                                                           
// Author: Sverre Briseid
//                                                                           
// Revision: $Id: ParametricSurfacePropertySheet.h,v 1.2 2007-05-02 14:39:23 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARAMETRICSURFACEPROPERTYSHEET_H
#define _PARAMETRICSURFACEPROPERTYSHEET_H


#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/viewlib/ui_RectangularSurfacePropertySheet_form.h"
#include "GoTools/geometry/ParamSurface.h"
//#include "GoTools/viewlib/gvParametricSurfaceTesselator2.h"
#include "GoTools/viewlib/gvParametricSurfacePaintable.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"

#include <QObject>

class gvData;
//class ParametricSurfaceTesselator;
//class gvParametricSurfacePaintable;

/** Documentation ...
    etc
 */

class ParametricSurfacePropertySheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
    ParametricSurfacePropertySheet()
    {}

    ParametricSurfacePropertySheet(Go::ParametricSurfaceTesselator* tess,
				   gvParametricSurfacePaintable* pable,
				   std::shared_ptr<Go::ParamSurface>& surf)
	: tess_(tess), pable_(pable), form_(0), obs_(0), surf_(surf)
    {}

    virtual ~ParametricSurfacePropertySheet();

    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void apply();
 

private:
  Go::ParametricSurfaceTesselator* tess_;
  gvParametricSurfacePaintable* pable_;
  Ui::RectangularSurfacePropertySheet_form* form_;
  gvObserver* obs_;
  std::shared_ptr<Go::ParamSurface> surf_;
};


#endif // _PARAMETRICSURFACEPROPERTYSHEET_H

