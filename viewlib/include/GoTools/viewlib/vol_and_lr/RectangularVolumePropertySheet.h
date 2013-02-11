//===========================================================================
//                                                                           
// File: RectangularVolumePropertySheet.h                                    
//                                                                           
// Created: Thu Jul  5 15:03:15 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _RECTANGULARVOLUMEPROPERTYSHEET_H
#define _RECTANGULARVOLUMEPROPERTYSHEET_H



#include "GoTools/viewlib/vol_and_lr/ui_RectangularVolumePropertySheet_form.h"
#include "GoTools/viewlib/vol_and_lr/gvRectangularVolumePaintable.h"

#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/RectangularVolumeTesselator.h"

#include <QObject>

class gvData;
//class RectangularSurfaceTesselator;
class gvRectangularVolumePaintable;
// class RectangularSurfacePropertySheet_form;
//class ParamSurface;

/** Documentation ...
    etc
 */

class RectangularVolumePropertySheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
    RectangularVolumePropertySheet(Go::RectangularVolumeTesselator* tess,
				   gvRectangularVolumePaintable* pable,
				   shared_ptr<Go::ParamVolume>& vol)
  : tess_(tess), pable_(pable), form_(), obs_(0), vol_(vol)
    {}

    virtual ~RectangularVolumePropertySheet();

    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void apply();
 

private:
    Go::RectangularVolumeTesselator* tess_;
    gvRectangularVolumePaintable* pable_;
    Ui::RectangularVolumePropertySheet_form* form_;
    gvObserver* obs_;
    shared_ptr<Go::ParamVolume> vol_;
};


#endif // _RECTANGULARVOLUMEPROPERTYSHEET_H

