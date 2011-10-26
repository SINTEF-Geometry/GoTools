//===========================================================================
//                                                                           
// File: SurfaceResolutionSheet.h                                        
//                                                                           
// Created: Mon Jan  7 13:26:03 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SurfaceResolutionSheet.h,v 1.2 2007-05-02 14:39:24 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SURFACERESOLUTIONSHEET_H
#define _SURFACERESOLUTIONSHEET_H


#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "ui_SurfaceResolutionSheet_form.h"

#include <QObject>

class gvData;
//class RectangularSurfaceTesselator;
class gvRectangularSurfacePaintable;

/** Documentation ...
    etc
 */

class SurfaceResolutionSheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
SurfaceResolutionSheet(int ures = 20, int vres = 20) // The default res is 20x20.
     : form_(0), ures_(ures), vres_(vres), def_high_res_(200) //, obs_(0)
    {}
    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void ok();

signals:
    void return_value(int, int); // Return new resolution.


private:
    Ui::SurfaceResolutionSheet_form* form_;
    int ures_;
    int vres_;
    int def_high_res_; // Used in setHighRes().
/*     gvObserver* obs_; */

private slots:

    // Set default high resolution (optimized mode an advantage).
    // Accessed through button.
    void setHighRes();

};


#endif // _SPLINESURFACEPROPERTYSHEET_H

