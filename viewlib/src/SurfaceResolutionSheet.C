//===========================================================================
//                                                                           
// File: SurfaceResolutionSheet.C                                        
//                                                                           
// Created: Wed Oct 23 14:49:10 2002                                         
//                                                                           
// Author: Sverre Briseid
//                                                                           
// Revision:
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/SurfaceResolutionSheet.h"
#include "GoTools/viewlib/ui_SurfaceResolutionSheet_form.h"

// #include <q3groupbox.h>
#include <QPushButton>
#include <QSlider>
#include <QCheckBox>

using namespace Ui;

//===========================================================================
void SurfaceResolutionSheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
    form_ = new Ui::SurfaceResolutionSheet_form();//parent);

    QDialog* w = new QDialog();
    form_->setupUi(w);

    w->resize(form_->box->size());

    connect(form_->OkButton, SIGNAL(clicked()),
	    this, SLOT(ok()));
    connect(form_->CancelButton, SIGNAL(clicked()),
	    w, SLOT(close()));
    connect(form_->button200x200, SIGNAL(clicked()),
	    this, SLOT(setHighRes()));

    form_->UresSlider->setValue(ures_);
    form_->VresSlider->setValue(vres_);

    w->show();
}


//===========================================================================
void SurfaceResolutionSheet::ok()
//===========================================================================
{
    ures_ = form_->UresSlider->value();
    vres_ = form_->VresSlider->value();
//     obs_->observedChanged();

    emit return_value(ures_, vres_);
//     form_->close();
//     this->close();
}


//=========================================================================== 
void SurfaceResolutionSheet::setHighRes()
//=========================================================================== 
{
    form_->UresSlider->setValue(def_high_res_);
    form_->VresSlider->setValue(def_high_res_);

    ok();
}
