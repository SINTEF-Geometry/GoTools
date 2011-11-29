//===========================================================================
//                                                                           
// File: CurveResolutionSheet.C                                          
//                                                                           
// Created: Wed Oct 23 14:49:10 2002                                         
//                                                                           
// Author: Sverre Briseid
//                                                                           
// Revision: $Id: CurveResolutionSheet.C,v 1.3 2007-05-02 14:39:25 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/CurveResolutionSheet.h"
#include "GoTools/viewlib/ui_CurveResolutionSheet_form.h"

//#include <q3groupbox.h>
#include <QPushButton>
#include <QSlider>
#include <QCheckBox>

using namespace Ui;

//===========================================================================
void CurveResolutionSheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
    obs_ = obs;
    form_ = new Ui::CurveResolutionSheet_form();//parent);

    QDialog* dial = new QDialog();
    form_->setupUi(dial);

    dial->resize(form_->box->size());

    QObject::connect(form_->OkButton, SIGNAL(clicked()),
		     this, SLOT(ok()));
    QObject::connect(form_->CancelButton, SIGNAL(clicked()),
		     dial, SLOT(close()));
    connect(form_->button5000, SIGNAL(clicked()),
	    this, SLOT(setHighRes()));

    form_->ResSlider->setValue(res_);

    dial->show();
}


//===========================================================================
void CurveResolutionSheet::ok()
//===========================================================================
{
    res_ = form_->ResSlider->value();
//     obs_->observedChanged();

    emit return_value(res_);
//     form_->close();
}


//=========================================================================== 
void CurveResolutionSheet::setHighRes()
//===========================================================================
{
    form_->ResSlider->setValue(def_high_res_);

    ok();
}
