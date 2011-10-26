//===========================================================================
//                                                                           
// File: SplineCurvePropertySheet.C                                          
//                                                                           
// Created: Thu Jan 31 13:06:10 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SplineCurvePropertySheet.C,v 1.3 2007-05-02 14:39:26 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/SplineCurvePropertySheet.h"
#include "ui_SplineCurvePropertySheet_form.h"
#include "GoTools/tesselator/CurveTesselator.h"
#include "GoTools/viewlib/gvCurvePaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

// #include <q3groupbox.h>
#include <QPushButton>
#include <QSlider>
#include <QCheckBox>
#include <QSpinBox>

using namespace Ui;
using namespace Go;

//===========================================================================
void SplineCurvePropertySheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
    obs_ = obs;
    form_ = new Ui::SplineCurvePropertySheet_form();//parent);

    QWidget* w = new QWidget();
    form_->setupUi(w);

    w->resize(form_->box->size());

    connect(form_->ApplyButton, SIGNAL(clicked()),
	    this, SLOT(apply()));
    connect(form_->CloseButton, SIGNAL(clicked()),
	    w, SLOT(close()));

    form_->VisibleCheck->setChecked(pable_->visible());
    int res;
    tess_->getRes(res);
    form_->ResSlider->setValue(res);

    gvColor col = pable_->getNormalColor();
    form_->redSpinBox->setValue(int(255.0*col.rgba[0]));
    form_->greenSpinBox->setValue(int(255.0*col.rgba[1]));
    form_->blueSpinBox->setValue(int(255.0*col.rgba[2]));
    form_->alphaSpinBox->setValue(int(255.0*col.rgba[3]));

    w->show();
}


//===========================================================================
void SplineCurvePropertySheet::apply()
//===========================================================================
{
    pable_->setVisible(form_->VisibleCheck->isChecked());
    tess_->changeRes(form_->ResSlider->value());
    obs_->observedChanged();
    gvColor new_col((float)form_->redSpinBox->value()/255.0f,
		    (float)form_->greenSpinBox->value()/255.0f,
		    (float)form_->blueSpinBox->value()/255.0f,
		    (float)form_->alphaSpinBox->value()/255.0f);
    pable_->setColor(new_col);
}

