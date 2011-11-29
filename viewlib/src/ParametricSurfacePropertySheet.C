//===========================================================================
//                                                                           
// File: ParametricSurfacePropertySheet.C                                        
//                                                                           
// Created: Mon Jan  7 13:28:22 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ParametricSurfacePropertySheet.C,v 1.3 2007-05-02 14:39:25 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/ParametricSurfacePropertySheet.h"
#include "GoTools/viewlib/ui_RectangularSurfacePropertySheet_form.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include "GoTools/viewlib/gvParametricSurfacePaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

// #include <q3groupbox.h>
#include <QPushButton>
#include <QSlider>
#include <QCheckBox>

using namespace Ui;


//===========================================================================
ParametricSurfacePropertySheet::~ParametricSurfacePropertySheet()
//===========================================================================
{
   if (form_)
   {
      delete form_;
      form_=NULL;
   }
}


//===========================================================================
void ParametricSurfacePropertySheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
    obs_ = obs;
    form_ = new Ui::RectangularSurfacePropertySheet_form();//NULL); //parent);

    QWidget* w = new QWidget();
    form_->setupUi(w);

    w->resize(form_->box->size());

    QObject::connect(form_->ApplyButton, SIGNAL(clicked()),
		     this, SLOT(apply()));
    QObject::connect(form_->CloseButton, SIGNAL(clicked()),
		     w, SLOT(close()));

    form_->VisibleCheck->setChecked(pable_->visible());
    int ures, vres;
    tess_->getRes(ures, vres);
    form_->UresSlider->setValue(ures);
    form_->VresSlider->setValue(vres);

    w->show();
}


//===========================================================================
void ParametricSurfacePropertySheet::apply()
//===========================================================================
{
    int ures = form_->UresSlider->value();
    int vres = form_->VresSlider->value();
    if (form_->TurnOrientationCheck->isChecked()) {
	surf_->turnOrientation();
	form_->TurnOrientationCheck->setChecked(false);
    }
    pable_->setVisible(form_->VisibleCheck->isChecked());
    tess_->changeRes(ures, vres);
    obs_->observedChanged();
}

