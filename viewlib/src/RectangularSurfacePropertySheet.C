//===========================================================================
//                                                                           
// File: RectangularSurfacePropertySheet.C                                        
//                                                                           
// Created: Mon Jan  7 13:28:22 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectangularSurfacePropertySheet.C,v 1.3 2007-05-02 14:39:25 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/RectangularSurfacePropertySheet.h"
#include "ui_RectangularSurfacePropertySheet_form.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/viewlib/gvRectangularSurfacePaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

// #include <q3groupbox.h>
#include <QPushButton>
#include <QSlider>
#include <QCheckBox>

using namespace Ui;


//===========================================================================
RectangularSurfacePropertySheet::~RectangularSurfacePropertySheet()
//===========================================================================
{
   if (form_)
   {
      delete form_;
      form_=NULL;
   }
}

//===========================================================================
void RectangularSurfacePropertySheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
   if (form_)
   {
      delete form_;
      form_=NULL;
   }
   
    obs_ = obs;
    form_ = new Ui::RectangularSurfacePropertySheet_form();//NULL);
    //insertChild(form_);

    QWidget* w = new QWidget();
    form_->setupUi(w);

    w->resize(form_->box->size());

    QObject::connect(form_->ApplyButton, SIGNAL(pressed()),
	    this, SLOT(apply()));
    connect(form_->CloseButton, SIGNAL(pressed()),
	    w, SLOT(close()));

//     std::cout << "I'm here!" << std::endl;

    form_->VisibleCheck->setChecked(pable_->visible());
    int ures, vres;
    tess_->getRes(ures, vres);
    form_->UresSlider->setValue(ures);
    form_->VresSlider->setValue(vres);

    w->show();
}



//===========================================================================
void RectangularSurfacePropertySheet::apply()
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

