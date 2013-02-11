//===========================================================================
//                                                                           
// File: RectangularVolumePropertySheet.C                                    
//                                                                           
// Created: Thu Jul  5 17:09:35 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================



#include "GoTools/viewlib/vol_and_lr/RectangularVolumePropertySheet.h"
#include "GoTools/viewlib/vol_and_lr/ui_RectangularVolumePropertySheet_form.h"
#include "GoTools/viewlib/vol_and_lr/gvRectangularVolumePaintable.h"

#include "GoTools/trivariate/RectangularVolumeTesselator.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

// #include <q3groupbox.h>
#include <QPushButton>
#include <QSlider>
#include <QCheckBox>

using namespace Ui;


//===========================================================================
RectangularVolumePropertySheet::~RectangularVolumePropertySheet()
//===========================================================================
{
   if (form_)
   {
      delete form_;
      form_=NULL;
   }
}

//===========================================================================
void RectangularVolumePropertySheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
   if (form_)
   {
      delete form_;
      form_=NULL;
   }
   
    obs_ = obs;
    form_ = new Ui::RectangularVolumePropertySheet_form();//NULL);
    //insertChild(form_);

    QDialog* w = new QDialog();
    form_->setupUi(w);

    w->resize(form_->box->size());

    QObject::connect(form_->ApplyButton, SIGNAL(pressed()),
	    this, SLOT(apply()));
    connect(form_->CloseButton, SIGNAL(pressed()),
	    w, SLOT(close()));

//     std::cout << "I'm here!" << std::endl;

#if 0
    form_->VisibleCheck->setChecked(pable_->visible());
#endif

    int res;
    tess_->getRes(res);
    form_->ResSlider->setValue(res);

    w->show();
}



//===========================================================================
void RectangularVolumePropertySheet::apply()
//===========================================================================
{
    int res = form_->ResSlider->value();
#if 0
    int vres = form_->VresSlider->value();
    if (form_->TurnOrientationCheck->isChecked()) {
	surf_->turnOrientation();
	form_->TurnOrientationCheck->setChecked(false);
    }
    pable_->setVisible(form_->VisibleCheck->isChecked());
#endif
    tess_->changeRes(res);
    obs_->observedChanged();
}

