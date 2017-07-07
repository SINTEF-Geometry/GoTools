/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/viewlib/vol_and_lr/RectangularVolumePropertySheet.h"
#include "GoTools/viewlib/vol_and_lr/ui_RectangularVolumePropertySheet_form.h"
#include "GoTools/viewlib/vol_and_lr/gvRectangularVolumePaintable.h"

#include "GoTools/trivariate/RectangularVolumeTesselator.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

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

