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

#include "GoTools/viewlib/CurveResolutionSheet.h"
#include "GoTools/viewlib/ui_CurveResolutionSheet_form.h"

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
