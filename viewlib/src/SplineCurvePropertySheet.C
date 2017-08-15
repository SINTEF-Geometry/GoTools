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

#include "GoTools/viewlib/SplineCurvePropertySheet.h"
#include "GoTools/viewlib/ui_SplineCurvePropertySheet_form.h"
#include "GoTools/tesselator/CurveTesselator.h"
#include "GoTools/viewlib/gvCurvePaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

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

