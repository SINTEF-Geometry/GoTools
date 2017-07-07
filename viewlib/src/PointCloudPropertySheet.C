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

#include "GoTools/viewlib/PointCloudPropertySheet.h"
#include "GoTools/viewlib/gvPointCloudPaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

#include "GoTools/viewlib/ui_PointCloudPropertySheet_form.h"

#include <QPushButton>
#include <QSlider>
#include <QCheckBox>
#include <QLCDNumber>

#include <cmath>

namespace {
    const double MAX_POINT_SIZE = 10.0;
}

using namespace Ui;

//===========================================================================
PointCloudPropertySheet::PointCloudPropertySheet(gvPointCloudPaintable* pable)
    : pable_(pable), form_(), obs_(0)
{}


//===========================================================================
void PointCloudPropertySheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
    obs_ = obs;
    form_ = new Ui::PointCloudPropertySheet_form();//parent);

    QWidget* w = new QWidget();
    form_->setupUi(w);

    w->resize(form_->box->size());

    QObject::connect(form_->ApplyButton, SIGNAL(clicked()),
		     this, SLOT(apply()));
    QObject::connect(form_->CloseButton, SIGNAL(clicked()),
		     w, SLOT(close()));
// 		     form_, SLOT(close()));

    form_->VisibleCheck->setChecked(pable_->visible());
    form_->IdCheck->setChecked(pable_->getPaintId());
    double fraction = pable_->fractionRendered();
    int percentage = int(floor(fraction * 100.0) + 0.5);
    double ps = pable_->pointSize();
    int pspercent = int(ps/MAX_POINT_SIZE*100.0 + 0.5);
    form_->RenderSlider->setValue(percentage);
    form_->RenderLCDNumber->display(percentage);
    form_->PointsizeSlider->setValue(pspercent);
    form_->PointsizeLCDNumber->display(pspercent);

    w->show();
}



//===========================================================================
void PointCloudPropertySheet::apply()
//===========================================================================
{
    pable_->setVisible(form_->VisibleCheck->isChecked());
    pable_->setPaintId(form_->IdCheck->isChecked());
    int percentage = form_->RenderSlider->value();
    pable_->setFractionRendered(double(percentage)/100.0);
    double pointsize =
	double(form_->PointsizeSlider->value()) * MAX_POINT_SIZE / 100.0;
    pable_->setPointSize(pointsize);
    obs_->observedChanged();
}


