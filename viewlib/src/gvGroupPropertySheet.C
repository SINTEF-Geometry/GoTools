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

#include "GoTools/viewlib/gvGroupPropertySheet.h"
#include "GoTools/viewlib/ui_gvGroupPropertySheet_form.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/viewlib/gvRectangularSurfacePaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

#include <QPushButton>
#include <QSlider>
#include <QCheckBox>
#include <QLineEdit>
#include <QSpinBox>
#include <QComboBox>

using namespace std;


//===========================================================================
gvGroupPropertySheet::gvGroupPropertySheet(vector<int>& members, QString& def_name)
    : members_(members), def_name_(def_name)
//===========================================================================
{}

//===========================================================================
void gvGroupPropertySheet::createSheet(QWidget* parent, gvObserver* obs)
//===========================================================================
{
    obs_ = obs;

    form_ = new Ui::gvGroupPropertySheet_form();//parent);

    QWidget* w = new QWidget();
    form_->setupUi(w);
    w->resize(form_->box->size());

    form_->GroupName->setText(def_name_); // Initializing.

    connect(form_->ApplyButton, SIGNAL(clicked()),
	    this, SLOT(accept()));
    connect(form_->CloseButton, SIGNAL(clicked()),
	    w, SLOT(close()));
// 	    form_, SLOT(close()));

    w->show();
}

//===========================================================================
vector<int> gvGroupPropertySheet::getMembers()
//===========================================================================
{
    return members_;
}

//===========================================================================
void gvGroupPropertySheet::accept()
//===========================================================================
{
    // Here we are to read values off the form.
    QString name = form_->GroupName->text();
    obs_->observedChanged();

    emit value_changed(members_, name);
//     this->close();
//     form_->close();
}
