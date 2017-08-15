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

#include "GoTools/viewlib/gvResolutionDialog.h"
#include <QLayout>
#include <QPushButton>
#include <QSlider>
#include <QLabel>



//===========================================================================
gvResolutionDialog::gvResolutionDialog(int current_res_u, int current_res_v,
				       int minimum_res, int maximum_res,
				       QWidget * parent,
				       const char * name,
				       bool modal,
				       Qt::WindowFlags f)
//===========================================================================
  : QDialog(parent),//, name, modal, f),
      ures_(current_res_u), vres_(current_res_v)
{
    // Make boxes
    //QVBoxLayout* vertbox = new QVBoxLayout(this);
    QHBoxLayout* uresbox = new QHBoxLayout();//vertbox);
    QHBoxLayout* vresbox = new QHBoxLayout();//vertbox);
    QHBoxLayout* buttonbox = new QHBoxLayout();//vertbox);
    // Make u box contents
//     uslide_ = new QSlider(minimum_res, maximum_res, 1, current_res_u,
// 			  Qt::Horizontal, this);
    uslide_ = new QSlider(Qt::Horizontal, this);
    uslide_->setRange(minimum_res, maximum_res);
    uslide_->setSingleStep(1);
    uslide_->setSliderPosition(current_res_u);
    QLabel* ulabel = new QLabel(this);
    ulabel->setNum(current_res_u);
    connect(uslide_, SIGNAL(valueChanged(int)),
	    ulabel, SLOT(setNum(int)));
    uresbox->addWidget(uslide_);
    uresbox->addWidget(ulabel);

    // Make v box contents
//     vslide_ = new QSlider(minimum_res, maximum_res, 1, current_res_v,
// 			  Qt::Horizontal, this);
//     QLabel* vlabel = new QLabel(this);
//     vlabel->setNum(current_res_v);
//     connect(vslide_, SIGNAL(valueChanged(int)),
// 	    vlabel, SLOT(setNum(int)));
//     vresbox->add(vslide_);
//     vresbox->add(vlabel);
    vslide_ = new QSlider(Qt::Horizontal, this);
    vslide_->setRange(minimum_res, maximum_res);
    vslide_->setSingleStep(1);
    vslide_->setSliderPosition(current_res_v);
    QLabel* vlabel = new QLabel(this);
    vlabel->setNum(current_res_v);
    connect(vslide_, SIGNAL(valueChanged(int)),
	    vlabel, SLOT(setNum(int)));
    vresbox->addWidget(vslide_);
    vresbox->addWidget(vlabel);

    // Make button box contents
    QPushButton* applybutton = new QPushButton("Change", this);
    connect(applybutton, SIGNAL(clicked()),
	    this, SLOT(apply()));
    QPushButton* okbutton = new QPushButton("Close", this);
    connect(okbutton, SIGNAL(clicked()),
	    this, SLOT(accept()));
    buttonbox->addWidget(applybutton);
    buttonbox->addWidget(okbutton);
}


//===========================================================================
gvResolutionDialog::~gvResolutionDialog()
//===========================================================================
{
}

//===========================================================================
void gvResolutionDialog::apply()
//===========================================================================
{
    if ((ures_ != uslide_->value()) || (vres_ != vslide_->value())) {
	ures_ = uslide_->value();
	vres_ = vslide_->value();
	emit valuesChanged(ures_, vres_);
    }
}

//===========================================================================
void gvResolutionDialog::accept()
//===========================================================================
{
    apply();
    QDialog::accept();
}
