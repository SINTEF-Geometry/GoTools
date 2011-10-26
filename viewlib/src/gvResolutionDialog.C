//===========================================================================
//                                                                           
// File: gvResolutionDialog.C                                                
//                                                                           
// Created: Wed Jul  4 09:56:27 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvResolutionDialog.C,v 1.2 2007-05-02 14:39:26 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/viewlib/gvResolutionDialog.h"
// #include <q3vbox.h>
// #include <q3hbox.h>
#include <QLayout>
#include <QPushButton>
#include <QSlider>
#include <QLabel>
//Added by qt3to4:
// #include <Q3VBoxLayout>
// #include <Q3HBoxLayout>



//===========================================================================
gvResolutionDialog::gvResolutionDialog(int current_res_u, int current_res_v,
				       int minimum_res, int maximum_res,
				       QWidget * parent,
				       const char * name,
				       bool modal,
				       Qt::WFlags f)
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
