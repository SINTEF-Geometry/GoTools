//===========================================================================
//                                                                           
// File: PointCloudPropertySheet.C                                           
//                                                                           
// Created: Wed May 29 14:50:08 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PointCloudPropertySheet.C,v 1.3 2007-05-02 14:39:25 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================



#include "GoTools/viewlib/PointCloudPropertySheet.h"
#include "GoTools/viewlib/gvPointCloudPaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

#include "ui_PointCloudPropertySheet_form.h"

// #include <q3groupbox.h>
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


