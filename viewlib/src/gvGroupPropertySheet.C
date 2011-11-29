//===========================================================================
//                                                                           
// File: gvGroupPropertySheet.C                                        
//                                                                           
// Created: Wed Aug 28 12:22:31 2002                                         
//                                                                           
// Author: Sverre Briseid <sbr@sintef.no>
//                                                                           
// Revision:
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/gvGroupPropertySheet.h"
#include "GoTools/viewlib/ui_gvGroupPropertySheet_form.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/viewlib/gvRectangularSurfacePaintable.h"

#include "GoTools/viewlib/gvData.h"
#include "GoTools/viewlib/gvObserver.h"

// #include <q3groupbox.h>
#include <QPushButton>
#include <QSlider>
#include <QCheckBox>
// #include <q3table.h>
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
