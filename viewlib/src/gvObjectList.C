//===========================================================================
//                                                                           
// File: gvObjectList.C                                                      
//                                                                           
// Created: Tue Jul  3 13:56:53 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvObjectList.C,v 1.4 2007-05-23 15:44:48 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/viewlib/gvObjectList.h"
#include "GoTools/viewlib/gvData.h"
#include <QCheckBox>
#include <QLayout>
// #include <q3scrollview.h>
//Added by qt3to4:
// #include <Q3VBoxLayout>
// #include <Qt3Support/q3buttongroup.h>
#include <Qt3Support>

//===========================================================================
gvObjectList::gvObjectList(gvData& data,
			   QWidget* parent, const char* name, Qt::WFlags f)
//===========================================================================
//     : QWidget(parent, name, f),
    : QWidget(parent, f),
      data_(data),
      numobj_(0),
      bg_(0),
      lay1_(0),
      lay2_(0),
      scroll_area_(0)

{
    data.registerObserver(this);
    buildGUI();
    setCorrectButtonStates();
}


//===========================================================================
gvObjectList::~gvObjectList()
//===========================================================================
{
}


//===========================================================================
void gvObjectList::observedChanged()
//===========================================================================
{
   if (data_.numObjects() != numobj_)
      buildGUI();
   setCorrectButtonStates();
}


//===========================================================================
void gvObjectList::buildGUI()
//===========================================================================
{
    if (bg_) delete bg_;
    if (scroll_area_) delete scroll_area_;
    numobj_ = data_.numObjects();
//      std::cout << numobj_ << " objects." << std::endl;
//     bg_ = new QButtonGroup("Objects in model", this);
//     bg_ = new QButtonGroup(this);
    bg_ = new Q3ButtonGroup(this);
//     bg_->setFixedWidth(110);
//     bg_->setMaximumHeight(1000);

    scroll_area_ = new QScrollArea(this);
//     scroll_area_->setGeometry( 100, 100 , 1000 , 1000);
    scroll_area_->setFixedWidth(130);

    if (!lay1_)
      lay1_ = new QVBoxLayout(this);//, 2);
    lay1_->addWidget(scroll_area_);//bg_);
    lay2_ = new QVBoxLayout(bg_);//, 1);
    for (int i = 0; i < numobj_; ++i) {
	QString s = "Object " + QString::number(i);
	QCheckBox *cb=new QCheckBox(s, bg_);
	if (data_.object(i).get()==NULL)
	   cb->setHidden(true);
	lay2_->addWidget(cb);
    }

    scroll_area_->setWidget(bg_);
//     scroll_area_->setAlignment(Qt::AlignRight);
    connect(bg_, SIGNAL(clicked(int)),
	    this, SLOT(clicked(int)));
    bg_->show();
}

//===========================================================================
void gvObjectList::clicked(int id)
//===========================================================================
{
    // Get the state of button number id
//     bool bstate = bg_->find(id)->isOn();
    bool bstate = bg_->find(id)->isChecked();
    data_.setSelectedStateObject(id, bstate);
}


//===========================================================================
void gvObjectList::setCorrectButtonStates()
//===========================================================================
{
    for (int i = 0; i < numobj_; ++i) {
	if (data_.object(i).get()==NULL)
	{
	  bg_->find(i)->setHidden(true);
	  continue;
	}
	bool bstate = bg_->find(i)->isChecked();//On();
	if (bstate != data_.getSelectedStateObject(i))
	    bg_->find(i)->toggle();
    }
}
