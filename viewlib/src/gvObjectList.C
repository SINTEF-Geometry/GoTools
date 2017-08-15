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

#include "GoTools/viewlib/gvObjectList.h"
#include "GoTools/viewlib/gvData.h"
#include <QCheckBox>
#include <QLayout>

//===========================================================================
gvObjectList::gvObjectList(gvData& data,
			   QWidget* parent, const char* name, Qt::WindowFlags f)
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
    bg_ = new QButtonGroup(this); // Used for accessing the values, not connected with the GUI.
    
    scroll_area_ = new QScrollArea(this);
    QWidget * scroll_area_widget_contents = new QWidget;
    scroll_area_->setFixedWidth(130);

    if (!lay1_)
        lay1_ = new QVBoxLayout(this);//, 2);
    lay1_->addWidget(scroll_area_);//bg_);
    lay2_ = new QVBoxLayout();//scroll_area_widget);//, 1);
    scroll_area_widget_contents->setLayout(lay2_);
    for (int i = 0; i < numobj_; ++i) {
	QString s = "Object " + QString::number(i);
	QCheckBox *cb=new QCheckBox(s);//, bg_);
        cb->setMinimumHeight(16);
	if (data_.object(i).get()==NULL)
	   cb->setHidden(true);
        lay2_->addWidget(cb);
	bg_->addButton(cb, i);
    }
    lay2_->insertStretch(-1, 1); // Remove space between elements (in the list of GeomObject's).

    scroll_area_->setWidgetResizable(true);
    scroll_area_->setWidget(scroll_area_widget_contents);
        
    connect(bg_, SIGNAL(buttonClicked(int)),
            this, SLOT(clicked(int)));
    bg_->setExclusive(false);

}

//===========================================================================
void gvObjectList::clicked(int id)
//===========================================================================
{
    // Get the state of button number id
//     bool bstate = bg_->find(id)->isOn();
    QAbstractButton* button = bg_->button(id);
    bool bstate = button->isChecked();//->find(id)->isEnabled();//isChecked();
    // MESSAGE("clicked(): id: " << id << ", is down: " << bstate);
    data_.setSelectedStateObject(id, bstate);
}


//===========================================================================
void gvObjectList::setCorrectButtonStates()
//===========================================================================
{    
    for (int i = 0; i < numobj_; ++i) {
	if (data_.object(i).get()==NULL)
	{
            QAbstractButton* button = bg_->button(i);
            button->hide();
            continue;
	}
        QAbstractButton* button = bg_->button(i);
        bool checked = button->isChecked();//->find(id)->isEnabled();//isChecked();
        // MESSAGE("setCorrectButtonStates(): Button number " << i << " checked status: " << checked);
        if (checked != data_.getSelectedStateObject(i))
        {
            button->toggle();
        }
    }
}
