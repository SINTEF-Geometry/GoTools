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

#ifndef _GVOBJECTLIST_H
#define _GVOBJECTLIST_H

// Qt includes
#include <QWidget>
// #include <qbuttongroup.h>
#include <Q3ButtonGroup>
#include <QLayout>
#include <QPushButton>
//Added by qt3to4:
// #include <Q3VBoxLayout>
#include <QScrollArea>

class gvData;
#include "GoTools/viewlib/gvObserver.h"

/** Documentation ...
    etc
 */

class gvObjectList : public QWidget, public gvObserver
{

Q_OBJECT

public:
    gvObjectList(gvData& data,
		 QWidget* parent=0, const char* name=0, Qt::WFlags f=0);
    virtual ~gvObjectList();
    virtual void observedChanged();
    void buildGUI();

protected slots:
    void clicked(int id);

private:
    gvData& data_;
    int numobj_;
//     QButtonGroup* bg_;
    Q3ButtonGroup* bg_;
    QVBoxLayout* lay1_;
    QVBoxLayout* lay2_;

    QScrollArea* scroll_area_;

    void setCorrectButtonStates();
};


#endif // _GVOBJECTLIST_H

