//===========================================================================
//                                                                           
// File: gvObjectList.h                                                      
//                                                                           
// Created: Tue Jul  3 13:54:11 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvObjectList.h,v 1.3 2007-05-02 14:39:24 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

