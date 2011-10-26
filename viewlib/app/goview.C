//===========================================================================
//                                                                           
// File: goview.C                                                                
//                                                                           
// Created: Thu Apr 26 13:15:50 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: goview.C,v 1.4 2009-01-29 12:57:46 jnygaard Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <utility>
#include <QApplication>

#include "GoTools/viewlib/DefaultDataHandler.h"
#include "GoTools/viewlib/gvApplication.h"


/** An application for viewing spline surfaces and other geometrical objects.
 *
 *  The program is based on the QT toolkit, OpenGL and Sintef's Go toolkit.
 *  
 */

int less_largest_x_cntr=0, less_largest_y_cntr=0, less_smallest_y_cntr=0;

int main(int argc, char** argv)
{
//     // Use the generic app object
    QApplication theapp(argc, argv);

    // Create our main widget
    std::auto_ptr<DataHandler> dh(new DefaultDataHandler);
    gvApplication* appwidget = new gvApplication(dh, NULL, argv[0]);

    appwidget->resize(500, 530);
//    appwidget->resize(1300, 1100);
    appwidget->show();
//     theapp.setMainWidget(appwidget);
    theapp.setActiveWindow(appwidget);

    int return_value=theapp.exec();
    return return_value;
}

