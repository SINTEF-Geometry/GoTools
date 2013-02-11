//===========================================================================
//                                                                           
// File: goview_vol_and_lr.C                                                 
//                                                                           
// Created: Fri Feb  8 11:27:13 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <utility>
#include <QApplication>

#include "GoTools/viewlib/vol_and_lr/DataHandlerVolAndLR.h"
#include "GoTools/viewlib/gvApplication.h"


/** An application for viewing spline surfaces and other geometrical objects.
 *
 *  The program is based on the QT toolkit, OpenGL and Sintef's Go toolkit.
 *
 *  This version adds support for volumes and LRSplines.
 *  
 */

//int less_largest_x_cntr=0, less_largest_y_cntr=0, less_smallest_y_cntr=0;

int main(int argc, char** argv)
{
//     // Use the generic app object
    QApplication theapp(argc, argv);

    // Create our main widget
    std::auto_ptr<DataHandler> dh(new DataHandlerVolAndLR);
    gvApplication* appwidget = new gvApplication(dh, NULL, argv[0]);

    appwidget->resize(500, 530);
//    appwidget->resize(1300, 1100);
    appwidget->show();
//     theapp.setMainWidget(appwidget);
    theapp.setActiveWindow(appwidget);

    int return_value=theapp.exec();
    return return_value;
}
