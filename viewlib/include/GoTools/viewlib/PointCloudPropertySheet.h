//===========================================================================
//                                                                           
// File: PointCloudPropertySheet.h                                           
//                                                                           
// Created: Wed May 29 14:48:33 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PointCloudPropertySheet.h,v 1.2 2007-05-02 14:39:23 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _POINTCLOUDPROPERTYSHEET_H
#define _POINTCLOUDPROPERTYSHEET_H


#include "GoTools/viewlib/gvPropertySheet.h"
#include "ui_PointCloudPropertySheet_form.h"

#include <QObject>

class gvData;
// class PointCloudPropertySheet_form;
class gvPointCloudPaintable;

/** Documentation ...
    The property sheet for PointCloud.
 */

class PointCloudPropertySheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
PointCloudPropertySheet(gvPointCloudPaintable* pable);
    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void apply();
 

private:
    gvPointCloudPaintable* pable_;
    Ui::PointCloudPropertySheet_form* form_;
    gvObserver* obs_;
};



#endif // _POINTCLOUDPROPERTYSHEET_H

