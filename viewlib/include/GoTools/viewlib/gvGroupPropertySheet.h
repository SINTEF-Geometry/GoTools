//===========================================================================
//                                                                           
// File: gvGroupPropertySheet.h                                        
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

#ifndef _GVGROUPPROPERTYSHEET_H
#define _GVGROUPPROPERTYSHEET_H



#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"

#include <QObject>
// #include <q3table.h>
#include <QComboBox>
#include <vector>
#include <QLabel>
#include "ui_gvGroupPropertySheet_form.h"

class gvData;
//class RectangularSurfaceTesselator;
class gvRectangularSurfacePaintable;
// class gvGroupPropertySheet_form;

/** Representation of several objects connected into one logical group.
 */

class gvGroupPropertySheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
    gvGroupPropertySheet(std::vector<int>& members, QString& def_name);

    virtual void createSheet(QWidget* parent, gvObserver* obs);

    std::vector<int> getMembers();

public slots:
    void accept();

signals:
 void value_changed(std::vector<int>& members, QString name);


private:
    Ui::gvGroupPropertySheet_form* form_;
    std::vector<int> members_; // Save members according to position in relevant vector.
    QString def_name_;
    gvObserver* obs_;
};


#endif // _GVGROUPPROPERTYSHEET_H

