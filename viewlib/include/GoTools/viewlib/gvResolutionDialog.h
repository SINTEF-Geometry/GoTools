//===========================================================================
//                                                                           
// File: gvResolutionDialog.h                                                
//                                                                           
// Created: Wed Jul  4 09:47:53 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: gvResolutionDialog.h,v 1.2 2007-05-02 14:39:25 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVRESOLUTIONDIALOG_H
#define _GVRESOLUTIONDIALOG_H


#include <QDialog>
class QSlider;

/** Documentation ...
    etc
 */

class gvResolutionDialog : public QDialog
{

Q_OBJECT

public:
    gvResolutionDialog(int current_res_u, int current_res_v,
		       int minimum_res, int maximum_res,
		       QWidget * parent=0,
		       const char * name=0,
		       bool modal=FALSE,
		       Qt::WFlags f=0);
    virtual ~gvResolutionDialog();


signals:
    void valuesChanged(int, int);

protected slots:
    void apply();
    virtual void accept();

protected:
    int ures_;
    int vres_;
    QSlider* uslide_;
    QSlider* vslide_;
};



#endif // _GVRESOLUTIONDIALOG_H

