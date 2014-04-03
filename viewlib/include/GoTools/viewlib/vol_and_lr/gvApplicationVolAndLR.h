//===========================================================================
//                                                                           
// File: gvApplicationVolAndLR.h                                             
//                                                                           
// Created: Thu Jun 20 08:26:56 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GVAPPLICATIONVOLANDLR_H
#define _GVAPPLICATIONVOLANDLR_H


#include "GoTools/viewlib/gvApplication.h"


class gvApplicationVolAndLR : public gvApplication
{

Q_OBJECT


public:
    gvApplicationVolAndLR(std::auto_ptr<DataHandler> dh,
			  QWidget * parent=0,
			  const char * name=0,
			  Qt::WFlags f=0);

    virtual ~gvApplicationVolAndLR();

public slots:
    void translate_to_origin(); // All selected objects are translated by the center of their bounding box.

protected:
    void buildExtraGUI();

};

#endif // _GVAPPLICATIONVOLANDLR_H

