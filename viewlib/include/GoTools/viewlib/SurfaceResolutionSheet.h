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

#ifndef _SURFACERESOLUTIONSHEET_H
#define _SURFACERESOLUTIONSHEET_H


#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/viewlib/ui_SurfaceResolutionSheet_form.h"

#include <QObject>

class gvData;
//class RectangularSurfaceTesselator;
class gvRectangularSurfacePaintable;

/** Documentation ...
    etc
 */

class SurfaceResolutionSheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
SurfaceResolutionSheet(int ures = 20, int vres = 20) // The default res is 20x20.
     : form_(0), ures_(ures), vres_(vres), def_high_res_(200) //, obs_(0)
    {}
    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void ok();

signals:
    void return_value(int, int); // Return new resolution.


private:
    Ui::SurfaceResolutionSheet_form* form_;
    int ures_;
    int vres_;
    int def_high_res_; // Used in setHighRes().
/*     gvObserver* obs_; */

private slots:

    // Set default high resolution (optimized mode an advantage).
    // Accessed through button.
    void setHighRes();

};


#endif // _SPLINESURFACEPROPERTYSHEET_H

