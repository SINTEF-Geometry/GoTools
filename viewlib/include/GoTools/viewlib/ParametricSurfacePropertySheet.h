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

#ifndef _PARAMETRICSURFACEPROPERTYSHEET_H
#define _PARAMETRICSURFACEPROPERTYSHEET_H


#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/viewlib/ui_RectangularSurfacePropertySheet_form.h"
#include "GoTools/geometry/ParamSurface.h"
//#include "GoTools/viewlib/gvParametricSurfaceTesselator2.h"
#include "GoTools/viewlib/gvParametricSurfacePaintable.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"

#include <QObject>

class gvData;
//class ParametricSurfaceTesselator;
//class gvParametricSurfacePaintable;

/** Documentation ...
    etc
 */

class ParametricSurfacePropertySheet : public QObject, public gvPropertySheet
{

Q_OBJECT

public:
    ParametricSurfacePropertySheet()
    {}

    ParametricSurfacePropertySheet(Go::ParametricSurfaceTesselator* tess,
				   gvParametricSurfacePaintable* pable,
				   shared_ptr<Go::ParamSurface>& surf)
	: tess_(tess), pable_(pable), form_(0), obs_(0), surf_(surf)
    {}

    virtual ~ParametricSurfacePropertySheet();

    virtual void createSheet(QWidget* parent, gvObserver* obs);

public slots:
    void apply();
 

private:
  Go::ParametricSurfaceTesselator* tess_;
  gvParametricSurfacePaintable* pable_;
  Ui::RectangularSurfacePropertySheet_form* form_;
  gvObserver* obs_;
  shared_ptr<Go::ParamSurface> surf_;
};


#endif // _PARAMETRICSURFACEPROPERTYSHEET_H

