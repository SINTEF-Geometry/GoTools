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

#include "GoTools/viewlib/DefaultDataHandler.h"

#include "GoTools/viewlib/vol_and_lr/DataHandlerVolAndLR.h"
#include "GoTools/viewlib/vol_and_lr/RectangularVolumePropertySheet.h"
#include "GoTools/viewlib/vol_and_lr/gvRectangularVolumePaintable.h"

#include "GoTools/viewlib/gvRectangularSurfacePaintable.h"
#include "GoTools/viewlib/RectangularSurfacePropertySheet.h"

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/RectangularVolumeTesselator.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"

#include <GoTools/geometry/GoTools.h>
#include "GoTools/geometry/Factory.h"

using namespace Go;
using std::vector;
// using std::shared_ptr;
// using dynamic_pointer_cast;

//===========================================================================
DataHandlerVolAndLR::~DataHandlerVolAndLR()
  //===========================================================================
{
}

//===========================================================================
DataHandlerVolAndLR::DataHandlerVolAndLR()
  //===========================================================================
{
    // Create the default factory
    GoTools::init();
    Registrator<SplineVolume> r700;
    Registrator<LRSplineSurface> r293;

    //      cout << "Registering should be processed by now." << endl;
}


//===========================================================================
void DataHandlerVolAndLR::create(shared_ptr<GeomObject> obj,
				const gvColor& col, int id)
  //===========================================================================
{
    //    cout << "DataHandlerVolAndLR::create... " << flush;
    ClassType type = obj->instanceType();

    // Make unbounded elementary objects bounded with "canonical" parameter
    // bounds. This is ugly and arbitrary. We need this hack in order to
    // tesselate.
    switch (type) {
    case Class_LRSplineSurface:
      {
	const ParamSurface& sf
	  = dynamic_cast<const ParamSurface&>(*obj);
	if (sf.dimension() == 1)
	{
	    LRSplineSurface& lr_sf
		= dynamic_cast<LRSplineSurface&>(*obj);
#if 1
	    MESSAGE("Setting parameter domain to the unit square!");
	    lr_sf.setParameterDomain(0.0, 1.0, 0.0, 1.0);
#endif
	    MESSAGE("Lifting lrspline_sf from 1D to 3D.");
	    lr_sf.to3D();
	}

#if 1
	{
	    LRSplineSurface& lr_sf
		= dynamic_cast<LRSplineSurface&>(*obj);
	    MESSAGE("Translating the surface to the origin!");
	    BoundingBox bd_box = lr_sf.boundingBox();
	    Point mid_pt = 0.5*(bd_box.low() + bd_box.high());
	    Point transl_pt = -mid_pt;
	    std::cout << "transl_pt: (" << transl_pt[0] << ", " << transl_pt[1] << ", " << transl_pt[2] << std::endl;
	    lr_sf.translate(transl_pt);
	}
#endif

	shared_ptr<RectangularSurfaceTesselator> te(new RectangularSurfaceTesselator(sf));
	shared_ptr<gvRectangularSurfacePaintable> pa
	  (new gvRectangularSurfacePaintable(*(te->getMesh()), col, id));
	shared_ptr<ParamSurface> psf = 
	  dynamic_pointer_cast<ParamSurface, GeomObject>(obj);
	shared_ptr<gvPropertySheet> ps(new RectangularSurfacePropertySheet(te.get(), pa.get(), 
									   psf));
	tesselator_ = te;
	paintable_ = pa;
	property_sheet_ = ps;
	break;
      }
    case Class_SplineVolume:
      {
//	MESSAGE("SplineVolume support coming soon!");

	const SplineVolume& sv
	  = dynamic_cast<const SplineVolume&>(*obj);

	shared_ptr<RectangularVolumeTesselator> te(new RectangularVolumeTesselator(sv));
	shared_ptr<gvRectangularVolumePaintable> pa
	  (new gvRectangularVolumePaintable(*(te->getMesh()), col, id));
	shared_ptr<ParamVolume> pvol = 
	  dynamic_pointer_cast<ParamVolume, GeomObject>(obj);
	shared_ptr<gvPropertySheet> ps(new RectangularVolumePropertySheet(te.get(), pa.get(), 
									  pvol));
	tesselator_ = te;
	paintable_ = pa;
	property_sheet_ = ps;
	break;

      }
    default:
	DefaultDataHandler::create(obj, col, id);
    }
  //    cout << "finished" << endl;
}

