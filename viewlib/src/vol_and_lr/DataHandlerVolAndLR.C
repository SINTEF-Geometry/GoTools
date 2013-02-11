//===========================================================================
//                                                                           
// File: DataHandlerVolAndLR.C                                               
//                                                                           
// Created: Fri Feb  8 16:23:50 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================



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

