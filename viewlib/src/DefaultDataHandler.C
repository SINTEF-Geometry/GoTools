//===========================================================================
//                                                                           
// File: DefaultDataHandler.C                                                
//                                                                           
// Created: Fri Jan  4 14:57:25 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: DefaultDataHandler.C,v 1.2 2009-01-14 13:03:37 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/viewlib/DefaultDataHandler.h"

#include "GoTools/tesselator/CurveTesselator.h"
#include "GoTools/tesselator/RectangularSurfaceTesselator.h"
#include "GoTools/tesselator/NoopTesselator.h"
#include "GoTools/tesselator/LineCloudTesselator.h"
#include "GoTools/tesselator/RectGridTesselator.h"

#include "GoTools/viewlib/gvCurvePaintable.h"
#include "GoTools/viewlib/gvRectangularSurfacePaintable.h"
#include "GoTools/viewlib/gvPointCloudPaintable.h"
#include "GoTools/viewlib/gvLineCloudPaintable.h"
#include "GoTools/viewlib/gvQuadsPaintable.h"

#include "GoTools/viewlib/gvParametricSurfacePaintable.h"
#include "GoTools/tesselator/ParametricSurfaceTesselator.h"
#include "GoTools/viewlib/ParametricSurfacePropertySheet.h"

#include "GoTools/viewlib/gvPropertySheet.h"
#include "GoTools/viewlib/SplineCurvePropertySheet.h"
#include "GoTools/viewlib/RectangularSurfacePropertySheet.h"
#include "GoTools/viewlib/PointCloudPropertySheet.h"

#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/GoTools.h"

using namespace Go;
using std::vector;
// using std::shared_ptr;
// using dynamic_pointer_cast;

//===========================================================================
DefaultDataHandler::~DefaultDataHandler()
  //===========================================================================
{
}

//===========================================================================
DefaultDataHandler::DefaultDataHandler()
  //===========================================================================
{
    // Create the default factory
    GoTools::init();

    //      cout << "Registering should be processed by now." << endl;
}


//===========================================================================
void DefaultDataHandler::create(shared_ptr<GeomObject> obj,
				const gvColor& col, int id)
  //===========================================================================
{
  //    cout << "DefaultDataHandler::create... " << flush;
  ClassType type = obj->instanceType();
  switch (type)
    {
    case Class_SplineCurve:
    case Class_CurveOnSurface:
    case Class_BoundedCurve:
    case Class_Ellipse:
      {
	const ParamCurve& cv
	  = dynamic_cast<const ParamCurve&>(*obj);
	shared_ptr<CurveTesselator> te(new CurveTesselator(cv));
	shared_ptr<gvCurvePaintable> pa(new gvCurvePaintable(*(te->getMesh()), col, id));
	shared_ptr<gvPropertySheet> ps(new SplineCurvePropertySheet(te.get(), pa.get()));
	tesselator_ = te;
	paintable_ = pa;
	property_sheet_ = ps;
	break;
      }
    case Class_SplineSurface:
    case Class_Plane:
    case Class_Cylinder:
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
    case Class_BoundedSurface:
      {
	// We analyze the loops to see if they are valid.
	BoundedSurface& bd_sf = dynamic_cast<BoundedSurface&>(*obj);
	bd_sf.analyzeLoops();
	int valid_state = 0;
	if (!bd_sf.isValid(valid_state)) {
	    MESSAGE("Bounded surface not valid, state: " <<
		    valid_state << ". Trying to fix.");
	    shared_ptr<BoundedSurface> bd_sf_ptr =
		dynamic_pointer_cast<BoundedSurface>(obj);
	    BoundedUtils::fixInvalidBoundedSurface(bd_sf_ptr);
	    int state = 0;
	    bool sf_ok = bd_sf.isValid(state);
	    if (!sf_ok)
		MESSAGE("Failed fixing bd_sf!");
	}
	const ParamSurface& sf
	  = dynamic_cast<const ParamSurface&>(*obj);
	shared_ptr<ParametricSurfaceTesselator> te(new ParametricSurfaceTesselator(sf));
	shared_ptr<gvParametricSurfacePaintable> pa
	  (new gvParametricSurfacePaintable(*(te->getMesh()), col, id));
	shared_ptr<ParamSurface> psf = 
	  dynamic_pointer_cast<ParamSurface, GeomObject>(obj);
	shared_ptr<gvPropertySheet> ps(new ParametricSurfacePropertySheet(te.get(), pa.get(), 
									  psf));
	tesselator_ = te;
	paintable_ = pa;
	property_sheet_ = ps;
	break;
      }
//     case Class_Plane:
// 	{
// 	    tesselator_ = shared_ptr<Tesselator>();
// 	    paintable_ = shared_ptr<gvPaintable>();
// 	    property_sheet_ = shared_ptr<gvPropertySheet>();
// 	    break;
// 	}
    case Class_PointCloud:
      {
	const PointCloud3D& cl
	  = dynamic_cast<const PointCloud3D&>(*obj);
	shared_ptr<Tesselator> te(new NoopTesselator);
	shared_ptr<gvPointCloudPaintable> pa(new gvPointCloudPaintable(cl, col, id));
	shared_ptr<gvPropertySheet> ps(new PointCloudPropertySheet(pa.get()));
	tesselator_ = te;
	paintable_ = pa;
	property_sheet_ = ps;
	break;
      }
    case Class_LineCloud:
      {
	const LineCloud& cl
	  = dynamic_cast<const LineCloud&>(*obj);
	shared_ptr<LineCloudTesselator> te(new LineCloudTesselator(cl));
	shared_ptr<gvLineCloudPaintable> pa(new gvLineCloudPaintable(te->getRenderCloud(), col, id));
	shared_ptr<gvPropertySheet> ps;
	tesselator_ = te;
	paintable_ = pa;
	property_sheet_ = ps;
	break;
      }
    case Class_RectGrid:
      {
	const RectGrid& rg
	  = dynamic_cast<const RectGrid&>(*obj);
	shared_ptr<RectGridTesselator> te(new RectGridTesselator(rg));
	shared_ptr<gvQuadsPaintable> pa(new gvQuadsPaintable(*(te->getMesh()), col, id));
	shared_ptr<gvPropertySheet> ps;
	tesselator_ = te;
	paintable_ = pa;
	property_sheet_ = ps;
	break;
      }
    default:
      THROW("No such type code: " << type);
    }
  //    cout << "finished" << endl;
}

