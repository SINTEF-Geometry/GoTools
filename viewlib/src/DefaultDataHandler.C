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

// #include "GoTools/lrsplines2D/LRSplineSurface.h"

// #include "GoTools/trivariate/SplineVolume.h"
// #include "GoTools/trivariate/RectangularVolumeTesselator.h"
// #include "GoTools/viewlib/volume/RectangularVolumePropertySheet.h"
// #include "GoTools/viewlib/volume/gvRectangularVolumePaintable.h"

#include <GoTools/geometry/GoTools.h>
#include "GoTools/geometry/Factory.h"

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
//    Registrator<SplineVolume> r700;
//    Registrator<LRSplineSurface> r293;

    //      cout << "Registering should be processed by now." << endl;
}


//===========================================================================
void DefaultDataHandler::create(shared_ptr<GeomObject> obj,
				const gvColor& col, int id)
  //===========================================================================
{
    //    cout << "DefaultDataHandler::create... " << flush;
    if (obj.get() == NULL)
    {
	return;
    }

    ClassType type = obj->instanceType();

    // Make unbounded elementary objects bounded with "canonical" parameter
    // bounds. This is ugly and arbitrary. We need this hack in order to
    // tesselate.
    switch (type) {
        case Class_Line:
            {
                shared_ptr<Line> line 
                    = dynamic_pointer_cast<Line, GeomObject>(obj);
                if (!line->isBounded()) {
                    line->setParamBounds(-1.0, 1.0);
                }
                break;
            }
        case Class_Hyperbola:
            {
                shared_ptr<Hyperbola> hyperbola 
                    = dynamic_pointer_cast<Hyperbola, GeomObject>(obj);
                if (!hyperbola->isBounded()) {
                    hyperbola->setParamBounds(-1.0, 1.0);
                }
                break;
            }
        case Class_Parabola:
            {
                shared_ptr<Parabola> parabola 
                    = dynamic_pointer_cast<Parabola, GeomObject>(obj);
                if (!parabola->isBounded()) {
                    parabola->setParamBounds(-1.0, 1.0);
                }
                break;
            }
        case Class_Plane:
            {
                shared_ptr<Plane> plane 
                    = dynamic_pointer_cast<Plane, GeomObject>(obj);
                if (!plane->isBounded()) {
                    plane->setParameterBounds(-1.0, -1.0, 1.0, 1.0);
                }
                break;
            }
        case Class_Cylinder:
            {
                shared_ptr<Cylinder> cylinder
                    = dynamic_pointer_cast<Cylinder, GeomObject>(obj);
                if (!cylinder->isBounded()) {
                    cylinder->setParamBoundsV(-1.0, 1.0);
                }
                break;
            }
        case Class_Cone:
            {
                shared_ptr<Cone> cone 
                    = dynamic_pointer_cast<Cone, GeomObject>(obj);
                if (!cone->isBounded()) {
                    cone->setParamBoundsV(-1.0, 1.0);
                }
                break;
            }
    }


  switch (type)
    {
    case Class_SplineCurve:
    case Class_CurveOnSurface:
    case Class_BoundedCurve:
    case Class_Line:
    case Class_Circle:
    case Class_Ellipse:
    case Class_Hyperbola:
    case Class_Parabola:
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
    case Class_Sphere:
    case Class_Cone:
    case Class_Torus:
    case Class_SurfaceOfRevolution:
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
    // case Class_LRSplineSurface:
    //   {
    // 	const ParamSurface& sf
    // 	  = dynamic_cast<const ParamSurface&>(*obj);
    // 	shared_ptr<RectangularSurfaceTesselator> te(new RectangularSurfaceTesselator(sf));
    // 	shared_ptr<gvRectangularSurfacePaintable> pa
    // 	  (new gvRectangularSurfacePaintable(*(te->getMesh()), col, id));
    // 	shared_ptr<ParamSurface> psf = 
    // 	  dynamic_pointer_cast<ParamSurface, GeomObject>(obj);
    // 	shared_ptr<gvPropertySheet> ps(new RectangularSurfacePropertySheet(te.get(), pa.get(), 
    // 									   psf));
    // 	tesselator_ = te;
    // 	paintable_ = pa;
    // 	property_sheet_ = ps;
    // 	break;
    //   }
    case Class_BoundedSurface:
      {
	// We analyze the loops to see if they are valid.
	BoundedSurface& bd_sf = dynamic_cast<BoundedSurface&>(*obj);
	bd_sf.analyzeLoops();
	int valid_state = 0;
	if (!bd_sf.isValid(valid_state)) {
	    MESSAGE("Bounded surface (id = " << id << ") not valid, state: " <<
		    valid_state << ". Trying to fix.");
	    shared_ptr<BoundedSurface> bd_sf_ptr =
		dynamic_pointer_cast<BoundedSurface>(obj);
	    double max_tol_mult = 1.3;
	    BoundedUtils::fixInvalidBoundedSurface(bd_sf_ptr, max_tol_mult);
	    int state = 0;
	    bool sf_ok = bd_sf.isValid(state);
	    if (!sf_ok)
		MESSAGE("Failed fixing bd_sf!");
	    else
		MESSAGE("Fixed bd_sf!");
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
//     case Class_SplineVolume:
//       {
// //	MESSAGE("SplineVolume support coming soon!");

// 	const SplineVolume& sv
// 	  = dynamic_cast<const SplineVolume&>(*obj);

// 	shared_ptr<RectangularVolumeTesselator> te(new RectangularVolumeTesselator(sv));
// 	shared_ptr<gvRectangularVolumePaintable> pa
// 	  (new gvRectangularVolumePaintable(*(te->getMesh()), col, id));
// 	shared_ptr<ParamVolume> pvol = 
// 	  dynamic_pointer_cast<ParamVolume, GeomObject>(obj);
// 	shared_ptr<gvPropertySheet> ps(new RectangularVolumePropertySheet(te.get(), pa.get(), 
// 									  pvol));
// 	tesselator_ = te;
// 	paintable_ = pa;
// 	property_sheet_ = ps;
// 	break;

//       }
    default:
      THROW("No such type code: " << type);
    }
  //    cout << "finished" << endl;
}

