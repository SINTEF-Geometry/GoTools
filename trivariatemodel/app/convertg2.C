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

// #include "GoTools/trivariate/SplineVolume.h"
// #include "GoTools/trivariate/CylinderVolume.h"
// #include "GoTools/trivariate/SphereVolume.h"
// #include "GoTools/trivariate/ConeVolume.h"
// #include "GoTools/trivariate/Parallelepiped.h"
#include "GoTools/trivariate/CurveOnVolume.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
// #include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
// #include "GoTools/geometry/Cylinder.h"
// #include "GoTools/geometry/Plane.h"
// #include "GoTools/geometry/SurfaceOfLinearExtrusion.h"
// #include "GoTools/geometry/Line.h"
// #include "GoTools/geometry/Cone.h"
// #include "GoTools/geometry/Circle.h"
// #include "GoTools/geometry/Sphere.h"
// #include "GoTools/geometry/Torus.h"
// #include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/igeslib/IGESconverter.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include <fstream>

//using namespace std;
using namespace Go;
using std::vector;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Usage : Input file (g2 or g22 format), File type (1=g2, 2=g22), Output file (g2)" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::string infile(argv[1]);

  int file_type = atoi(argv[2]);

  std::ofstream outfile(argv[3]);

  GoTools::init();
  Registrator<SurfaceOnVolume> r211;
  Registrator<CurveOnVolume> r111;

  vector<shared_ptr<GeomObject> > geom;
  if (file_type == 2)
    {
      VolumeModelFileHandler filehandler;
      vector<shared_ptr<ParamSurface> > sfs = filehandler.readSurface(infile.c_str());
      geom.insert(geom.end(), sfs.begin(), sfs.end());
    }
  else
    {
      std::ifstream is(infile);
      IGESconverter conv;
      conv.readgo(is);
      geom = conv.getGoGeom();
    }

  vector<shared_ptr<GeomObject> > geom2;
  for (size_t ki=0; ki<geom.size(); ++ki)
    {
      if (geom[ki]->instanceType() == Class_BoundedSurface)
	{
	  shared_ptr<BoundedSurface> bd_sf =
	    dynamic_pointer_cast<BoundedSurface, GeomObject>(geom[ki]);
	  if (bd_sf.get())
	    {
	      shared_ptr<ParamSurface> surf = bd_sf->underlyingSurface();
	      if (surf->instanceType() == Class_SurfaceOnVolume)
		{
		  shared_ptr<SurfaceOnVolume> vol_sf =
		    dynamic_pointer_cast<SurfaceOnVolume, GeomObject>(surf);
		  surf = vol_sf->spaceSurface();
		  vector<CurveLoop> bd_loops = bd_sf->allBoundaryLoops();
		  for (size_t kj=0; kj<bd_loops.size(); ++kj)
		    {
		      int nmb = bd_loops[kj].size();
		      for (int kr=0; kr<nmb; ++kr)
			{
			  shared_ptr<ParamCurve> cv = bd_loops[kj][kr];
			  shared_ptr<CurveOnSurface> sf_cv =
			    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(cv);
			  if (sf_cv.get())
			    {
			      shared_ptr<CurveOnVolume> vol_cv =
				dynamic_pointer_cast<CurveOnVolume, ParamCurve>(sf_cv->spaceCurve());
			      if (vol_cv.get())
				sf_cv->setSpaceCurve(vol_cv->spaceCurve());
			      sf_cv->setUnderlyingSurface(surf);
			    }
			}
		    }
		  bd_sf->replaceSurf(surf);
		}
	      geom2.push_back(bd_sf);
	    }
	}
      else if (geom[ki]->instanceType() == Class_CurveOnVolume)
	{
	  shared_ptr<CurveOnVolume> vol_cv =
	    dynamic_pointer_cast<CurveOnVolume, GeomObject>(geom[ki]);
	  geom2.push_back(vol_cv->spaceCurve());
	}
      else if (geom[ki]->instanceType() == Class_SurfaceOnVolume)
	{
	  shared_ptr<SurfaceOnVolume> vol_sf =
	    dynamic_pointer_cast<SurfaceOnVolume, GeomObject>(geom[ki]);
	  geom2.push_back(vol_sf->spaceSurface());
	}
      else
	geom2.push_back(geom[ki]);
    }

  for (size_t ki=0; ki<geom2.size(); ++ki)
    {
      geom2[ki]->writeStandardHeader(outfile);
      geom2[ki]->write(outfile);
    }
}


  
