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

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 5) {
    std::cout << "Input parameters : Input surface file(.g2), output file(.g2), boundary file (.csv) resolution" << std::endl;
    exit(-1);
  }
  // Read input arguments
  std::ifstream infile(argv[1]);
  std::ofstream outfile(argv[2]);
  std::ofstream outfile2(argv[3]);
  double res = atof(argv[4]);

  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;
  //Registrator<BoundedSurface> r210;

  // Read input surface
  ObjectHeader header;
  try {
    header.read(infile);
  }
  catch (...)
    {
      std::cerr << "Exiting" << std::endl;
      exit(-1);
    }
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  geom_obj->read(infile);
  
  shared_ptr<ParamSurface> sf = 
    dynamic_pointer_cast<ParamSurface, GeomObject>(geom_obj);
  if (!sf.get())
    {
      std::cerr << "Input file contains no surface" << std::endl;
      exit(-1);
    }
  // Fetch domain
  RectDomain domain = sf->containingDomain();
  double min_u = domain.umin();
  double max_u = domain.umax();
  double min_v = domain.vmin();
  double max_v = domain.vmax();

  // Compute number of samples
  int nmb_u = (int)((max_u - min_u)/(double)res) + 1;  
  int nmb_v = (int)((max_v - min_v)/(double)res) + 1;
  if ((nmb_u-1)*res < max_u - min_u)
    {
      // Adjust domain
      double mid = 0.5*(min_u + max_u);
      double nmb = (nmb_u-1)/2.0;
      min_u = mid - nmb*res;
      max_u = mid + nmb*res;
    }  
  if ((nmb_v-1)*res < max_v - min_v)
    {
      // Adjust domain
      double mid = 0.5*(min_v + max_v);
      double nmb = (nmb_v-1)/2.0;
      min_v = mid - nmb*res;
      max_v = mid + nmb*res;
    }

  // Compute grid
  double nodata_val = -9999;
  vector<double> val;
  shared_ptr<BoundedSurface> bd_sf = 
    dynamic_pointer_cast<BoundedSurface, ParamSurface>(sf);
  if (bd_sf.get())
    bd_sf->evalGrid(nmb_u, nmb_v, min_u, max_u, 
					 min_v, max_v,
					 val, nodata_val);
  else
    sf->evalGrid(nmb_u, nmb_v, min_u, max_u, min_v, max_v,
		 val, nodata_val);

  vector<double> data;
  int ki, kj, kr;
  if (sf->dimension() == 1)
    {
      double del_u = (max_u - min_u)/(double)(nmb_u-1);
      double del_v = (max_v - min_v)/(double)(nmb_v-1);
      data.reserve(3*val.size());
      for (kj=0, kr=0; kj<nmb_v; ++kj)
	for (ki=0; ki<nmb_u; ++ki, ++kr)
	  {
	    if (val[kr] == nodata_val)
	      continue;
	    data.push_back(min_u+ki*del_u);
	    data.push_back(min_v+kj*del_v);
	    data.push_back(val[kr]);
	  }
    }
  else
    {
      data.reserve(val.size());
      for (kr=0; kr<(int)val.size(); kr+=3)
	{
	  if (val[kr] == nodata_val)
	      continue;
	  data.insert(data.end(), val.begin()+kr, val.begin()+kr+3);
	}
    }
  
  
  vector<double> data2;
  if (bd_sf.get())
    {

      // Evaluate points along the surface boundary
      CurveLoop loop = bd_sf->outerBoundaryLoop();
      int nmb = loop.size();
      for (ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamCurve> cv = loop[ki];
	  shared_ptr<CurveOnSurface> sf_cv = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv);
	  double len;
	  if (sf_cv.get())
	    len = sf_cv->parameterCurve()->estimatedCurveLength();
	  else
	    len = cv->estimatedCurveLength();
	  int nsample = std::max((int)(len/res), 3);
	  double t1 = cv->startparam();
	  double t2 = cv->endparam();
	  double tdel = (t2 - t1)/(double)(nsample);
	  double tpar;
	  for (tpar=t1, kj=0; kj<nsample; ++kj, tpar+=tdel)
	    {
	      Point pos = cv->point(tpar);
	      if (sf->dimension() == 1)
		{
		  Point par = sf_cv->parameterCurve()->point(tpar);
		  data2.insert(data2.end(), par.begin(), par.end());
		  data2.push_back(pos[0]);
		}
	      else
		data2.insert(data2.end(), pos.begin(), pos.end());
	    }
	}
	 
      (void)outfile2.precision(15);
      outfile2 << "id,geo" << std::endl;
      outfile2 << "1, \"POLYGON ((";
      for (ki=0; ki<(int)data2.size(); ki+=3)
	{
	  //outfile2 << data2[ki-3] << "," << data2[ki-2] << "," << data2[ki-1] << ",";
	  outfile2 << data2[ki] << " " << data2[ki+1] << " " << data2[ki+2] << ",";
	}
      // outfile2 << data2[ki-3] << "," << data2[ki-2] << "," << data2[ki-1] << ",";
      // outfile2 << data2[0] << "," << data2[1] << "," << data2[2] << std::endl;
      outfile2 << data2[0] << " " << data2[1] << " " << data2[2] << "))\" " << std::endl;
    }

  data.insert(data.end(), data2.begin(), data2.end());
  PointCloud3D cloud(data.begin(), (int)data.size()/3);
  cloud.writeStandardHeader(outfile);
  cloud.write(outfile);
}
