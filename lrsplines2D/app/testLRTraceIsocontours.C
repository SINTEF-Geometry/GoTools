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

#include <fstream>
#include <iostream>
#include <chrono>

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"


using namespace std;
using namespace Go;

namespace {
  vector<double> contour_vals(const shared_ptr<ParamSurface>& surf, 
			      int num_contours);
}// end anonymous namespace


int main(int varnum, char* vararg[])
{

  if (varnum != 5 && varnum != 7)
    {
      std::cout << "<Surface> <number of planes> <tolerance> <maximum number of missing knots at finish> (min isoval, max isoval)" << std::endl;
      exit(-1);
    }
  
  // Create the default factory
  GoTools::init();
  Registrator<LRSplineSurface> r293;

  // Read surface
  ifstream is(vararg[1]);
  ObjectHeader header;
  header.read(is);
  shared_ptr<GeomObject> geom_obj(Factory::createObject(header.classType()));
  shared_ptr<ParamSurface> surf =
    dynamic_pointer_cast<ParamSurface,GeomObject>(geom_obj);
  if (!surf.get())
    {
      std::cout << "Object one is not a surface" << std::endl;
      exit(1);
    }
  surf->read(is);
  is.close();

  
  // Computing isocontours
  int nmb_iso = atoi(vararg[2]);
  double tol = atof(vararg[3]);
  int threshold_missing = atoi(vararg[4]);

  vector<double> isovals;
  if (varnum == 5)
    isovals = contour_vals(surf, nmb_iso); // 70
  else
    {
      isovals.resize(nmb_iso);
      double min = atof(vararg[5]);
      double max = atof(vararg[6]);
      double del = (nmb_iso == 1) ? 0 : (max - min)/(double)(nmb_iso-1);
      for (int ka=0; ka<nmb_iso; ++ka)
	isovals[ka] = min + (double)ka*del;
    }

  //const vector<double> isovals {1689.7636324695331};  // this isocontour causes topology problems with surface: data/256_lr_1d.g2
  auto t1 = chrono::high_resolution_clock::now();
  const vector<CurveVec> curves = LRTraceIsocontours(surf,
  						     isovals,
						     threshold_missing,
  						     tol);

  // const auto ssurf = lrsurf.asSplineSurface();
  // const vector<CurveVec> curves = SSurfTraceIsocontours(*ssurf,
  // 							isovals,
  // 							1e-5 * span,
  // 							include_3D,
  // 							use_sisl_marching);

  auto t2 = chrono::high_resolution_clock::now();
  cout << "Curves found in " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;
  
  ofstream os("curves.g2");
  ofstream os2("curves_2D.g2");
  //cout << "Number of curves found: " << endl;
  for (size_t i = 0; i != curves.size(); ++i) {
    //cout << "At height " << isovals[i] << ": " << curves[i].size() << " curves." << endl;
    for (auto cv : curves[i]) {
      if (cv.second.get())
	{
	  cv.second->writeStandardHeader(os);
	  cv.second->write(os);
	  cv.first->writeStandardHeader(os2);
	  cv.first->write(os2);
	}
    }
  }
  os.close();
  
  return 0;
}

namespace {


// =============================================================================
  vector<double> contour_vals(const shared_ptr<ParamSurface>& surf, 
			      int num_contours)
// =============================================================================
{
  BoundingBox bbox = surf->boundingBox();

  const double minval = bbox.low()[0];
  const double maxval = bbox.high()[0];
  
  vector<double> result(num_contours, 0);
  for (int i = 0; i != num_contours; ++i)
    result[i] = minval + ((maxval - minval)/(num_contours-1)) * i;

  return result;
}


} // end anonymous namespace
