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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines3D/LRSplinePlotUtils3D.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/PointCloud.h"     // For visualization purposes
#include "GoTools/geometry/LineCloud.h"      // For visualization purposes

#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: lrspline_in (.g2) refinement_in lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ifstream filein2(argv[2]);
  std::ofstream fileout(argv[3]);

  GoTools::init();
  ObjectHeader oh;
  oh.read(filein);
      
  shared_ptr<LRSplineVolume> vol;
  if (oh.classType() == Class_SplineVolume)
    {
      double knot_tol = 1.0e-9;
      shared_ptr<SplineVolume> tmp(new SplineVolume);
      tmp->read(filein);
      vol = shared_ptr<LRSplineVolume>(new LRSplineVolume(tmp.get(), knot_tol));
    }
  else if (oh.classType() == Class_LRSplineVolume)
    {
      vol = shared_ptr<LRSplineVolume>(new LRSplineVolume());
      vol->read(filein);
    }

   std::ofstream ofm0("ref_vol_mesh0.g2");
  vector<double> mid0, bd0;
  extractElementMidAndBoundary(vol, mid0, bd0);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud0(mid0.begin(), mid0.size()/3);

  elem_cloud0.writeStandardHeader(ofm0);
  elem_cloud0.write(ofm0);

  // Make line cloud from element boundaries
  LineCloud elem_lines0(bd0.begin(), bd0.size()/6);
  elem_lines0.writeStandardHeader(ofm0);
  elem_lines0.write(ofm0);


  int nmb_refs;
  filein2 >> nmb_refs;
  for (int ki=0; ki<nmb_refs; ++ki)
    {
      double parval, start1, start2, end1, end2;
      int dir;
      int mult;

      filein2 >> parval;
      filein2 >> start1;
      filein2 >> end1;
      filein2 >> start2;
      filein2 >> end2;
      filein2 >> dir;
      filein2 >> mult;
      std::cout << "Iteration no. " << ki << std::endl;
      Direction3D d = (dir==0) ? XDIR : ((dir==1) ? YDIR : ZDIR);
      vol->refine(d, parval, start1, end1, start2, end2, mult, true);
    }

   std::ofstream ofm("ref_vol_mesh.g2");
  vector<double> mid, bd;
  extractElementMidAndBoundary(vol, mid, bd);
  
  // Make point cloud from mid parameters
  PointCloud3D elem_cloud(mid.begin(), mid.size()/3);

  elem_cloud.writeStandardHeader(ofm);
  elem_cloud.write(ofm);

  // Make line cloud from element boundaries
  LineCloud elem_lines(bd.begin(), bd.size()/6);
  elem_lines.writeStandardHeader(ofm);
  elem_lines.write(ofm);

  vol->writeStandardHeader(fileout);
  vol->write(fileout);
  return 0;
}
