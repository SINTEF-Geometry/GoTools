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

#include <iostream>
#include <fstream>

#include "GoTools/geometry/GoTools.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"

using namespace std;
using namespace Go;


int main (int argc, char *argv[]) {

  if (argc != 9)
    {
      cout << "usage: ./segmentContours <input lrvol(.g2)> <output contours (.g2)> <par. dir.> <nmb sections> <min. par. val.> <max. par. val> <isoval> <tol>" << endl;
      return -1;
    }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  int dir = atoi(argv[3]);
  int nmb_sec = atoi(argv[4]);
  double minval = atof(argv[5]);
  double maxval = atof(argv[6]);
  double isoval = atof(argv[7]);
  double tol = atof(argv[8]);
  int threshold_missing = 100;

  GoTools::init();
  ObjectHeader oh;
  oh.read(ifs);
      
  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifs);

  std::cout << "LR Volume read" << std::endl;
  
  vector<double> iso;
  iso.push_back(isoval);
  double del = (nmb_sec == 1) ? 0 : (maxval - minval)/(double)(nmb_sec-1);
  vector<double> secs;
  for (int ka=0; ka<nmb_sec; ++ka)
    secs.push_back(minval + (double)ka*del);
      
  vector<shared_ptr<LRSplineSurface> > sfs;
  vol->constParamSurfaces(secs, dir, sfs);

  std::cout << "Constant parameter surfaces constructed" << std::endl;

  std::ofstream of("seg_sfs.g2");
  
  int ix[2];
  ix[0] = (dir+1) % 3;
  ix[1] = (dir+2) % 3;
  for (int ka=0; ka<nmb_sec; ++ka)
    {
      if (sfs[ka].get())
	{
	  sfs[ka]->writeStandardHeader(of);
	  sfs[ka]->write(of);
	}
      else
	{
	  std::cout << "No surface, ka= " << ka << std::endl;
	  continue;
	}
      
      const vector<CurveVec> curves = LRTraceIsocontours(*sfs[ka],
							 iso,
							 threshold_missing,
							 tol);

      for (size_t ki=0; ki<curves.size(); ++ki)
	{
	  for (size_t kj=0; kj<curves[ki].size(); ++kj)
	    {
	      int numc = curves[ki][kj].first->numCoefs();
	      std::vector<double>::const_iterator it = curves[ki][kj].first->coefs_begin();

	      vector<double> cf3(3*numc);
	      for (int kr=0; kr<numc; ++kr)
		{
		  int kh1, kh2;
		  for (kh1=0, kh2=0; kh1<3; ++kh1)
		    {
		      if (kh1 == dir)
			continue;
		    cf3[3*kr+ix[kh2]] = *it;
		    ++it;
		    ++kh2;
		    }
		  cf3[3*kr+dir] = secs[ka];
		}
	      
	      shared_ptr<SplineCurve> tmpcv(new SplineCurve(numc, curves[ki][kj].first->order(),
							    curves[ki][kj].first->knotsBegin(),
							    &cf3[0], 3));
	      
	      tmpcv->writeStandardHeader(ofs);
	      tmpcv->write(ofs);
	    }
	}
    }
}

		  
