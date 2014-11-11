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

#include "GoTools/lrsplines2D/LRTrimUtils.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 5) {
    std::cout << "Usage: point cloud (.g2), trim_out.g2, max recursion, div" << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  int max_rec = atoi(argv[3]);
  int nmb_div = atoi(argv[4]);

  ObjectHeader header;
  header.read(filein);
  PointCloud3D points;
  points.read(filein);

  BoundingBox box = points.boundingBox();
  Point low = box.low();
  Point high = box.high();
  Point mid = 0.5*(low + high);
  Vector3D vec(-mid[0], -mid[1], 0.0);
  points.translate(vec);

  std::ofstream of("translated_cloud.g2");
  points.writeStandardHeader(of);
  points.write(of);


  vector<double> points2(points.rawData(), 
			 points.rawData()+points.dimension()*points.numPoints());
  vector<vector<double> > seqs;
  LRTrimUtils trimutil(points2, 1);
  trimutil.computeTrimSeqs(max_rec, nmb_div, seqs);
  
  double udel, vdel;
  trimutil.getDomainLengths(udel, vdel);
  std::cout << "Minimum sub domains, diag = "<< trimutil.getDomainDiag();
  std::cout << ", udel = " << udel << ", vdel = " << vdel << std::endl;

  std::ofstream of2("translated_seq.g2");
  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      if (seqs[kr].size() <= 2)
	continue;
      of2 << "410 1 0 0" << std::endl;
      of2 << seqs[kr].size()/2-1 << std::endl;
      size_t ki, kj;
      for (ki=0; ki<seqs[kr].size()-4; ki+=2)
	{
	  for (kj=0; kj<2; ++kj)
	    of2 << seqs[kr][ki+kj] << "  ";
	  of2 << 0 << " ";
	  for (; kj<4; ++kj)
	    of2 << seqs[kr][ki+kj] << "  ";
	  of2 << 0 << std::endl;
	}
      of2 << seqs[kr][ki] << " " << seqs[kr][ki+1] << " " << 0 << " ";
      of2 << seqs[kr][0] << " " << seqs[kr][1] << " " << 0 << std::endl;
    }

  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      if (seqs[kr].size() <= 2)
	continue;
       fileout << "410 1 0 0" << std::endl;
      fileout << seqs[kr].size()/2-1 << std::endl;
      size_t ki, kj;
      for (ki=0; ki<seqs[kr].size()-4; ki+=2)
	{
	  for (kj=0; kj<2; ++kj)
	    fileout << seqs[kr][ki+kj]+mid[(int)kj] << "  ";
	  fileout << 0 << " ";
	  for (; kj<4; ++kj)
	    fileout << seqs[kr][ki+kj]+mid[(int)kj-2] << "  ";
	  fileout << 0 << std::endl;
	}
      fileout << seqs[kr][ki]+mid[0] << " " << seqs[kr][ki+1]+mid[1] << " " << 0 << " ";
      fileout << seqs[kr][0]+mid[0] << " " << seqs[kr][1]+mid[1] << " " << 0 << std::endl;
    }

  std::ofstream of3("extreme_params.g2");
  vector<double> umin_par, umax_par, vmin_par, vmax_par;
  double umin = 1.0e8;
  double umax = -1.0e8;
  double vmin = 1.0e8;
  double vmax = -1.0e8;
  double tol = 1.0e-6;
  for (size_t kr=0; kr<seqs.size(); ++kr)
    {
      size_t ki;
      for (ki=0; ki<seqs[kr].size(); ki+=2)
	{
	  if (fabs(seqs[kr][ki]-umin) < tol)
	    {
	      umin_par.insert(umin_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	  if (seqs[kr][ki] < umin)
	    {
	      umin = seqs[kr][ki];
	      umin_par.clear();
	      umin_par.insert(umin_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	  if (fabs(seqs[kr][ki]-umax) < tol)
	    {
	      umax_par.insert(umax_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	  if (seqs[kr][ki] > umax)
	    {
	      umax = seqs[kr][ki];
	      umax_par.clear();
	      umax_par.insert(umax_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	  if (fabs(seqs[kr][ki+1]-vmin) < tol)
	    {
	      vmin_par.insert(vmin_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	  if (seqs[kr][ki+1] < vmin)
	    {
	      vmin = seqs[kr][ki+1];
	      vmin_par.clear();
	      vmin_par.insert(vmin_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	  if (fabs(seqs[kr][ki+1]-vmax) < tol)
	    {
	      vmax_par.insert(vmax_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	  if (seqs[kr][ki+1] > vmax)
	    {
	      vmax = seqs[kr][ki+1];
	      vmax_par.clear();
	      vmax_par.insert(vmax_par.end(), seqs[kr].begin()+ki,
			      seqs[kr].begin()+ki+2);
	    }
	}
    }
  for (size_t ki=0; ki<umin_par.size(); ki+=2)
    {
      of3 << "400 1 0 0" << std::endl;
      of3 << "1" << std::endl;
      of3 << umin_par[ki] << " " << umin_par[ki+1] << " " << 0.0 << std::endl;
    }
  for (size_t ki=0; ki<vmin_par.size(); ki+=2)
    {
      of3 << "400 1 0 0" << std::endl;
      of3 << "1" << std::endl;
      of3 << vmin_par[ki] << " " << vmin_par[ki+1] << " " << 0.0 << std::endl;
    }
  for (size_t ki=0; ki<umax_par.size(); ki+=2)
    {
      of3 << "400 1 0 0" << std::endl;
      of3 << "1" << std::endl;
      of3 << umax_par[ki] << " " << umax_par[ki+1] << " " << 0.0 << std::endl;
    }
  for (size_t ki=0; ki<vmax_par.size(); ki+=2)
    {
      of3 << "400 1 0 0" << std::endl;
      of3 << "1" << std::endl;
      of3 << vmax_par[ki] << " " << vmax_par[ki+1] << " " << 0.0 << std::endl;
    }
}
