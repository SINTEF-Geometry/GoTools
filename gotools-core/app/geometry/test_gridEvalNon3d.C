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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/Point.h"
#include <fstream>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{

  if (argc != 6)
    {
      cout << "Usage: " << argv[0] << " surfaceinfile surface3doutfile points3doutfile num_u num_v" << endl;
      exit(-1);
    }

  ifstream filein(argv[1]);
  ALWAYS_ERROR_IF(filein.bad(), "Bad or no curvee input filename");
  ObjectHeader head;
  filein >> head;
  if (head.classType() != SplineSurface::classType()) {
    THROW("Not a spline surface");
  }

  SplineSurface sf;
  filein >> sf;

  ofstream fileoutsurf(argv[2]);
  ALWAYS_ERROR_IF(fileoutsurf.bad(), "Bad surface output filename");

  ofstream fileoutpts(argv[3]);
  ALWAYS_ERROR_IF(fileoutpts.bad(), "Bad points output filename");

  int num_u = atoi(argv[4]);
  int num_v = atoi(argv[5]);

  vector<double> pts, param_u, param_v;

  sf.gridEvaluator(num_u, num_v, pts, param_u, param_v);

  vector<double> coefs3d;
  vector<Point> pts3d;
  int dim = sf.dimension();
  bool rational = sf.rational();

  int ctrl_pts = sf.numCoefs_u() * sf.numCoefs_v();
  vector<double>::const_iterator it = sf.ctrl_begin();

  for (int i = 0; i < ctrl_pts; ++i)
    {
      if (dim <= 3)
	for (int j = 0; j < 3; ++j)
	  {
	    if (j>=dim)
	      coefs3d.push_back(0.0);
	    else
	      {
		coefs3d.push_back(*it);
		++it;
	      }
	  }
      else
	{
	  for (int j = 0; j < 3; ++j, ++it)
	    coefs3d.push_back(*it);
	  it += (dim-3);
	}
      if (rational)
	{
	  coefs3d.push_back(*it);
	  ++it;
	}
    }

  int pts_pos = 0;
  for (int i = 0; i < num_u*num_v; ++i)
    {
      double x, y, z;
      if (dim == 0)
	x = 0.0;
      else
	x = pts[pts_pos];
      if (dim <= 1)
	y = 0.0;
      else
	y = pts[pts_pos+1];
      if (dim <= 2)
	z = 0.0;
      else
	z = pts[pts_pos+2];
      pts_pos += dim;
      pts3d.push_back(Point(x, y, z));
    }

  SplineSurface sf3d(sf.basis_u(), sf.basis_v(), coefs3d.begin(), 3, rational);

  sf3d.writeStandardHeader(fileoutsurf);
  sf3d.write(fileoutsurf);

  fileoutpts << "400 1 0 4 255 255 0 255" << endl;
  fileoutpts << pts3d.size() << endl;
  for (int i = 0; i < (int)pts3d.size(); ++i)
    fileoutpts << pts3d[i] << endl;
}

