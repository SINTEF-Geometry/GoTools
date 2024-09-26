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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;
using namespace Go;

int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};


int main (int argc, char *argv[]) {

  if (argc != 6) {
    cout << "usage: ./hightAccuracy <input 4d pt cloud(.g2)> <input lrspline(.g2)> <output hight> <output distance> <no levels>" << endl;
    return -1;
  }

  ifstream ifs1(argv[1]);
  ifstream ifs2(argv[2]);
  ofstream ofs1(argv[3]);
  ofstream ofs2(argv[4]);

  int level = atoi(argv[5]);

  int num_pts;
  ifs1 >> num_pts;

  GoTools::init();

  vector<double> pc4d;

  ObjectHeader oh;
  oh.read(ifs2);

  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifs2);

  double minh = 1.0e8, maxh = -1.0e8;
  double mind = 1.0e8, maxd = -1.0e8;
  double avd = 0.0;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs1 >> p0 >> p1 >> p2 >> q0;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
      minh = std::min(minh, q0);
      maxh = std::max(maxh, q0);

      Point pos;
      vol->point(pos, p0, p1, p2);
      double dist = q0 - pos[0];
      mind = std::min(mind, dist);
      maxd = std::max(maxd, dist);
      avd += fabs(dist)/(double)num_pts;
    }
  //avd /= (double)num_pts;
  std::cout << "Min height: " << minh << ", max height: " << maxh << std::endl;
  std::cout << "Max d: " << maxd << ", mind: " << mind << ", average: " << avd << std::endl;

  int ki;
  vector<double> limits_h(level+1);
  vector<vector<double> > level_h(level+2);

  // Set hight and distance levels 
  double delh = (maxh-minh)/(double)level;
  limits_h[0] = minh;
  for (ki=1; ki<=level; ++ki)
    {
      limits_h[ki] = ki*delh;
    }

  vector<double> limits_d(2*level+1);
  vector<vector<double> > level_d(2*level+2);

  double mleveld = std::max(fabs(mind), fabs(maxd));
  double deld = mleveld/(double)level;
  limits_d[level] = 0;
  for (ki=1; ki<=level; ++ki)
    {
      limits_d[level-ki] = -ki*deld;
      limits_d[level+ki] = ki*deld;
    }

  int kj;
  double *curr = &pc4d[0];
  for (int ix=0; ix!=num_pts; ++ix, curr+=4)
    {
      Point pos;
      vol->point(pos, curr[0], curr[1], curr[2]);
      double dist = curr[3] - pos[0];
 
     // Find classifications
      for (kj=0; kj<limits_d.size(); ++kj)
	if (dist < limits_d[kj])
	  {
	    level_d[kj].push_back(curr[0]);
	    level_d[kj].push_back(curr[1]);
	    level_d[kj].push_back(curr[2]);
	    break;
	  }
      if (kj == limits_d.size())
	{
	  level_d[kj].push_back(curr[0]);
	  level_d[kj].push_back(curr[1]);
	  level_d[kj].push_back(curr[2]);
	}

      for (kj=0; kj<limits_h.size(); ++kj)
	if (curr[3] < limits_h[kj])
	  {
	    level_h[kj].push_back(curr[0]);
	    level_h[kj].push_back(curr[1]);
	    level_h[kj].push_back(curr[2]);
	    break;
	  }
      if (kj == limits_h.size())
	{
	  level_h[kj].push_back(curr[0]);
	  level_h[kj].push_back(curr[1]);
	  level_h[kj].push_back(curr[2]);
	}

    }

  // Write to files
  for (ki=0; ki<level_d.size(); ++ki)
    {
      std::cout << "Level [";
      if (ki == 0)
	std::cout << "-inf,";
      else
	std::cout << limits_d[ki-1];
      if (ki < level_d.size()-1)
	std::cout << ", " << limits_d[ki] << "]: ";
      else
	std::cout << ", inf]: ";
      std::cout << level_d[ki].size()/3 << std::endl;
      if (level_d[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_d[ki].begin(), level_d[ki].size()/3);

      double cc[3];
      if (ki <= level)
	{
	  cc[0] = ((level-ki)*colors[0][0] + ki*colors[1][0])/level;
	  cc[1] = ((level-ki)*colors[0][1] + ki*colors[1][1])/level;
	  cc[2] = ((level-ki)*colors[0][2] + ki*colors[1][2])/level;
	}
      else
	{
	  cc[0] = ((ki-level-1)*colors[2][0] + 
		   (2*level-ki+1)*colors[1][0])/level;
	  cc[1] = ((ki-level-1)*colors[2][1] + 
		   (2*level-ki+1)*colors[1][1])/level;
	  cc[2] = ((ki-level-1)*colors[2][2] + 
		   (2*level-ki+1)*colors[1][2])/level;
	}

      ofs2 << "400 1 0 4 " << cc[0] << " " << cc[1];
      ofs2 << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(ofs2);
    }

  for (ki=0; ki<level_h.size(); ++ki)
    {
      if (level_h[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_h[ki].begin(), level_h[ki].size()/3);

      double cc[3];
      cc[0] = (ki*colors[2][0] + (level-ki)*colors[1][0])/level;
      cc[1] = (ki*colors[2][1] + (level-ki)*colors[1][1])/level;
      cc[2] = (ki*colors[2][2] + (level-ki)*colors[1][2])/level;

      ofs1 << "400 1 0 4 " << cc[0] << " " << cc[1];
      ofs1 << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(ofs1);
    }

}
