#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;



int colors[3][3] = {
  {0, 255, 0},
  {255, 255, 255},
  {255, 0, 0},
};

int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: surface in (.g2), point cloud (.g2), points_out.g2, max level, nmb _levels" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  
  double max_level = atof(argv[4]);
  int nmb_level = atoi(argv[5]);
  double min_level = -max_level;

  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int ki, kj;
  double *curr;
  double dist;
  vector<double> limits(2*nmb_level+1);
  vector<vector<double> > level_points(2*nmb_level+2);

  // Set distance levels 
  double del = max_level/(double)nmb_level;
  limits[nmb_level] = 0;
  for (ki=1; ki<=nmb_level; ++ki)
    {
      limits[nmb_level-ki] = -ki*del;
      limits[nmb_level+ki] = ki*del;
    }

  // For each point, classify according to distance
  for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
    {
      // Evaluate
      Point pos;
      sf1->point(pos, curr[0], curr[1]);
      dist = curr[2]-pos[0];

      // Find classification
      for (kj=0; kj<limits.size(); ++kj)
	if (dist < limits[kj])
	  {
	    level_points[kj].push_back(curr[0]);
	    level_points[kj].push_back(curr[1]);
	    level_points[kj].push_back(curr[2]);
	    break;
	  }
      if (kj == limits.size())
	{
	  level_points[kj].push_back(curr[0]);
	  level_points[kj].push_back(curr[1]);
	  level_points[kj].push_back(curr[2]);
	}
    }

  // Write to file
  for (ki=0; ki<level_points.size(); ++ki)
    {
      if (level_points[ki].size() == 0)
	continue;

      // Make point cloud
      PointCloud3D level_cloud(level_points[ki].begin(), level_points[ki].size()/3);

      double cc[3];
      if (ki <= nmb_level)
	{
	  cc[0] = ((nmb_level-ki)*colors[0][0] + ki*colors[1][0])/nmb_level;
	  cc[1] = ((nmb_level-ki)*colors[0][1] + ki*colors[1][1])/nmb_level;
	  cc[2] = ((nmb_level-ki)*colors[0][2] + ki*colors[1][2])/nmb_level;
	}
      else
	{
	  cc[0] = ((ki-nmb_level-1)*colors[2][0] + 
		   (3*nmb_level-ki-1)*colors[1][0])/nmb_level;
	  cc[1] = ((ki-nmb_level-1)*colors[2][1] + 
		   (3*nmb_level-ki-1)*colors[1][1])/nmb_level;
	  cc[2] = ((ki-nmb_level-1)*colors[2][2] + 
		   (3*nmb_level-ki-1)*colors[1][2])/nmb_level;
	}

      fileout << "400 1 0 4 " << cc[0] << " " << cc[1];
      fileout << " " << cc[2] << " 255" << std::endl;
      level_cloud.write(fileout);
    }
}

      




 
