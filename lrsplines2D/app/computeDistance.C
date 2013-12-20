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



int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: surface in (.g2) point cloud (.g2) lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream fileout(argv[3]); 
  
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

  int ki;
  double *curr;
  double dist;
  double maxdist = 0.0;
  double mindist = 0.0;
  double avdist = 0.0;
  for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
    {
      // Evaluate
      Point pos;
      sf1->point(pos, curr[0], curr[1]);
      dist = curr[2]-pos[0];

      fileout << curr[0] << "  " << curr[1] << "  " << curr[2];
      fileout << "  " << dist << std::endl;

      //dist = fabs(dist);
      maxdist = std::max(maxdist, dist);
      mindist = std::min(mindist, dist);
      avdist += fabs(dist);
    }

  avdist /= (double)nmb_pts;

  std::cout << "Max dist: " << maxdist << "Max dist below: " << mindist;
  std::cout << ", average dist: " << avdist << std::endl;
}


