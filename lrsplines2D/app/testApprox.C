#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Usage: point cloud (.g2) lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]);
  
  ObjectHeader header;
  header.read(filein);
  PointCloud3D points;
  points.read(filein);

  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  const double AEPSGE = 1.0e-3;
  int dim = 1;
  LRSurfApprox approx(data, 1, AEPSGE, false, false);

  double maxdist, avdist; // will be set below
  int nmb_out_eps;        // will be set below
  int max_iter = 4;
  shared_ptr<LRSplineSurface> surf = approx.getApproxSurf(maxdist, avdist, nmb_out_eps, max_iter);

  std::cout << "Maxdist= " << maxdist << ", avdist= " << avdist;
  std::cout << ", nmb out= " << nmb_out_eps << std::endl;

  if (surf.get())
    {
      surf->to3D();
      surf->writeStandardHeader(fileout);
      surf->write(fileout);
    }
}

