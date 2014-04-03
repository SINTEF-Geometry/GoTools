#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineUtils.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
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
  if (argc != 4) {
    std::cout << "Usage: surface in (.g2), point cloud (.g2), param_points_out.raw" << std::endl;
    return -1;
  }

  std::ifstream sfin(argv[1]);
  std::ifstream ptsin(argv[2]);
  std::ofstream os(argv[3]); 
  
  ObjectHeader header1;
  header1.read(sfin);
  shared_ptr<LRSplineSurface> sf1(new LRSplineSurface());
  sf1->read(sfin);

  if (sf1->dimension() != 3)
    return -1; 

  ObjectHeader header2;
  header2.read(ptsin);
  PointCloud3D points;
  points.read(ptsin);

  BoundingBox box = sf1->boundingBox();
  Point low = box.low();
  Point high = box.high();
  Point mid = 0.5*(low + high);
   bool translate = true;
  if (translate)
    {
      Vector3D vec(-mid[0], -mid[1], -mid[2]);
      points.translate(vec);
      sf1->translate(-mid);
    }
  
  int nmb_pts = points.numPoints();
  vector<double> data(points.rawData(), points.rawData()+3*nmb_pts);

  int dim = sf1->dimension();
  int maxiter = 4;
  double aeps = 0.001;

  // double umin = sf1->paramMin(XFIXED);
  // double umax = sf1->paramMax(XFIXED);
  // double vmin = sf1->paramMin(YFIXED);
  // double vmax = sf1->paramMax(YFIXED);

  int ki, kj;
  double *curr;
  double dist;

  double maxdist = 0.0;
  double avdist = 0.0;

  std::streamsize prev = os.precision(15);
  os << nmb_pts << std::endl;

  // For each point, project onto surface
  for (ki=0, curr=&data[0]; ki<nmb_pts; ++ki, curr+=3)
    {
      // Get seed
      Point curr_pt(curr, curr+dim);
      LRBSpline2D *bspline = LRSplineUtils::mostComparableBspline(sf1.get(), curr_pt);
      Point seed = bspline->getGrevilleParameter();

      // Perform closest point
      double upar, vpar;
      Point close_pt;
      sf1->closestPoint(curr_pt, upar, vpar, close_pt,
			dist, aeps, maxiter, NULL, seed.begin());

      maxdist = std::max(maxdist, dist);
      avdist += fabs(dist);

      os << upar << " " << vpar << " ";
      if (translate)
	os << curr[0]+mid[0] << " " << curr[1]+mid[1] << " " << curr[2]+mid[2];
      else
	os << curr[0] << " " << curr[1] << " " << curr[2];
      os << std::endl;
    }

  avdist /= (double)nmb_pts;
  std::cout << "Max dist: " << maxdist << ", average dist: " << avdist << std::endl;
}

      




 
