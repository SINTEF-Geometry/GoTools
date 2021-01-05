#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/PointCloud.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#include <vector>

using std::vector;

int main (int argc, char *argv[]) {

  if (argc != 4) {
    std::cout << "usage: ./scale_intensity <input 4d pt cloud> <output point cloud> <multiplication factor> " << std::endl;
    return -1;
  }

  std::ifstream ifs(argv[1]);
  std::ofstream ofs(argv[2]);
  double fac = atof(argv[3]);

  int num_pts;
  ifs >> num_pts;

  vector<double> pc4d;

  double extent[6];
  extent[0] = extent[2] = extent[4] = 1.0e8; //std::numeric_limits<double>::max();
  extent[1] = extent[3] = extent[5] = -1.0e8; //std::numeric_limits<double>::lowest();
  double minv=1.0e8, maxv=-1.0e8;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs >> p0 >> p1 >> p2 >> q0;
      q0 *= fac;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
    }

  // Write to output
  ofs << num_pts << std::endl;
  for (int ix=0; ix!=num_pts; ++ix)
    {
      for (int ka=0; ka<4; ++ka)
	ofs << pc4d[4*ix+ka] << " ";
      ofs << std::endl;
    }
}
