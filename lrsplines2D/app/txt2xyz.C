#include "GoTools/geometry/PointCloud.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 6) {
    std::cout << "Usage: infile outfile.g2 multx multy multz " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ofstream fileout(argv[2]); 

  double multx = atof(argv[3]);
  double multy = atof(argv[4]);
  double multz = atof(argv[5]);

  vector<double> pts;
  double av_xy = 0.0, av_z = 0.0;
  double max_xy = 0.0, max_z = 0.0;

  double xmin = HUGE, ymin = HUGE, zmin = HUGE;
  double xmax = 0.0, ymax = 0.0, zmax = 0.0;

  int nmb = 0;
  while (!filein.eof())
    {
      double x, y, z, eps_xy, eps_z;
      filein >> x >> y >> z >> eps_z >> eps_xy;
      pts.push_back(x*multx);
      pts.push_back(y*multy);
      pts.push_back(z*multz);

      xmin = std::min(xmin, x);
      ymin = std::min(ymin, y);
      zmin = std::min(zmin, z);
      xmax = std::max(xmax, x);
      ymax = std::max(ymax, y);
      zmax = std::max(zmax, z);

      av_xy += eps_xy;
      av_z += eps_z;
      max_xy = std::max(max_xy, eps_xy);
      max_z = std::max(max_z, eps_z);
      nmb++;
    }
  av_xy /= (double)nmb;
  av_z /= (double)nmb;

  std::cout << "Average xy: " << av_xy <<", max xy: " << max_xy << std::endl;
  std::cout << "Average z: " << av_z << ", max z: " << max_z << std::endl;

  std::cout << "Bound: " << xmin << " " << xmax << " " << ymin;
  std::cout << " " << ymax << " " << zmin << " " << zmax << std::endl;

  PointCloud3D points(pts.begin(), pts.size()/3);
  points.writeStandardHeader(fileout);
  points.write(fileout);
}
