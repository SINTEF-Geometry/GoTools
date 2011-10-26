#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>

using namespace Go;
using namespace std;
using std::shared_ptr;

int main(int argc, char* argv[] )
{
    if (argc != 2)
      cout << "Usage: " << "infile " << endl;

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  ObjectHeader head;
  is >> head;

  // Read volume from file
  SplineVolume vol;
  is >> vol;

  ofstream result("closest_points.g2");

  double eps = 1.0e-6;
  int more = 1;
  double xyz[3];
  double u_par, v_par, w_par;
  double dist;
  Point pnt, clo_pnt;
  while (more)
  {
      std::cout << "Give point: ";
      std::cin >>  xyz[0];
      std::cin >>  xyz[1];
      std::cin >>  xyz[2];

      pnt = Point(xyz[0], xyz[1], xyz[2]);

      vol.closestPoint(pnt, u_par, v_par, w_par, clo_pnt, dist, eps);
      std::cout << "Par: " << u_par << " " << v_par << " " << w_par;
      std::cout << ", point: " << clo_pnt[0] << " " << clo_pnt[1] << " " << clo_pnt[2];
      std::cout << ", dist: " << dist << std::endl;
      
      result << "400 1 0 4 0 255 0 255 \n";
      result << "1 \n";
      result << pnt[0]  << " " << pnt[1] << " " << pnt[2] << "\n";
      result << "400 1 0 4 255 0 0 255 \n";
      result << "1 \n";
      result << clo_pnt[0]  << " " << clo_pnt[1] << " " << clo_pnt[2] << "\n";

      std::cout << "More? ";
      std::cin >> more;
  }
}
