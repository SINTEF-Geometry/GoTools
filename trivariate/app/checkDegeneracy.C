#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
  if (argc != 3 && argc != 5)
      cout << "Usage: " << "infile, tolerance, (swap pardir1, swap pardir2)" << endl;

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  ObjectHeader head;
  is >> head;

  // Read volume from file
  SplineVolume vol;
  is >> vol;

  double tol = atof(argv[2]);
  int dir1=-1, dir2=-1;
  if (argc > 3)
  {
      dir1 = atoi(argv[3]);
      dir2 = atoi(argv[4]);
  }

  if (dir1 >= 0 && dir2 >=0 && dir1 != dir2)
      vol.swapParameterDirection(dir1, dir2);

  int is_degen[6];
  vol.checkDegeneracy(tol, is_degen);
  for (int ki=0; ki<6; ++ki)
  {
      std::cout << "Degenerate: " << is_degen[ki] << std::endl;
      if (is_degen[ki])
      {
	  bool b, r, t, l;
	  int type;
	  bool deg;
	  deg = vol.isDegenerate(ki, type, b, r, t, l, tol);
	  std::cout << "b=" << b <<", r=" << r <<", t=" << t << ", l=" << l << std::endl;
      }
  }
}
