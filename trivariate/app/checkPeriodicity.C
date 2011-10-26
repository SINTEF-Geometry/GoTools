#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>

using namespace Go;
using namespace std;
using std::shared_ptr;

int main(int argc, char* argv[] )
{
  if (argc != 3 && argc != 5)
      cout << "Usage: " << "infile eps, (swap pardir1, swap pardir2)" << endl;

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  ObjectHeader head;
  is >> head;

  // Read volume from file
  SplineVolume vol;
  is >> vol;

  double eps = atof(argv[2]);

  int dir1=-1, dir2=-1;
  if (argc > 3)
  {
      dir1 = atoi(argv[3]);
      dir2 = atoi(argv[4]);
  }

  if (dir1 >= 0 && dir2 >=0 && dir1 != dir2)
      vol.swapParameterDirection(dir1, dir2);

  // For all parameter directions
  for (int ki=0; ki<3; ki++)
  {
      int period = vol.volumePeriodicity(ki, eps);
      std::cout << "Periodicity: " << period << std::endl;
  }
}


