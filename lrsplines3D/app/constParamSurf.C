#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"

using namespace std;
using namespace Go;


int main (int argc, char *argv[]) {

  if (argc != 5 && argc != 6) {
    cout << "usage: ./constParamSurf <input lrvol(.g2)> <output lrsurf (.g2)> <par. dir.> <par. val.> (<nmb sfs>)" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  int dir = atoi(argv[3]);
  double val = atof(argv[4]);
  int nmb = 1;
  if (argc == 6)
    nmb = atoi(argv[5]);

  GoTools::init();

  ObjectHeader oh;
  for (int ki=0; ki<nmb; ++ki)
    {
      oh.read(ifs);
      
      shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
      vol->read(ifs);

      Direction3D dir2 = (dir == 0) ? XDIR : ((dir == 1) ? YDIR : ZDIR);
      double t1 = vol->paramMin(dir2);
      double t2 = vol->paramMax(dir2);
      if (val >= t1 && val <= t2)
	{
	  shared_ptr<LRSplineSurface> lrsurf(vol->constParamSurface(val, dir));
      
	  lrsurf->writeStandardHeader(ofs);
	  lrsurf->write(ofs);
	}
    }
}
