#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRVolStitch.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;

int main (int argc, char *argv[]) {

  if (argc != 7) {
    cout << "usage: ./testMBA <input lrsplines(.g2)> <output lrsplines(.g2)> <nmb u> <nmb v> <nmb w> <cont>" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  int nmb_u = atoi(argv[3]);
  int nmb_v = atoi(argv[4]);
  int nmb_w = atoi(argv[5]);
  int cont = atoi(argv[6]);
  double eps = 0.001;
  
  int nmb = nmb_u*nmb_v*nmb_w;

  // Read spline volumes
  vector<shared_ptr<LRSplineVolume> > vols(nmb);
  ObjectHeader header;
  for (int ki=0; ki<nmb; ++ki)
    {
      header.read(ifs);
      vols[ki] = shared_ptr<LRSplineVolume>(new LRSplineVolume());
      vols[ki]->read(ifs);
    }
  
  LRVolStitch stitch;
  stitch.stitchRegVols(vols, nmb_u, nmb_v, nmb_w, eps, cont);
  
  for (size_t ki=0; ki<vols.size(); ++ ki)
    {
      vols[ki]->writeStandardHeader(ofs);
      vols[ki]->write(ofs);
    }
}
 
