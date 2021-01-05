#include <iostream>
#include <fstream>

#include "GoTools/lrsplines3D/LRVolApprox.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;

int main (int argc, char *argv[]) {

  if (argc != 5 && argc != 6 && argc != 7) {
    cout << "usage: ./testMBA <input 4d pt cloud(.g2)> <output lrspline(.g2)> <tolerance> <levels> (<scaling factor q>) (<init MBA (0/1)>)" << endl;
    return -1;
  }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  double epsge = atof(argv[3]);
  int levels = atoi(argv[4]);
  int initMBA = 1;
  double facq = 1.0;
  if (argc >= 6)
    facq = atof(argv[5]);
  if (argc == 7)
    initMBA = atoi(argv[6]);

  // ObjectHeader oh;
  // oh.read(ifs);

  int num_pts;
  ifs >> num_pts;

  vector<double> pc4d;

  double domain[6];
  domain[0] = domain[2] = domain[4] = 1.0e8;
  domain[1] = domain[3] = domain[5] = -1.0e8;
  double mba_level = 0.0; //0.5*(minval+maxval); //0.0;
  double minval = std::numeric_limits<double>::max();
  double maxval = std::numeric_limits<double>::lowest();
  for (int ix=0; ix!=num_pts; ++ix)
    {
      double p0, p1, p2, q0;
      ifs >> p0 >> p1 >> p2 >> q0;
      q0 *= facq;
      pc4d.push_back(p0);
      pc4d.push_back(p1);
      pc4d.push_back(p2);
      pc4d.push_back(q0);
      
      domain[0] = std::min(domain[0], p0);
      domain[1] = std::max(domain[1], p0);
      domain[2] = std::min(domain[2], p1);
      domain[3] = std::max(domain[3], p1);
      domain[4] = std::min(domain[4], p2);
      domain[5] = std::max(domain[5], p2);

      mba_level += q0/(double)num_pts;
      minval = std::min(minval, q0);
      maxval = std::max(maxval, q0);
    }

  std::cout << "Domain: [" << domain[0] << "," << domain[1] << "]x[" << domain[2];
  std::cout << "," << domain[3] << "]x[" << domain[4] << "," << domain[5] << "]" << std::endl;
  std::cout << "Range: [" << minval << "," << maxval << "]" << std::endl;
  int dim = 1;
  int ncoef = 6;
  int order = 3;
  LRVolApprox vol_approx(ncoef, order, ncoef, order, ncoef, order,
			 pc4d, dim, domain, epsge, mba_level);
  vol_approx.setInitMBA(initMBA);
  
  double max, average, av_all;
  int num_out;
  cout << "Starting approximation..." << endl;

  shared_ptr<LRSplineVolume> result = vol_approx.getApproxVol(max,av_all,average,num_out,levels);

  result->writeStandardHeader(ofs);
  result->write(ofs);

  //result->setParameterDomain(0.0,1.0,0.0,1.0,0.0,1.0);

  //SplineVolume* tp_vol = result->asSplineVolume();

  // result->writeStandardHeader(ofs);
  // result->write(ofs);

  // tp_vol->writeStandardHeader(ofs);
  // tp_vol->write(ofs);

  cout << "It is done...\n"
       << "max = " << max << " "
       << "av_all = " << av_all << " "
       << "average = " << average << " "
       << "num_out ="  << num_out << " "
       << "levels = " << levels << " "
       << "num_elements = " << result->numElements() << " "
       << endl;

}
