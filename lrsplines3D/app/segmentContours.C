#include <iostream>
#include <fstream>

#include "GoTools/geometry/GoTools.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"

using namespace std;
using namespace Go;


int main (int argc, char *argv[]) {

  if (argc != 9)
    {
      cout << "usage: ./segmentContours <input lrvol(.g2)> <output contours (.g2)> <par. dir.> <nmb sections> <min. par. val.> <max. par. val> <isoval> <tol>" << endl;
      return -1;
    }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  int dir = atoi(argv[3]);
  int nmb_sec = atoi(argv[4]);
  double minval = atof(argv[5]);
  double maxval = atof(argv[6]);
  double isoval = atof(argv[7]);
  double tol = atof(argv[8]);
  int threshold_missing = 100;

  GoTools::init();
  ObjectHeader oh;
  oh.read(ifs);
      
  shared_ptr<LRSplineVolume> vol(new LRSplineVolume());
  vol->read(ifs);

  std::cout << "LR Volume read" << std::endl;
  
  vector<double> iso;
  iso.push_back(isoval);
  double del = (nmb_sec == 1) ? 0 : (maxval - minval)/(double)(nmb_sec-1);
  vector<double> secs;
  for (int ka=0; ka<nmb_sec; ++ka)
    secs.push_back(minval + (double)ka*del);
      
  vector<shared_ptr<LRSplineSurface> > sfs;
  vol->constParamSurfaces(secs, dir, sfs);

  std::cout << "Constant parameter surfaces constructed" << std::endl;

  std::ofstream of("seg_sfs.g2");
  
  int ix[2];
  ix[0] = (dir+1) % 3;
  ix[1] = (dir+2) % 3;
  for (int ka=0; ka<nmb_sec; ++ka)
    {
      if (sfs[ka].get())
	{
	  sfs[ka]->writeStandardHeader(of);
	  sfs[ka]->write(of);
	}
      else
	{
	  std::cout << "No surface, ka= " << ka << std::endl;
	  continue;
	}
      
      const vector<CurveVec> curves = LRTraceIsocontours(*sfs[ka],
							 iso,
							 threshold_missing,
							 tol);

      for (size_t ki=0; ki<curves.size(); ++ki)
	{
	  for (size_t kj=0; kj<curves[ki].size(); ++kj)
	    {
	      int numc = curves[ki][kj].first->numCoefs();
	      std::vector<double>::const_iterator it = curves[ki][kj].first->coefs_begin();

	      vector<double> cf3(3*numc);
	      for (int kr=0; kr<numc; ++kr)
		{
		  int kh1, kh2;
		  for (kh1=0, kh2=0; kh1<3; ++kh1)
		    {
		      if (kh1 == dir)
			continue;
		    cf3[3*kr+ix[kh2]] = *it;
		    ++it;
		    ++kh2;
		    }
		  cf3[3*kr+dir] = secs[ka];
		}
	      
	      shared_ptr<SplineCurve> tmpcv(new SplineCurve(numc, curves[ki][kj].first->order(),
							    curves[ki][kj].first->knotsBegin(),
							    &cf3[0], 3));
	      
	      tmpcv->writeStandardHeader(ofs);
	      tmpcv->write(ofs);
	    }
	}
    }
}

		  
